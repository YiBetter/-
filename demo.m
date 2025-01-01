clc;
clear;
close all;

%% 用于语音信号传输的无线通信系统
% 读取WAV文件
[file, path] = uigetfile('*.wav', '选择一个音频文件');
if isequal(file, 0)
    disp('用户未选择文件');
    return;
end
audioFile = fullfile(path, file);

tic;
% 读取音频信号
[signal, fs] = audioread(audioFile);
% 取单通道音频信号，如果是立体声，选择一个通道
if size(signal, 2) > 1
    signal = signal(:, 1);  % 选择左声道
end

% 定义量化位数
quantizationBits = 6;  % 6位量化
levels = 2^quantizationBits;  % 量化级数

% 对音频信号进行量化处理
signalMin = min(signal);
signalMax = max(signal);
signalRange = signalMax - signalMin;

% 归一化音频信号
normalizedSignal = (signal - signalMin) / signalRange;  

% 量化音频信号（将信号离散化到6位量化级别）
quantizedSignal = round(normalizedSignal * (levels - 1));

% 信源编码（长度为6）
sourceEncodedData = signalToSourceCode(quantizedSignal, 6);

% 信道编码（(7,3)汉明码）
channelEncodedData = hammingEncode(sourceEncodedData);

% 调制（使用BPSK调制）
carrierFreq = 2e9; % 载波频率2GHz
modulatedSignal = bpskModulate(channelEncodedData, carrierFreq, fs);

% 无线信道（加噪 + 衰落）
SNR = 40;
p = 0.01;
% 信道模拟
receivedSignal = simulateChannel(modulatedSignal, SNR, p);

% 解调
demodulatedSignal = bpskDemodulate(receivedSignal, carrierFreq, fs);

% 信道译码（汉明码解码）
[channelDecodedData,twoBitsError] = hammingDecode(demodulatedSignal);
if twoBitsError==true
    disp('传输过程中在某一个分组内(7,3)分组码发生了2比特错误');
end

% 信源译码（将二进制数组转化为6位带符号二进制数）
receivedSignal = sourceDecode(channelDecodedData, 6);

% 恢复原始量化信号，将数字信号转换回模拟(音频)信号
receivedSignal = receivedSignal/(levels-1);
recoveredSignal = receivedSignal*signalRange + signalMin;

% 绘制结果
figure;
subplot(2, 1, 1);
plot(signal);
title('原始音频信号');
subplot(2, 1, 2);
plot(recoveredSignal);
title('恢复后的音频信号(SNR=40dB, p=0.01)');

% 播放恢复后的音频
sound(recoveredSignal, fs);
toc;
%% --- 辅助函数 ---

% 信源编码：将信号按指定长度编码（例如6位）
function encodedData = signalToSourceCode(data, codeLength)
    encodedData = zeros(1,length(data)*6);
    for i = 1:length(data)
        encodedData(i*6-5:i*6) = dec2bin(data(i), codeLength) - '0';
    end
end

% 汉明编码（7,3）：将3位数据编码为7位
function encodedData = hammingEncode(data)
    encodedData = zeros(1,length(data)*7/3);
    for i = 1:length(data)/3
        block = data((i-1)*3+1:i*3); % 取3位数据
        % (7,3)汉明码编码
        P1=mod(block(1) + block(2),2);
        P2=mod(block(1) + block(3),2);
        P3=mod(block(2) + block(3),2);
        P4=mod(P1+P2+P3+block(1)+block(2)+block(3),2);
        encodedBlock = [P1,P2,block(1),P3,block(2),block(3),P4];
        encodedData(i*7-6:i*7) = encodedBlock;
    end
end

% 汉明解码（7,3）：解码并检测错误
function [decodedData,flag] = hammingDecode(data)
    flag = false;
    decodedData = zeros(1,length(data)*3/7);
    for i = 1:length(data)/7
        block = data((i-1)*7+1:i*7);
        % 校验并解码
        S1 = mod(block(1) + block(3) + block(5),2);
        S2 = mod(block(2) + block(3) + block(6),2);
        S3 = mod(block(4) + block(5) + block(6),2);
        S4 = mod(block(1)+block(2)+block(3)+block(4)+block(5)+block(6)+block(7),2);
        if S4==0 % 发现两位错或无错
            if S1+S2+S3~=0 % 发现两位错向上报告，无法纠正
                flag = true;
            end
        else % 发现一位错并纠正
            errorPos = S1+S2*2+S3*4;
            if errorPos==0
                block(7) = mod(block(7)+1,2);
            else
                block(errorPos) = mod(block(errorPos) + 1, 2); % 修正错误;
            end
            %disp(['第',num2str(i),'块出错']);
        end
        decodedData(i*3-2:i*3) = [block(3),block(5),block(6)];
    end
end

% BPSK调制
function modulatedSignal = bpskModulate(data, carrierFreq, fs)
    t = (0:length(data)-1) / fs;
    carrier = cos(2 * pi * carrierFreq * t);
    
    % 将0映射为-1，1映射为+1
    dataBPSK = 2 * data - 1;  % 将0和1映射到-1和+1

    modulatedSignal = carrier .* dataBPSK;
end

% --- 信道模拟函数 ---
function receivedSignal = simulateChannel(signal, SNR, fadingProb)
    % 模拟信道传输的函数，包括Rayleigh衰落和AWGN噪声

    % 计算信号的平均功率
    signalPower = mean(abs(signal).^2); % 平均功率计算公式
    
    % Rayleigh 衰落方差
    fadingVar = signalPower / 2; % 设置 Rayleigh 衰落的方差为信号平均功率的一半
    
    % Rayleigh 衰落
    fadingMask = rand(size(signal)) > fadingProb; % 深度衰落概率为 fadingProb
    rayleighFading = sqrt(fadingVar) * abs(randn(size(signal)) + 1j * randn(size(signal))) / sqrt(2);
    
    % 施加衰落
    fadedSignal = signal .* rayleighFading .* fadingMask;

    % 添加AWGN噪声
    receivedSignal = awgn(fadedSignal, SNR, 'measured');
end

% BPSK解调
function demodulatedSignal = bpskDemodulate(receivedSignal, carrierFreq, fs)
    t = (0:length(receivedSignal)-1) / fs;
    carrier = cos(2 * pi * carrierFreq * t);
    
    % 乘以载波并低通滤波
    demodulatedSignal = receivedSignal .* carrier;
    
    % 将低通滤波的截止频率设置为 Nyquist 范围内
    cutoffFreq = min(fs / 2 - 1, carrierFreq / 2); % 确保小于奈奎斯特频率
    demodulatedSignal = lowpass(demodulatedSignal, cutoffFreq, fs);
    
    % 判决为0或1
    demodulatedSignal(demodulatedSignal > 0) = 1;
    demodulatedSignal(demodulatedSignal <= 0) = 0;
end

% 信源解码：恢复原始信号
function recoveredSignal = sourceDecode(data, codeLength)
    recoveredSignal = zeros(1,length(data)/codeLength);
    for i = 1:length(data)/codeLength
        block = data((i-1)*codeLength+1:i*codeLength);
        recoveredSignal(i) = bin2dec(num2str(block));
    end
end
