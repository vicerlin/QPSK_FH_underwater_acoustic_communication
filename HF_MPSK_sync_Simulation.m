close all 
clear;clc;
%   QPSK EPSK + HP(跳频) 调制 20230821gu
% 原始二进制数据
% data = [0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1];
bits=[];
bitstring0 = textread('date.txt','%s');                 %480bit
bitstring0 = char(bitstring0(1));
for i = 1 : strlength(bitstring0)
    member = str2num(bitstring0(i));
    bits = [bits member];
end
data = [];
for i= 1 : 10
    data = [data bits];               %循环加长数据流
end

M = 4;
M_bit = log2(M);
msg_length = length(data);
Np = msg_length / M_bit;
fc = 10000;
lower_freq = 5000;                                                          %换能器扫频范围
upper_freq = 15000;
T = 0.002;
fs = 96000;
t = 0 : 1/fs : T-1/fs;
bps = M_bit / T
%%跳频参数
hop_pattern = [8 3 1 6 9 5 2 7 4 10];  
Bandwith = upper_freq - lower_freq;
dhf = Bandwith/length(hop_pattern); 
fb = lower_freq - dhf; 
df = 1/T;                                                                   %频率分辨率
fk=hop_pattern.*dhf;
Fh0=[];
if Np <= length(hop_pattern)                                                %循环跳频图案表
    for i = 1 : Np
        Fh0 = [Fh0 fk(i)+fb];
    end
else
    k = fix(Np/(length(hop_pattern)));
    for j = 1 : k
        for i = 1 : length(hop_pattern)
            Fh0 = [Fh0 fk(i)+fb];
        end
    end
    if mod(Np , length(hop_pattern)) > 0
        for i = 1 : mod(Np,length(hop_pattern))
            Fh0 = [Fh0 fk(i)+fb];
        end
    end
end

% 将二进制数据转换为符号（星座点）
symbols = bi2de(reshape(data, M_bit, [])', 'left-msb');
% 8PSK调制映射表
% constellation = 0.5*exp(1j *((0:M-1)*2*pi/M +pi/4));
constellation = exp(1j *(0:M-1)*2*pi/M);
% 星座点映射
modulatedSymbols = [];
for i = 1 : Np
    modulatedSymbols0 = constellation(symbols(i) + 1).*exp(1j*2*pi*Fh0(i)*t);
    modulatedSymbols0 = real(modulatedSymbols0);                            %发射信号取实部进行发射
    modulatedSymbols = [modulatedSymbols modulatedSymbols0];
end
% figure;
% plot(modulatedSymbols);

%% 添加同步头
symbol_sz=length(modulatedSymbols);%一帧信号长
coarse_interval=10*2048;%粗细同步间隔
sync_f1=lower_freq;
sync_f2=upper_freq;
sync_len=1700;
fine_interval=6*2048-(2048-sync_len);%细同步和信号间隔
head_begin=50*2048+3000;%开头空白
tail_end=100*2048;%尾巴空白
Nsync=0:sync_len-1;
ychirp=chirp(Nsync/fs,sync_f1,sync_len/fs,sync_f2);
ychrip_b=chirp(Nsync/fs,sync_f2,sync_len/fs,sync_f1);
sync_b=[ychrip_b, zeros(1,2048-sync_len)];
sync=[ychirp, zeros(1,2048-sync_len)];
sync=0.7*sync/max(abs(sync));
nsync=[];%循环同步头
for i=1:5
    nsync=[nsync sync];
end
% 信号拼接
QPSK_FH_singal_tx=[nsync,zeros(1,coarse_interval),sync,zeros(1,fine_interval),modulatedSymbols,zeros(1,fine_interval),sync_b,zeros(1,tail_end)];
% 输出音频文件
fname = 'test.wav';                          % 设定文件名称 注意格式 QPSK_FH_singal_1Kbps
audiowrite(fname,QPSK_FH_singal_tx,fs); 


%% 解调

    %对同步头和信号进行相关操作，截取信号
%     floder='20230620';
%     [y,yfs]=audioread(fullfile(floder,'QPSK_FH_singal_1Kbps.wav'));%2101_094648_1th_FH01.wav
    [y,yfs]=audioread('test.wav');%2101_094648_1th_FH01.wav
    fir_bwth = [lower_freq-500 upper_freq+500]*2/yfs;
    bpf = fir1(256,fir_bwth);
    rx_sig_bpf = conv(y,bpf);
    filtered_signal = rx_sig_bpf(129:end-128);
%     % 设计线性相位带通滤波器2 线性度佳
%     N = 50; % 滤波器阶数（奇数）
%     Frequencies_Hz_normalized= [0, (lower_freq-500)/fs*2, (upper_freq+500)/fs*2 ,1]; 
%     Gain_dB=[0,1,1,0];
%     filter_coefficients=firpm(N,Frequencies_Hz_normalized,Gain_dB);
%     % 应用滤波器
%     filtered_signal=filter(filter_coefficients, 1, y);
    sync_f = fliplr(sync);%将数组从左向右翻转
    corr1 = conv(filtered_signal,sync_f);
    [~,loc]= findpeaks(abs(corr1),'MinPeakDistance',25000,'MinPeakHeight',0.4);
    start_pos=loc(2)+fine_interval;
    my_sig = filtered_signal( start_pos+1 : start_pos+symbol_sz);
    my_sig =my_sig.';
    figure(2)
    plot(my_sig);spectrogram(my_sig,256,128,256,yfs,'yaxis')
    title('截取帧信号')  
    hold on;


%% 解调器 ，采用积分器解调
ReceiveData = [];
for i = 1 : Np
    sig = my_sig(round((i-1)*T*fs+1) : round(i*T*fs));
    data1 = sig .* exp(-1j*2*pi*Fh0(i)*t);
    %%积分
    ReceiveData0 = 2/T * trapz(t,data1);
    I(i) = real(ReceiveData0);
    Q(i) = imag(ReceiveData0);
    ReceiveData = [ReceiveData ReceiveData0]; 

end
scatterplot(ReceiveData);

%欧式距离来判断所在象限
demodulatedSymbols = zeros(size(msg_length));
for i = 1:Np
    [~, index] = min(abs(ReceiveData(i) - constellation));
    demodulatedSymbols(i) = index - 1;
end

% 解调后的二进制数据
decodedData = de2bi(demodulatedSymbols, M_bit, 'left-msb')';
decodedData = decodedData(:)';

%计算误码率
error_bin=find(decodedData-data);
error_rate1=length(error_bin)/length(data);
fprintf('原始误码率是 %f\n',error_rate1);


