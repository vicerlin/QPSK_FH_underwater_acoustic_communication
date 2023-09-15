clear;close all;clc;
addpath('conv_coding_function\')
bits=[];
bitstring0=textread('date.txt','%s');                 %1440bit
bitstring0=char(bitstring0(1));
for i = 1:strlength(bitstring0)
    member=str2num(bitstring0(i));
    bits =[bits member];
end

data=[];
for i=1:10
    data =[data bits];               %循环加长数据流
end
bitstream=2*data-1;                              %将单极性转换成双极性
Nb=length(bitstream);                               %总数据量bit
T=0.002;                                            %比特周期2ms
fs=96000;                                           %采样频率96k

I=[];Q=[];                                                
 %数据流并转串，奇数进I路,偶数进Q路
for i=1:Nb
    if mod(i,2)~=0
        I=[I,bitstream(i)];                           %I、Q矩阵 Nb/2长度数组
    else
        Q=[Q,bitstream(i)];                           %每个码元的持续时间是输入码元的2倍
    end
end

%% QPSK编码
M = 4;
M_bit = log2(M);
Np = Nb / M_bit;
Tb = 0.002;                                                % Tb>=T  设置码元时间为2ms,传输速率为2/2=1Kbps
fc = 10000;                                               %调制的中频率为10kHZ
lower_freq = 5000;                                          %换能器扫频范围
upper_freq = 15000;
Ffs=fs;                                                 %采样频率96kHZ
tt=0:1/Ffs:Tb-1/Ffs;                                        %每个码元内的采样时间点
t_sig=0:1/Ffs:Np*Tb-1/Ffs;                                 %整段信号的采样时间点

I_carries=[];Q_carries=[];
for i=1:Nb/2
    I_carries=[I_carries,I(i)*cos(2*pi*fc*tt)];             %将数据流调制在在载波中
    Q_carries=[Q_carries,Q(i)*sin(2*pi*fc*tt)];
end
QPSK_singal=I_carries-Q_carries;                            %含N/2个元素(码元)


%% 跳频调制
hop_pattern=[8 3 1 6 9 5 2 7 4 10];  
fk=hop_pattern.*2/T;                                       %跳频的频偏,正交频带，混沌序列y=x^4+x
Fh=[];
fb=lower_freq-2/T;                                         %10k跳频基于的频率
if Nb/2<=length(hop_pattern)                               %循环跳频图案
    for i=1:Nb/2
        Fh=[Fh,fk(i)*ones(1,Tb*Ffs)+fb];
    end
else
    k=fix(Nb/(2*length(hop_pattern)));
    for j=1:k
        for i=1:length(hop_pattern)
            Fh=[Fh,fk(i)*ones(1,Tb*Ffs)+fb];
        end
    end
    if mod(Nb/2,length(hop_pattern))>0
        for i=1:mod(Nb/2,length(hop_pattern))
            Fh=[Fh,fk(i)*ones(1,Tb*Ffs)+fb];
        end
    end
end
disp(max(Fh))
disp(min(Fh))
I_FH=[];Q_FH=[];
window = blackman(length(tt))';
for i=1:Nb/2
%     disp(k);
%     disp(Fh(Tb*Ffs*i));
    I_FH=[I_FH,I(i)*cos(2*pi*Fh(Tb*Ffs*i)*tt).*window];
    Q_FH=[Q_FH,Q(i)*sin(2*pi*Fh(Tb*Ffs*i)*tt).*window];
end

QPSK_FH_singal=sqrt(2)/2*(I_FH-Q_FH);
QPSK_FH_singal=QPSK_FH_singal.*0.8;

figure(1)
subplot(2,1,1);
plot(t_sig,Fh);xlim([0,Nb/2*Tb]);legend("frequency hopping pattern");
subplot(2,1,2);
plot(t_sig,QPSK_FH_singal);legend("FH Signal with QPSK");xlabel('time /s','FontSize',8,'rotation',0,'Position',[Tb*Nb/2,-2.5,0]);




%% 添加同步头
symbol_sz=length(QPSK_FH_singal);%一帧信号长
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
%% 信号拼接
QPSK_FH_singal_tx=[nsync,zeros(1,coarse_interval),sync,zeros(1,fine_interval),QPSK_FH_singal,zeros(1,fine_interval),sync_b,zeros(1,tail_end)];
% 输出音频文件
fname = 'QPSK_FH_singal_1Kbps.wav';                          % 设定文件名称 注意格式 QPSK_FH_singal_1Kbps
%  audiowrite(fname,QPSK_FH_singal_tx,Ffs);                       % 输出文件


%% 积分解调
    %对同步头和信号进行相关操作
    [y,yfs]=audioread('QPSK_FH_singal_1Kbps.wav');%2101_094648_1th_FH01.wav
%     y = QPSK_FH_singal_tx;
%     yfs = fs;
    fir_bwth = [lower_freq-500 upper_freq+500]*2/yfs;
    bpf = fir1(256,fir_bwth);
    rx_sig_bpf = conv(y,bpf);
    filtered_signal = rx_sig_bpf(129:end-128);
    sync_f = fliplr(sync);%将数组从左向右翻转
    corr1 = conv(filtered_signal,sync_f);
    [~,loc]= findpeaks(abs(corr1),'MinPeakDistance',25000,'MinPeakHeight',0.4);
    start_pos=loc(2)+fine_interval;
    my_sig = filtered_signal( start_pos+1 : start_pos+symbol_sz);
    my_sig =my_sig.'; %如果129行注释去掉即y = QPSK_FH_singal_tx; 此句需要注释
    figure(2)
    plot(my_sig);spectrogram(my_sig,256,128,256,yfs,'yaxis')
    title('截取帧信号')  
    hold on;

%% 解调器 ，采用积分器解调
Fh0=[]; t=tt; msg_length = Nb;
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
ReceiveData = [];
for i = 1 : Np
    sig = my_sig(round((i-1)*T*fs+1) : round(i*T*fs));
    data1 = sqrt(2).*sig .* sin(2*pi*Fh0(i)*t);
    data2 = -1*sqrt(2).*sig .* cos(2*pi*Fh0(i)*t);
    %%积分
    a = 4/T * trapz(t,data1);
    b = 4/T * trapz(t,data2);
    ReceiveData0 = a + sqrt(-1)*b;
    ReceiveData = [ReceiveData ReceiveData0]; 

end
constellation = exp(1j *([0 1 3 2]*2*pi/M+pi/4));

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
