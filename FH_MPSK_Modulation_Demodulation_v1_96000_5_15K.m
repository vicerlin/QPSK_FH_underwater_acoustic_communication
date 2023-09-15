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


%% 升余弦滤波器
ps_type = 'sqrt';
ps_grp_delay = 8;
beta = 0.25;
B=1000;%带宽
rate = Ffs*Tb;%这个rate是1变多的T个数
b= 2*rcosine(1, rate, ps_type, beta, ps_grp_delay);
a = 1;
N_bias = ps_grp_delay*rate;%过滤波器头尾多出来的点数
T=Tb;%码元周期

    %% 截取信号
    %对同步头和信号进行相关操作
    [y,yfs]=audioread('QPSK_FH_singal_1Kbps.wav');%2101_094648_1th_FH01.wav
    bpf = fir1(256,[4500,15500]*2/yfs);
    rx_sig_bpf = conv(y,bpf);
    rx_sig_bpf = rx_sig_bpf(129:end-128);
    rx_sig_bpf = rx_sig_bpf(501:end);%%为什么要再往后面取500个点,不对会同步头的寻找定位造成影响
    sync_f = fliplr(sync);%将数组从左向右翻转
    corr1 = conv(rx_sig_bpf,sync_f);
    [~,loc]= findpeaks(abs(corr1),'MinPeakDistance',25000,'MinPeakHeight',0.4);
    start_pos=loc(2)+fine_interval;
    my_sig = rx_sig_bpf( start_pos+1 : start_pos+symbol_sz);
    my_sig =my_sig.';
    figure(3)
    plot(my_sig);spectrogram(my_sig,256,128,256,yfs,'yaxis')
    title('截取帧信号')  
    hold on;


%% 解跳--同步截取码元长度后解扩
cut_receive=[];
filter_mix=[];
sig_despread=[];
I_r=[];Q_r=[];I_receive=[];Q_receive=[];
fpass=20000;  
%subplot(2,1,2)
for i=1:Nb/2
    cut_receive=my_sig(round(1+(i-1)*Tb*Ffs) : round(Tb*Ffs*i));
%     mix_cut_receive=cut_receive.*cos(2*pi*Fh(Tb*Ffs*i)*tt);
    I_cut_receive=cut_receive.*cos(2*pi*Fh(Tb*Ffs*i)*tt);      %I路乘以cos
    Q_cut_receive=-1*cut_receive.*sin(2*pi*Fh(Tb*Ffs*i)*tt);   %减去Q路乘以sin直接到基带
    I_receive=[I_receive I_cut_receive];
    Q_receive=[Q_receive Q_cut_receive];%混频，下变频

%     filter_mix_frequency=lowpass(mix_cut_receive,fpass,Ffs);                       %低通滤波器过滤高频分量
%     mix_frequency=abs(fftshift(fft(mix_cut_receive,length(tt))));                  %fft之后取模abs()，fft数据变换具有对称性,并fftshift（）将零频分量搬运到频谱中心                                
%     mix_frequency=mix_frequency(length(tt)/2+1:length(tt));
%     f=(length(tt)/2:length(tt)-1)*Ffs/length(tt)-Ffs/2;                           %将横坐标转化成频率显示
%     plot(f,mix_frequency); title('混频后频谱')                                    %Nyquist频率为fs/2=200HZ。整个频谱图是以Nyquist频率为对称轴的，只需考察0~Nyquist频率范围内的幅频特性
%     hold on; 
%     sig_despread=[sig_despread,filter_mix_frequency];
    %数据拼接得到解扩信号，下变频完成  
    %绘制
end
    I_baseband10=conv(b,I_receive);
    I_baseband=I_baseband10(1+N_bias:end-N_bias);
    Q_baseband10=conv(b,Q_receive);
    Q_baseband=Q_baseband10(1+N_bias:end-N_bias);
%% 判决
for i=1:Nb/2
    I_output=I_baseband((i-1)*rate+rate/2);
    I_r=[I_r I_output];
    if I_output>0 %抽样判决，大于0为1，否则为-1
        I_recover(i)=1;
    else
        I_recover(i)=-1;
    end
    Q_output=Q_baseband((i-1)*rate+rate/2);
    Q_r=[Q_r Q_output];
    if Q_output>0
        Q_recover(i)=1;
    else
        Q_recover(i)=-1;
    end
end
%% 星座图
message=[];
IQ_r=[];
for i=1:Nb/2
    IQ_r(i)=(I_r(i)+sqrt(-1)*Q_r(i));
    message=[message,IQ_r(i)];
end
scatterplot(message)
title('星座图')
%% 并串变换
demo_ori_bin=[];
for i=1:Nb
    if mod(i,2)~=0
        demo_ori_bin=[demo_ori_bin,I_recover((i-1)/2+1)];%奇数取I路信息
    else
        demo_ori_bin=[demo_ori_bin,Q_recover(i/2)];%偶数取Q路信息
    end
end
%%原始误码率计算
error_bin=find(demo_ori_bin-bitstream);
error_rate1=length(error_bin)/length(bitstream);
fprintf('原始误码率是 %f\n',error_rate1);

demo_ori_bin=(demo_ori_bin+1)/2;%双极性变单极性
string=bintostr(demo_ori_bin,length(demo_ori_bin));



