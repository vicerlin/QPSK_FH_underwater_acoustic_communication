close all 
clear ;clc;
%   QPSK EPSK + HP(跳频) 调制 20230821gu
% 原始二进制数据
data = [0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1];
M = 8;
M_bit = log2(M);
msg_length = length(data);
Np = msg_length / M_bit;
fc = 10000;
T = 0.002;
fs = 400000;
t = 0 : 1/fs : T-1/fs;
bps = M_bit / T
%%跳频参数
hop_pattern = [8 3 1 6 9 5 2 7 4 10];  
lower_freq = 5000;                                          %换能器扫频范围
upper_freq = 15000;
Bandwith = upper_freq - lower_freq;
dhf = Bandwith/length(hop_pattern); 
fb = lower_freq - dhf; 
df = 1/T;
fk=hop_pattern.*dhf;
Fh0=[];
if Np <= length(hop_pattern)                                               %循环跳频图案表
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
constellation = exp(1j * (0:M-1) * 2 * pi/M);
% 星座点映射
modulatedSymbols = [];
for i = 1 : Np
    modulatedSymbols0 = constellation(symbols(i) + 1).*exp(1j*2*pi*Fh0(i)*t);
    modulatedSymbols0 = real(modulatedSymbols0); %发射信号取实部进行发射
    modulatedSymbols = [modulatedSymbols modulatedSymbols0];
end
figure;
plot(modulatedSymbols);



%% 解调

% 接收信号加入噪声
noisySymbols = modulatedSymbols; %+ 0.4 * (randn(size(modulatedSymbols)) + 1j * randn(size(modulatedSymbols)));

% 解调器 ，采用积分器解调
ReceiveData = [];
for i = 1 : Np
    sig = noisySymbols((i-1)*T*fs+1 : i*T*fs);
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


