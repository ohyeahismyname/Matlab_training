clc
clear

SSB_5G_NR;
SNR_in_dB = 20 ; %dB
frame_num =10 ;
QAM=16;
xxx=(0:QAM-1);
yyy=qammod(xxx,QAM);
Eavg=mean(abs(yyy).^2);
NF=1/sqrt(Eavg);
q_bit=log2(QAM);
Tx_num=1;
Rx_num=1;

DMRS =-0.7071-0.7071*1i;
ans_matrix=zeros(1,frame_num);
% sync_offsets = zeros(1, frame_num);
%% SSB
PBCH=[zeros(702,1);qammod(randi([0,QAM-1],240,1),QAM,'gray');zeros(702,1)];% 把SSB存進去
PSS=[zeros(759,1);PSS;zeros(758,1)];
SSS=[PBCH(1:750);zeros(8,1);SSS;zeros(9,1);PBCH(895:end) ];
% SSS=[zeros(702,1);qammod(randi([0,QAM-1],48,1),QAM,'gray');zeros(8,1);SSS;zeros(9,1);qammod(randi([0,QAM-1],48,1),QAM,'gray');zeros(702,1)];


%%自相關
PSS_GB=[zeros(202,1);PSS(1:822);0;PSS(823:end);zeros(201,1)];
syn_freq=ifft(ifftshift(PSS_GB,1))*sqrt(2048) ;
%% main


for frame=1:frame_num
    SNR = 10^( SNR_in_dB/10);
    No  = 10^(-SNR_in_dB/10);
    n=sqrt(No/2)*(randn(1,1228800+200)+1i*randn(1,1228800+200));% 前後各加100個0當空白資料

    %訊號
    decimal_data=randi([0,QAM-1],12,14);
    ex_decimal_data=repmat(decimal_data,137,40);
    gray_moded = qammod(ex_decimal_data,QAM)*NF;

    %訊號加入SSB
    gray_moded(:,5)=PSS;
    gray_moded(:,6)=PBCH;
    gray_moded(:,7)=SSS;
    gray_moded(:,8)=PBCH;

    %加入GB
    GBhead=zeros(202,560);
    GBtail=zeros(201,560);
    DC=zeros(1,560);
    blank_space=zeros(1,100);

    TX=[GBhead;gray_moded(1:822,:);DC;gray_moded(823:1644,:);GBtail];
    X = ifft(ifftshift(TX,1))*sqrt(2048) ;

    %加入CP
    X_CP=zeros(1,1228800);
    start=1;
    for symbol=1:560
        if mod(symbol,28) ==0
            X_CP(1,start:start+208)=[X(2048-208+1:2048),symbol];
            X_CP(1,start+208:start+2048+208-1)=[X(1:2048,symbol)];
            start=start+2048+208;
        else
            X_CP(1,start:start+144)=[X(2048-144+1:2048),symbol];
            X_CP(1,start+144:start+2048+144-1)=[X(1:2048,symbol)];
            start=start+2048+144;
        end

    end

    %過通道
    PowerdB=[-2 -8 -10 -12 -15 -18];
    Total_H_Power= sum(10.^(PowerdB/10)); %總能量為1
    Ntap= 6;%通道數量
    H_Channel= sqrt(10.^(PowerdB/10));%db轉增益
%     H_Channel= H_Channel.*     (sqrt(1/(2*Tx_num))  * randn(1,Ntap)+1i*randn(1,Ntap));
    Y=conv(X_CP,H_Channel);
    after_channel_Y=Y(1:end-5);
    before_noise_Y=[blank_space,after_channel_Y,blank_space];
    before_synY=before_noise_Y+n;

    %     corr=xcorr(before_synY,syn_freq); %%相關係數
    corr=xcorr(syn_freq,before_synY);
    %     [corr_result, lags] = xcorr(syn_freq, before_synY,15);
    %     corr=xcorr(syn_freq,before_synY,15);
    [maxinum, startpositon]=max(corr);  %%找出訊號起始點
    [leng,non] = size(corr);            %%整個數列
    ans_matrix(frame)=ceil(leng/2)+1238112-100-startpositon-1247025;
end

%%畫圖
figure;
y_axis=zeros(1,31);
for i = 1:31
    syn=i-16 ;
    y_axis(i)=length(find(ans_matrix==syn));
end

bar(-15:15, y_axis);
title('Synchronization 不同的偏移量次數')
xlabel('偏移量')
ylabel('次數')

