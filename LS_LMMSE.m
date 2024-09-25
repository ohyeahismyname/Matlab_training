clc
clear
close all

SNR_in_dB = 0:5:30;
SNR_weight = 40;
window=6;
frame_num=50;
QAM=16;
xxx=(0:QAM-1);
yyy=qammod(xxx,QAM);
Eavg=mean(abs(yyy).^2);
NF=1/sqrt(Eavg);
q_bit=log2(QAM);
Tx_num=1;
Rx_num=1;

MSE_LS_dB=zeros(1,length(SNR_in_dB));
MSE_LMMS_dB=zeros(1,length(SNR_in_dB));

DMRS=-0.7071-0.7071*1i ;

W=weight_matrix(window,SNR_weight);

for count=1:7
    MSE_LS=0;
    MSE_LMMSE=0;

    for frame=1:frame_num
        SNR=10^( SNR_in_dB(count)/10);
        No=10^(-SNR_in_dB(count)/10);
        n=sqrt(No/2)*(randn(1,1228800)+1i*randn(1,1228800));
        	

        decimal_data=randi([0,QAM-1],12,14);
        gray_moded = qammod(decimal_data,QAM)*NF;

        %%加入DMRS
        gray_moded(2:2:12,3)=DMRS;
        ex_gray_mode=repmat(gray_moded,137,40);

        %加入GB
        GBhead=zeros(202,560);
        GBtail=zeros(201,560);
        DC=zeros(1,560);

        TX=[GBhead;ex_gray_mode(1:822,:);DC;ex_gray_mode(823:1644,:);GBtail];
        X = ifft(ifftshift(TX,1))*sqrt(2048) ;

        %加入CP
        X_CP=zeros(1,1228800);
        start=1;
        
%         for symbol=1:560
%             if mod(symbol,28) ==0
%                 X_CP(1,start:start+208)=[X(2048-208+1:2048),symbol];
%                 X_CP(1,start+208:start+2048+208-1)=[X(1:2048,symbol)];
%                 start=start+2048+208;
%             else
%                 X_CP(1,start:start+144)=[X(2048-144+1:2048),symbol];
%                 X_CP(1,start+144:start+2048+144-1)=[X(1:2048,symbol)];
%                 start=start+2048+144;
%             end
% 
%         end

        %過通道
        PowerdB=[-2 -8 -10 -12 -15 -18];
        Ntap= 6;%通道數量
        Total_H_Power= sum(10.^(PowerdB/10)); %總能量為1
        
        H_Channel= sqrt(10.^(PowerdB/10));%db轉增益
        H_Channel= H_Channel.*(sqrt(1/(2*Tx_num))* randn(1,Ntap)+1i*randn(1,Ntap));
        y=conv(X_CP,H_Channel);
        newy=y(1:end-5);
        noise_y=newy+n;

        %remove CP
        yhasnoCP=zeros(2048,560);
        start2=1;
        for symbol=1:560
            if mod(symbol,28) ==0
                yhasnoCP(:,symbol)=noise_y(1,start2+208:start2+2048+208-1);
                start2=start2+2048+208;
            else
                yhasnoCP(:,symbol)=noise_y(1,start2+144:start2+2048+144-1);
                start2=start2+2048+144;
            end
        end

        %remove GB
        Y = fftshift(fft(yhasnoCP/sqrt(2048) ) ,1);
        YhasnoGB=[Y(203:1024,:); Y(1026:1847,:)];



        %估測結果
        Y_LS=YhasnoGB(2:2:1644,3:14:560);%接收到應該是DMRS的點 822 40
        H_LS=Y_LS/DMRS;

        %         X_LS=DMRS*eye(822);%放入的每個點 822 822
        %         H_LS=inv(X_LS)*Y_LS; %822 40

        h=[H_Channel,zeros(1,2042)];
        H=fftshift(fft(h));
        H_Data=[H(1,203:1024),H(1,1026:1847)].';
        H_frame=repmat(H_Data,1,560);

        %         H=[H_Channel,zeros(1,34)];
        %         R_H_H=xcorr2(H,H');
        %         R_H_HS=xcorr2(H,H_LS');
        %         R_HS_HS=xcorr2(H_LS,H_LS);
        %
        %         H=[H_Channel,zeros(1,2042)];
        %         R_H_H=xcorr2(H,H');
        %         W=R_H_HS*inv(R_HS_HS);


        H_LMMSE=W * H_LS;
        H_LS_R=H_frame(2:2:1644,3:14:560);
        H_LMMSE_R=H_frame(1:1:1644,3:14:560);
        MSE_LS=MSE_LS+sum(abs(H_LS_R-H_LS),'all');
        MSE_LMMSE=MSE_LMMSE+sum(abs(H_LMMSE_R-H_LMMSE),'all');
        
    end
    MSE_LS_dB(count)=10*log(MSE_LS/(1644000));
    MSE_LMMS_dB(count)=10*log(MSE_LMMSE/(3288000));
end




plot(SNR_in_dB,MSE_LS_dB,'r-','LineWidth',2);
hold on;
plot(SNR_in_dB,MSE_LMMS_dB,'b-','LineWidth',2);
hold on;
grid on
title('5G-NR SISO-OFDM MSE of ZF');
xlabel('SNR (dB)');
ylabel('SNR(dB) from MSE');
legend('LS MSE','LMMSE MSE');