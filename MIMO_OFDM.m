clc
clear
close all

%% 設定區
frame_num=5;
SNR_in_dB = 0:5:40;
BER_SNR=zeros(1,length(SNR_in_dB));
BER_SNR2=zeros(1,length(SNR_in_dB));

%% 不變設定區

Tx_num=4;
Rx_num=4;
QAM=16;
q_bit=log2(16);
xxx=(0:QAM-1);
yyy=qammod(xxx,QAM);
Eavg=mean(abs(yyy).^2);
NF=1/sqrt(Eavg);

%%

for a=1:length(SNR_in_dB)
    BER1=0;
    BER2=0;
    SNR = 10^( SNR_in_dB(a)/10);
    No  = 10^(-SNR_in_dB(a)/10);
    for frame=1:frame_num
        fprintf("SNR:%d分之%d \t frame:%d分之%d\n",length(SNR_in_dB),a,frame_num,frame);
        n=sqrt(No / (2) ) *( randn (  Tx_num,1228800 )   +1i*randn( Tx_num,1228800 )   );

        %% 產生資料
        decimal_data=randi([0,QAM-1],1644,560,Tx_num);  %產生資料
        binary_data=dec2bin(decimal_data,q_bit);        %轉二進制
        ex_decimal_data=qammod(decimal_data,QAM,'gray')*NF; %normalized

        GBhead=zeros(202,560,Tx_num);
        GBtail=zeros(201,560,Tx_num);
        DC=zeros(1,560,Tx_num);
        TX=[GBhead;ex_decimal_data(1:822,:,:);DC;ex_decimal_data(823:1644,:,:);GBtail];%2048*560*n
        X = ifft(ifftshift(TX,1))*sqrt(2048) ;%2048*560*n

        xhasCP=zeros(1,1228800,Tx_num);     %空白 1*1228800*n
        start1=1;
        for symbol=1:560
            if mod(symbol,28)-1
                xhasCP(1, start1:start1+2048+144-1,:)=[X(2048-144+1:2048,symbol,:);X(:,symbol,:)];
                start1 = start1+2048+144;
            else
                xhasCP(1, start1:start1+2048+208-1,:)=[X(2048-208+1:2048,symbol,:);X(:,symbol,:)];
                start1 = start1+2048+208;
            end
        end

        %% 通道
        PowerdB=[-2 -8 -10 -12 -15 -18].';
        ex_PowerdB=repmat(PowerdB,1,Tx_num,Rx_num);
        Total_H_Power= sum(10.^(ex_PowerdB/10)); %總能量為1
        Ntap=6;
        H_Channel= sqrt(10.^(ex_PowerdB/10));% 6*n*n
        H_Channel= H_Channel.*(sqrt(1/(2*Tx_num))* (randn(Ntap,Tx_num,Rx_num)+1i*randn(Ntap,Tx_num,Rx_num)));

        %%
        y=zeros(Rx_num,1228805); %n*1228805
        for i=1:Tx_num           %1:4
            for j=1:Rx_num       %1:4
                y(j,:)=y(j,:)+conv(xhasCP(1,:,i)   ,H_Channel  (:,j,i)  );
            end
        end
        newy=y(:,1:end-5);          %扣除後面五筆資料
        y_noise=newy+n;             %加上雜訊
        %%
        yhasnoCP= zeros(2048,560,Rx_num);%空白2048*560*n
        start2=1;
        for symbol = 1:560
            if mod(symbol,28)-1
                for i=1:Rx_num
                    yhasnoCP(:,symbol,i) = y_noise(i,start2+144:start2+144+2048-1);
                end
                start2=start2+144+2048;
            else
                for i=1:Rx_num
                    yhasnoCP(:,symbol,i) = y_noise(i,start2+208:start2+208+2048-1);
                end
                start2=start2+208+2048;
            end
        end
        Y=fftshift(fft(yhasnoCP/sqrt(2048) ) ,1);%2048*560*n
        YhasnoGB=[Y(203:1024,:,:); Y(1026:1847,:,:)];%1644*560*n

        h=[H_Channel;zeros(2042,Tx_num,Rx_num)];%2048*n*n
        H=fftshift(fft(h,[],1),1);
        H_Data=[H(203:1024,:,:);H(1026:1847,:,:)];%1644*n*n

        H_frame=repmat(H_Data(:,:,:),1,1,1,560);%1644*n*n*560
        H_frame=permute(H_frame,[1,4,2,3]); %置換維度 1644*560*n*N

        %         x_hat2=zeros(1644,560,4);

        %通道估計
        x_hat=zeros(1644,560,4);            %空白估計值
        x_hat2=zeros(1644,560,4);

        time =zeros(1,1644*560);
        for carrier=1:1644
            for slot =1:560
                zf=zeros(Tx_num,1);
                h_temp=squeeze(H_frame(carrier,slot,:,:));  %變成二維
                y_temp=squeeze(YhasnoGB(carrier,slot,:));      %變成一維
                zf=zf+(inv(h_temp'*h_temp)*h_temp')*y_temp;      %inv(H'H)H'y
                x_hat(carrier,slot,:)=zf;

            end
        end
        x_hat2=ZF3(H_frame,YhasnoGB,Tx_num);

        x_hat=x_hat/NF;
        decimal_data_hat = qamdemod(x_hat,QAM,'gray');
        binary_data_hat=dec2bin(decimal_data_hat,q_bit);
        BER1=BER1+sum(sum(binary_data~= binary_data_hat));

        x_hat2=x_hat2/NF;
        decimal_data_hat2 = qamdemod(x_hat2,QAM,'gray');
        binary_data_hat2=dec2bin(decimal_data_hat2,q_bit);
        BER2=BER2+sum(sum(binary_data~= binary_data_hat2));

    end
    BER_SNR(1,a) =BER1/(1644*560*q_bit*frame_num*Tx_num);
    BER_SNR2(1,a) =BER2/(1644*560*q_bit*frame_num*Tx_num);
end


%% 畫圖
semilogy(SNR_in_dB,BER_SNR(1,:),'r-X','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR2(1,:),'b--O','LineWidth',2)
hold on
grid on;
title('BER of 5G NR MIMO OFDM');
xlabel('SNR (dB)');
ylabel('BER');
legend('ZF in matlab','ZF in c');
