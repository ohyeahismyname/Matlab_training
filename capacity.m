clc
clear
close all

SNR_in_dB = 0:5:40;
antenna=1:1:6;
QAM=16;
xxx=(0:QAM-1);
yyy=qammod(xxx,QAM);
Eavg=mean(abs(yyy).^2);
NF=1/sqrt(Eavg);
frame_num=100;
C_in_SNR= zeros(length(antenna),length(SNR_in_dB));

for count=1:length(antenna)
    Tx_num=(antenna(count));
    Rx_num=(antenna(count));

    for a=1:length(SNR_in_dB)
        Capacity=0;
        parfor frame=1:frame_num
		fprintf("MIMO:%d分之%d\t SNR:%d分之%d\t frame:%d分之%d\n",length(antenna),count,length(SNR_in_dB),a,frame_num,frame);
            SNR = 10^( SNR_in_dB(a)/10);
            No  = 10^(-SNR_in_dB(a)/10);
            n=sqrt(No   /   ( 2/Tx_num ) )  *( randn (  Tx_num,1228800 )+1i*randn( Tx_num,1228800    ));
            decimal_data=randi([0,QAM-1],1644,560,Tx_num); %產生資料
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

            %%通道
            PowerdB=[-2 -8 -10 -12 -15 -18].';
            ex_PowerdB=repmat(PowerdB,1,Tx_num,Rx_num);
            Total_H_Power= sum(10.^(ex_PowerdB/10)); %總能量為1
            Ntap=6;
            H_Channel= sqrt(10.^(ex_PowerdB/10));% 6*n*n
            H_Channel= H_Channel.*(sqrt(1/(2*Tx_num))* (randn(Ntap,Tx_num,Rx_num)+1i*randn(Ntap,Tx_num,Rx_num)));%6*n*n

            %%
            y=zeros(Rx_num,1228805); %n*1228805
            for i=1:Tx_num           %1:6
                for j=1:Rx_num       %1:6
                    y(j,:)=y(j,:)+conv(xhasCP(1,:,i),H_Channel  (:,j,i)  );
                end
            end
            newy=y(:,1:end-5);          %扣除後面五筆資料
            y_noise=newy+n;             %加上雜訊
            %%
            yhasnoCP= zeros(2048,560,Rx_num);%2048*560*n
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

            for Carrier= 1:1644
                for slot = 1:560
                    squareH=reshape(H_frame(Carrier,:,:,slot),Rx_num,Tx_num);
                    Capacity=Capacity+abs(log2 ( det (eye(Tx_num)+SNR*(1/Tx_num)*(squareH*squareH'))));
                end
            end

        end
        C_in_SNR(count,a)=Capacity/(1644*560*frame_num);
    end
end

%%畫圖

plot(SNR_in_dB,C_in_SNR(1,:),'r-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(2,:),'b-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(3,:),'g-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(4,:),'k-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(5,:),'m-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(6,:),'c-','LineWidth',2);
hold on;

grid on;
title('Capacity of MIMO OFDM');
xlabel('SNR (dB)');
ylabel('Capacity(bits/sec/Hz)');
legend('1x1','2x2','3x3','4x4','5x5','6x6');
