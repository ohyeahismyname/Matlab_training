clc
clear
close all

SNR_in_dB = 0:5:40;
SNR_weight = 45;
window=10;
N=2048;
frame_num=1;
QAM=16;
xxx=(0:QAM-1);
yyy=qammod(xxx,QAM);
Eavg=mean(abs(yyy).^2);

NF=1/sqrt(Eavg);
Tx_num=2;
Rx_num=2;

MSE_LS_dB=zeros(1,length(SNR_in_dB));
MSE_LMMS_dB=zeros(1,length(SNR_in_dB));
MSE_LI_dB=zeros(1,length(SNR_in_dB));

DMRS=0.7071+0.7071*1i ;

%% 自相關矩陣
R_H_HD=zeros(1644,822);
R_HD_HD=zeros(822,822);
DMRS_num=[204:2:1024,1027:2:1847];
real_num=[203:1:1024,1026:1:1847];

%R_HD_HD
for k=1:822
    for p=1:822
        if DMRS_num(k)==DMRS_num(p)
            R_HD_HD(p,k)=1;
        else
            %          if DMRS_num(k)~=real_num(p)
            R_HD_HD(p,k)=(1-exp(-1i*2*pi*window*(DMRS_num(p)-DMRS_num(k))/N))/(1i*2*pi*window*(DMRS_num(p)-DMRS_num(k))/N);
        end
    end
end

%R_H_HD
for k=1:1:822
    for p=1:1644
        if DMRS_num(k)==real_num(p)
            R_H_HD(p,k)=1;
        else
            R_H_HD(p,k)=(1-exp(-1*1i*2*pi*window*(real_num(p)-DMRS_num(k))/N))/(1i*2*pi*window*(real_num(p)-DMRS_num(k))/N);
        end
    end
end
SNR_Power=10^(SNR_weight/10);
%% main
for count=1:length(SNR_in_dB)
    MSE_LS=0;
    MSE_LMMSE=0;
    MSE_LI=0;

    for frame=1:frame_num
        fprintf("SNR:%d分之%d\t frame:%d分之%d\n",length(SNR_in_dB),count,frame_num,frame);
        SNR = 10^( SNR_in_dB(count)/10);
        No  = 10^(-SNR_in_dB(count)/10);
        n=sqrt(No  /  ( 2 ) )  *( randn (  Tx_num,1228800 )+1i*randn( Tx_num,1228800 ));
        decimal_data=randi([0,QAM-1],12,14,Tx_num); %產生資料
        binary_data=qammod(decimal_data,QAM,'gray')*NF;
        %插入DMRS
        binary_data(2:2:12,3,:)=DMRS;
        ex_decimal_data=repmat(binary_data,137,40,1);
        %CDM (-p)
        ex_decimal_data(2:4:1644,3:14:560,2) = -ex_decimal_data(2:4:1644,3:14:560,2);

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
        Total_H_Power= sum(10.^(ex_PowerdB./10)); %總能量為1
        Ntap=6;
        H_Channel= sqrt(10.^(ex_PowerdB./10))*sqrt(1/2);% 6*n*n
        H_Channel= H_Channel.*(sqrt(1/2)* (randn(Ntap,Tx_num,Rx_num)+1i*randn(Ntap,Tx_num,Rx_num)));%6*n*n

        %%
        y=zeros(Rx_num,1228805); %n*1228805
        for i=1:Tx_num           %1:6
            for j=1:Rx_num       %1:6
                y(j,:)=y(j,:)+conv(xhasCP(1,:,i) ,H_Channel  (:,j,i)  );
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
        H_frame=permute(H_frame,[1,4,2,3]);%1644*560*n*n

        %% LS
        Y_DMRS=YhasnoGB(2:2:1644,3:14:560,:);%接收到應該是DMRS的點 822 40
        H_LS=zeros(822,40,2,2);
        H_LS(1:2:822,:,1,1)= (1/2)*DMRS'*(Y_DMRS(1:2:822,:,1)+Y_DMRS(2:2:822,:,1));
        H_LS(1:2:822,:,1,2)= (1/2)*DMRS'*(Y_DMRS(2:2:822,:,1)-Y_DMRS(1:2:822,:,1));
        H_LS(1:2:822,:,2,1)= (1/2)*DMRS'*(Y_DMRS(1:2:822,:,2)+Y_DMRS(2:2:822,:,2));
        H_LS(1:2:822,:,2,2)= (1/2)*DMRS'*(Y_DMRS(2:2:822,:,2)-Y_DMRS(1:2:822,:,2));
        H_LS(2:2:822,:,:,:)=H_LS(1:2:822,:,:,:);%不知道為啥

        %% LMMSE
        H_LMMSE=zeros(1644,40,2,2);
        P0=eye(822);
        P1=eye(822);
        for i = 1:2:822
            P1(i,i) = -1;
        end

        H_LMMSE(:,:,1,1)= R_H_HD*P0'*inv(P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS'*Y_DMRS(:,:,1));
        H_LMMSE(:,:,2,1)= R_H_HD*P0'*inv(P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS'*Y_DMRS(:,:,2));
        H_LMMSE(:,:,1,2)= R_H_HD*P1'*inv(P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS'*Y_DMRS(:,:,1));
        H_LMMSE(:,:,2,2)= R_H_HD*P1'*inv(P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS'*Y_DMRS(:,:,2));

        %% 線性內差
        H_LI=zeros(1644,560,2,2);
        DMRS_ps=3:14:560;
        for symbol=1:560
            if symbol>=1 && symbol<=2  %前邊界
                foreward=11+symbol;
                backward=3-symbol;
                H_LI(:,symbol,:,:)=(foreward*H_LMMSE(:,1,:,:))/14+(backward*H_LMMSE(:,40,:,:))/14;
            end

            if symbol>=3 && symbol<=548 % 中間
                location=floor((symbol-3)/14)+1;
                foreward=symbol-DMRS_ps(location);
                backward=DMRS_ps(location+1)-symbol;
                H_LI(:,symbol,:,:)=(foreward*H_LMMSE(:,location+1,:,:))/14+backward*H_LMMSE(:,location,:,:)/14;
            end
            if symbol>=549 && symbol<=560 %後邊界
                foreward=symbol-549;
                backward=563-symbol;
                H_LI(:,symbol,:,:)=(foreward*H_LMMSE(:,1,:,:))/14+(backward*H_LMMSE(:,40,:,:))/14;
            end
        end


        %% 真實通道
        H_LS_R=H_frame(2:2:1644,3:14:560,:,:);
        H_LMMSE_R=H_frame(1:1:1644,3:14:560,:,:);

        %% 計算差
%         MSE_LS=MSE_LS+sum(abs(H_LS_R-H_LS).^2,'all');
        MSE_LMMSE=MSE_LMMSE+sum(abs(H_LMMSE_R-H_LMMSE).^2 ,'all');
%         MSE_LI=MSE_LI+sum(abs(H_frame-H_LI).^2,'all');
    end
    MSE_LS_dB(count)=10*log10(MSE_LS/(822*2*2*40*frame_num));
    MSE_LMMS_dB(count)=10*log10(MSE_LMMSE/(1644*40*2*2*frame_num));
    MSE_LI_dB(count)=10*log10(MSE_LI/(1644*560*2*2*frame_num));
end


plot(SNR_in_dB,MSE_LS_dB,'r-','LineWidth',2);
hold on;
plot(SNR_in_dB,MSE_LMMS_dB,'b-','LineWidth',2);
hold on;
plot(SNR_in_dB,MSE_LI_dB,'k-','LineWidth',2);
hold on;
grid on
title('5G-NR SISO-OFDM MSE of ZF');
xlabel('SNR (dB)');
ylabel('MSE');
legend('LS MSE','LMMSE MSE','LI MSE');