
clear
clc
close all

frame_num = 8;
SNR_in_dB = 0:5:40;        % 自己設訂雜訊大小
SNR_weight= 45;
window= 10;
DMRS_value=0.7071 + 0.7071*1i ;
Fs= 122.88*10^6;
delta_f= 1000;

%% 不變設定區

Tx_num=2;
Rx_num=2;
QAM=16;
q_bit=log2(16);
xxx=(0:QAM-1);
yyy=qammod(xxx,QAM);
Eavg=mean(abs(yyy).^2);
NF=1/sqrt(Eavg);

N=2048;
R_HD_HD	=zeros( 822,822);
R_H_HD	=zeros(1644,822);
DMRS_pos=[204:2:1024,1027:2:1847];
Real_pos=[203:1:1024,1026:1:1847];
for p=1:822
    for k=1:822%R_HD_HD 822*822
        if DMRS_pos(k)==DMRS_pos(p)
            R_HD_HD(k,p)=1;
        else
            R_HD_HD(k,p)=( 1-exp(-1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N))/( 1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N );
        end
    end
    for k=1:1644
        if Real_pos(k)==DMRS_pos(p)
            R_H_HD(k,p)=1;
        else
            R_H_HD(k,p)=(1-exp(-1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N))/( 1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N );
        end
    end
end
SNR_Power=10^(SNR_weight/10);

BER_SNR=zeros(1,length(SNR_in_dB));
BER_SNR2=zeros(1,length(SNR_in_dB));

for count=1:length(SNR_in_dB)
    SNR = 10^( SNR_in_dB(count)/10);
    No  = 10^(-SNR_in_dB(count)/10);
    BER = 0;                         
    BER2 = 0;
    fprintf("SNR : %d分之%d \n\n",length(SNR_in_dB),count);
    parfor frame=1:frame_num
        %輸入資料(含DMRS)
        decimal_data=randi([0,QAM-1],Tx_num*(1644*560-822*40),1);
        data_mod_L=qammod (decimal_data,QAM,'gray')*NF;
        binary_data=dec2bin(decimal_data,q_bit);
        data_mod=zeros(1644,560,2);
        index1=1;
        for symbol=1:560
            if mod(symbol-3,14)== 0
                data_range= reshape(data_mod_L(index1:index1+ 822*Tx_num-1),1, 822,Tx_num );
                data_mod(:,symbol,:) = reshape([ data_range  ; DMRS_value * ones(1,822,Tx_num)],1644,2);
                index1= index1 + 822*Tx_num;
            else
                data_mod(:,symbol,:) = reshape(data_mod_L(index1:index1+1644*Tx_num-1),1644,Tx_num );
                index1= index1 + 1644*Tx_num;
            end
        end
        %CDM
        data_mod(2:4:1644,3:14:560,2:2:Tx_num) = -data_mod(2:4:1644,3:14:560,2:2:Tx_num);
        %GB
        GBhead=zeros(202,560,Tx_num);
        GBtail=zeros(201,560,Tx_num);
        DC =zeros(1,560,Tx_num);
        X =[GBhead;data_mod(1:822,:,:) ;DC ;data_mod(823:end,:,:) ;GBtail];
        %IFFT
        x = ifft(ifftshift(X,1))*sqrt(2048);	%	2048 x 14*4*10
        %CP
        xhasCP = zeros(1,1228800,Tx_num);
        index1=1;
        for symbol=1:14*4*10
            if	mod(symbol,28)-1
                xhasCP(1, index1:index1+2048+144-1,:)=[ x(2048-144+1:2048,symbol,:) ; x(:,symbol,:)];
                index1 = index1+2048+144;
            else
                xhasCP(1, index1:index1+2048+208-1,:)=[ x(2048-208+1:2048,symbol,:) ; x(:,symbol,:)];
                index1 = index1+2048+208;
            end
        end
        %通道
        PowerdB = [-2 -8 -10 -12 -15 -18].';
        PowerdB_MIMO= repmat(PowerdB,1,Rx_num,Tx_num);
        H_Channel = sqrt(10.^(PowerdB_MIMO./10)).* sqrt(1/Tx_num );
        Ntap=6;
        H_Channel = H_Channel .* ( sqrt( 1/2 ).* ( randn(Ntap,Rx_num,Tx_num) + 1i*randn(Ntap,Rx_num,Tx_num) ) );
        %convolution
        y= zeros(Rx_num,1228800+5);
        for i= 1:Tx_num
            for j= 1:Rx_num
                y(j,:)=y(j,:)+conv(xhasCP(1,:,i),H_Channel(:,j,i) );
            end
        end
        newy=y(:,1:end-5); 
        %產生訊號
        n = sqrt(No/2) *( randn(Tx_num,1228800) + randn(Tx_num,1228800)*1i );% randn產生noise variance=No
        y_noise= newy+ n;
        %% without CFO
        %移除CP
        yhasnoCP= zeros(2048,560,Tx_num);
        index2 = 1;
        for symbol = 1:560
            if (mod(symbol,28)-1)
                for k = 1:Rx_num
                    yhasnoCP(:,symbol,k) = y_noise(k,index2+144:index2+144+2048-1);
                end
                index2= index2+144+2048;
            else
                for k = 1:Rx_num
                    yhasnoCP(:,symbol,k) = y_noise(k,index2+208:index2+208+2048-1);
                end
                index2= index2+208+2048;
            end
        end
        %FFT----以下為頻域
        Y_fft= fftshift( fft( yhasnoCP/sqrt(2048) ) ,1);
        %rm Guard Band
        Y= [ Y_fft( 203:1024,:,:) ; Y_fft(1026:1847,:,:) ];
        %取得LMMSE估測結果
        Y_DMRS= Y(2:2:1644 ,3:14:560,:);
        P0 = eye(822);
        P1 = ones(822,1);
        P1(1:2:822) = -1;
        P1 = diag(P1);
%         for i = 1:2:822
%             P1(i,i) = -1;
%         end
        H_LMMSE= zeros(1644,40,2,2);
        H_LMMSE(:,:,1,1)=R_H_HD* P0'*inv( P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS(:,:,1));
        H_LMMSE(:,:,2,1)=R_H_HD* P0'*inv( P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS(:,:,2));
        H_LMMSE(:,:,1,2)=R_H_HD* P1'*inv( P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS(:,:,1));
        H_LMMSE(:,:,2,2)=R_H_HD* P1'*inv( P0*R_HD_HD*P0'+ P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS(:,:,2));
        %線性內差
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
        %雜訊估測
        DMRS= data_mod(2:2:1644 , 3:14:560,:);
        DMRS_H= H_LMMSE (2:2:1644,:,:,:);
        DMRS_hat_Y= zeros(822,40,2);
        for carrier = 1:822
            for symbol = 1:40
                DMRS_hat_Y(carrier,symbol,:) = reshape(DMRS_H(carrier,symbol,:,:),2,2) * reshape(DMRS(carrier,symbol,:),2,1);
            end
        end
        Rx1_No=sum(abs( Y_DMRS(:,:,1)-DMRS_hat_Y(:,:,1) ).^2 ,'all' )/(40*822);
        Rx2_No=sum(abs( Y_DMRS(:,:,2)-DMRS_hat_Y(:,:,2) ).^2 ,'all' )/(40*822);
        NM_Y = zeros(1644,560,2);
        NM_H = zeros(1644,560,2,2);
        NM_Y(:,:,1)= Y(:,:,1)./ Rx1_No;
        NM_Y(:,:,2)= Y(:,:,2) ./ Rx2_No;
        NM_H(:,:,1,:)= H_LI(:,:,1,:)./ Rx1_No;
        NM_H(:,:,2,:)= H_LI(:,:,2,:)./ Rx2_No;
        X_hat_L=zeros(1644,560,Tx_num);
        % ZF
        for carrier=1:1644
            for slot =1:560
                zf=zeros(Tx_num,1);
                h_temp2=squeeze(NM_H(carrier,slot,:,:));  %變成二維
                y_temp2=squeeze(NM_Y(carrier,slot,:));      %變成一維
                zf=zf+(inv(h_temp2'*h_temp2)*h_temp2')*y_temp2;      %inv(H'H)H'y
                X_hat_L(carrier,slot,:)=zf;

            end

        end
        X_hat_L = X_hat_L/NF;
        %反解
        data_mod_L_hat = zeros(Tx_num*(1644*560-822*40),1);
        index1= 1;
        for symbol=1:560
            if mod( (symbol - 3) , 14) == 0
                data_mod_L_hat(index1:index1+ 822*Tx_num-1)	=reshape(X_hat_L(1:2:1644,symbol,:),822*Tx_num,1 );
                index1= index1 + 822*Tx_num;
            else
                data_mod_L_hat(index1:index1+1644*Tx_num-1)	= reshape(X_hat_L(1:1:1644,symbol,:),1644*Tx_num,1 );
                index1= index1 + 1644*Tx_num;
            end
        end

        %demod
        data_dec_L_hat=qamdemod(data_mod_L_hat,QAM,'gray');
        data_bin_L_hat=dec2bin(data_dec_L_hat,q_bit);
        %BER計算
        BER= BER+ sum(sum(data_bin_L_hat ~= binary_data ));
        %% 頻偏
        t = (10e-3/1228800):(10e-3/1228800):10e-3;
        y_noise= exp( 1i * 2 * pi * delta_f * t) .* y_noise;

        index1 = 1;
        CFO_sum= 0;
        CP_num= 0;
        for symbol = 1:560
            if (mod(symbol,28)-1)
                CFO_head=y_noise(:,index1: index1+143);
                CFO_tail=y_noise(:, index1+2048: index1+143+2048 );
                CFO_sum=CFO_sum+ sum(CFO_head .* (CFO_tail').','all');
                CP_num=CP_num + 144;
                index1=index1 + 144 + 2048;
            else		
                CFO_head=y_noise( :,index1: index1+207 );
                CFO_tail=y_noise( :,index1+2048: index1+207+2048 );
                CFO_sum=CFO_sum+ sum(CFO_head .* (CFO_tail').','all');
                CP_num=CP_num + 208;
                index1=index1+208+2048;
            end

        end
        %偏移乘回來
        hat_delta_f =-(Fs/2048)*angle(CFO_sum/CP_num)/(2*pi) ;
        y_noise=y_noise./exp(1i*2*pi*hat_delta_f*t);
        %移除CP
        yhasnoCP= zeros(2048 , 560, Tx_num);
        index1  = 1;
        for symbol = 1:560
            if (mod(symbol,28)-1)
                for k = 1:Rx_num
                    yhasnoCP(:,symbol,k) = y_noise(k,index1+144:index1+144+2048-1);
                end
                index1  = index1+144+2048;
            else
                for k = 1:Rx_num
                    yhasnoCP(:,symbol,k) = y_noise(k,index1+208:index1+208+2048-1);
                end
                index1  = index1 +208+2048;
            end
        end
        %FFT----以下為頻域
        Y_fft2 = fftshift( fft( yhasnoCP/sqrt(2048) ) ,1);
        %rm GB
        Y2= [ Y_fft2( 203:1024,:,:) ; Y_fft2( 1026:1847,:,:) ];
        Y_DMRS2=Y2(2:2:1644,3:14:560,:);
        H_LMMSE2=zeros(1644,40,2,2);
        H_LMMSE2(:,:,1,1)=R_H_HD*P0'*inv(P0*R_HD_HD*P0'+P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS2(:,:,1));
        H_LMMSE2(:,:,2,1)=R_H_HD*P0'*inv(P0*R_HD_HD*P0'+P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS2(:,:,2));
        H_LMMSE2(:,:,1,2)=R_H_HD*P1'*inv(P0*R_HD_HD*P0'+P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS2(:,:,1));
        H_LMMSE2(:,:,2,2)=R_H_HD*P1'*inv(P0*R_HD_HD*P0'+P1*R_HD_HD*P1'+(1/SNR_Power)*P0)*(DMRS_value'.*Y_DMRS2(:,:,2));
        %線性內差
        H_LI2=zeros(1644,560,2,2);
        DMRS_ps2=3:14:560;
        for symbol=1:560
            if symbol>=1 && symbol<=2  %前邊界
                foreward=11+symbol;
                backward=3-symbol;
                H_LI2(:,symbol,:,:)=(foreward*H_LMMSE2(:,1,:,:))/14+(backward*H_LMMSE2(:,40,:,:))/14;
            end

            if symbol>=3 && symbol<=548 % 中間
                location=floor((symbol-3)/14)+1;
                foreward=symbol-DMRS_ps2(location);
                backward=DMRS_ps2(location+1)-symbol;
                H_LI2(:,symbol,:,:)=(foreward*H_LMMSE2(:,location+1,:,:))/14+backward*H_LMMSE2(:,location,:,:)/14;
            end
            if symbol>=549 && symbol<=560 %後邊界
                foreward=symbol-549;
                backward=563-symbol;
                H_LI2(:,symbol,:,:)=(foreward*H_LMMSE2(:,1,:,:))/14+(backward*H_LMMSE2(:,40,:,:))/14;
            end
        end

        %雜訊估測
        DMRS2= data_mod(2:2:1644 , 3:14:560,:);
        DMRS_H2=H_LMMSE2(2:2:1644,:,:,:);
        DMRS_hat_Y2=zeros(822,40,2);
        for carrier=1:822
            for symbol=1:40
                DMRS_hat_Y2(carrier,symbol,:) = reshape(DMRS_H2(carrier,symbol,:,:),2,2) * reshape(DMRS2(carrier,symbol,:),2,1);
            end
        end
        Rx1_No2 = sum( abs( Y_DMRS2(:,:,1)  - DMRS_hat_Y2(:,:,1) ).^2 ,'all' )/(40*822);
        Rx2_No2 = sum( abs( Y_DMRS2(:,:,2)  - DMRS_hat_Y2(:,:,2) ).^2 ,'all' )/(40*822);

        NM_Y2 = zeros(1644,560,2);
        NM_H2 = zeros(1644,560,2,2);
        NM_Y2(:,:,1) = Y2(:,:,1)./ Rx1_No2;
        NM_Y2(:,:,2) = Y2(:,:,2)./ Rx2_No2;
        NM_H2(:,:,1,:) = H_LI2(:,:,1,:) ./ Rx1_No2;
        NM_H2(:,:,2,:) = H_LI2(:,:,2,:) ./ Rx2_No2;
        X_hat_L2=zeros(1644,560,Tx_num);
        % ZF
        for carrier=1:1644
            for slot =1:560
                zf=zeros(Tx_num,1);
                h_temp=squeeze(NM_H2(carrier,slot,:,:));  %變成二維
                y_temp=squeeze(NM_Y2(carrier,slot,:));      %變成一維
                zf=zf+(inv(h_temp'*h_temp)*h_temp')*y_temp;      %inv(H'H)H'y
                X_hat_L2(carrier,slot,:)=zf;
            end
        end
        X_hat_L2=X_hat_L2/NF;
        %反解資料
        data_mod_L_hat2 =zeros(Tx_num*(1644*560-822*40),1);
        index1= 1;
        for symbol=1:560
            if mod(symbol-3, 14) == 0
                data_mod_L_hat2(index1:index1+ 822*Tx_num-1)= reshape(X_hat_L2(1:2:1644,symbol,:), 822*Tx_num,1 );
                index1= index1 + 822*Tx_num;
            else
                data_mod_L_hat2(index1:index1+1644*Tx_num-1)= reshape(X_hat_L2(1:1:1644,symbol,:),1644*Tx_num,1 );
                index1	= index1 + 1644*Tx_num;
            end
        end
        %demod
        data_dec_L_hat2=qamdemod(data_mod_L_hat2,QAM,'gray');
        data_bin_L_hat2=dec2bin (data_dec_L_hat2,q_bit);     	
        %BER計算
        BER2=BER2+sum(sum(data_bin_L_hat2~= binary_data));
    end
    BER_SNR(1,count)=BER/(1644*560*frame_num*q_bit*Tx_num);
    BER_SNR2(1,count)=BER2/(1644*560*frame_num*q_bit*Tx_num);
end

semilogy(SNR_in_dB,BER_SNR(1,:),'b-','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR2(1,:),'r--','LineWidth',2)
hold on
grid on
title('BER of 5G NR MIMO OFDM(CFO)')
xlabel('SNR (dB)')
ylabel('BER')
legend('CFO','without CFO')
