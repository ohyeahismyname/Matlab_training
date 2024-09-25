
clear
close all
warning('off')

QAM=4;
SNR_in_dB = 0:1:16;
BER_SNR=zeros(1,length(SNR_in_dB));
BER_SNR2=zeros(1,length(SNR_in_dB));


xxx=(0:QAM-1);              % 0到M的數列
yyy=qammod(xxx,QAM);
Eavg = mean(abs(yyy).^2);     %average power
NF= 1/sqrt(Eavg);
q_bit = log2(QAM);        % 一個symbol可以傳幾個bit

for a=1:length(SNR_in_dB)
    SNR = 10^( SNR_in_dB(a)/10);
    No= 10^(-SNR_in_dB(a)/10);
    n = sqrt(No/2)*(randn(1,1228800)+1i*randn(1,1228800));
    BER=0;
    BER2=0;


    decimal_data = randi([0, 3], 12, 14);
    ex_decimal_data=repmat (decimal_data,137,40);
    binary_data=dec2bin(bin2gray(ex_decimal_data,'qam',QAM),q_bit);
    GBhead=zeros(202,560);
    GBtail=zeros(201,560);
    DC=zeros(1,560);
    gray_moded=qammod(ex_decimal_data,QAM,'gray')*NF;
    TX=[GBhead;gray_moded(1:822,:);DC;gray_moded(823:1644,:);GBtail];
    x = ifft(ifftshift(TX,1))*sqrt(2048) ;

    xhasCP=zeros(1,1228800);
    start=1;
    for symbol=1:560
        if mod(symbol,28) ==0
            xhasCP(1,start:start+208)=[x(2048-208+1:2048),symbol];
            xhasCP(1,start+208:start+2048+208-1)=[x(1:2048,symbol)];
            start=start+2048+208;
        else
            xhasCP(1,start:start+144)=[x(2048-144+1:2048),symbol];
            xhasCP(1,start+144:start+2048+144-1)=[x(1:2048,symbol)];
            start=start+2048+144;
        end

    end

    PowerdB=[-2 -8 -10 -12 -15 -18];
    Total_H_Power=sum(10.^(PowerdB/10));%總能量為一
    H_Channel=sqrt(10.^(PowerdB/10));%db轉增益
    %%
    y=conv(xhasCP,H_Channel);
    newy=y(1:end-5);
    noisey=newy+n;
    yhasnoCP=zeros(2048,560);

    start2=1;
    for symbol=1:560
        if mod(symbol,28) ==0
            yhasnoCP(:,symbol)=noisey(1,start2+208:start2+2048+208-1);
            start2=start2+2048+208;
        else
            yhasnoCP(:,symbol)=noisey(1,start2+144:start2+2048+144-1);
            start2=start2+2048+144;
        end
    end



    Y = fftshift( fft( yhasnoCP/sqrt(2048) ) ,1);
    YhasnoGB=[ Y( 203:1024,:) ; Y( 1026:1847,:) ];
    h=[H_Channel,zeros(1,2042)];
    H=fftshift(fft(h));%時域轉頻域
    H_Data=[H(1,202:1023),H(1,1025:1846)].';
    H_frame=repmat(H_Data,1,560); %frame的理想通道
    X_hat=(YhasnoGB./H_frame)/NF; %ZF detector

    decimal_data_hat = qamdemod(X_hat,QAM,'gray');
    binary_data_hat=dec2bin(bin2gray(decimal_data_hat,'qam',QAM),q_bit);

    BER=BER+sum(sum(binary_data~= binary_data_hat));
    BER_SNR(1,a) =BER/(1841280);
    %%

    yy=CONV(xhasCP,H_Channel);
    newyy=yy(1:end-5);
    noiseyy=newyy+n;
    yyhasnoCP=zeros(2048,560);
    start3=1;
    for symbol=1:560
        if mod(symbol,28) ==0
            yyhasnoCP(:,symbol)=noiseyy(1,start3+208:start3+2048+208-1);
            start3=start3+2048+208;
        else
            yyhasnoCP(:,symbol)=noiseyy(1,start3+144:start3+2048+144-1);
            start3=start3+2048+144;
        end
    end
    YY = fftshift( fft( yyhasnoCP/sqrt(2048) ) ,1);
    YYhasnoGB=[ YY( 203:1024,:) ; YY( 1026:1847,:) ];
    h=[H_Channel,zeros(1,2042)];
    H=fftshift(fft(h));%時域轉頻域
    H_Data=[H(1,202:1023),H(1,1025:1846)].';
    H_frame=repmat(H_Data,1,560); %frame的理想通道
    XX_hat=(YYhasnoGB./H_frame)/NF; %ZF detector

    decimal_data_hat2 = qamdemod(XX_hat,QAM,'gray');
    binary_data_hat2=dec2bin(bin2gray(decimal_data_hat2,'qam',QAM),q_bit);


    BER2=BER2+sum(sum(binary_data~= binary_data_hat2));
    BER_SNR2(1,a) =BER2/(1841280);

end
semilogy(SNR_in_dB,BER_SNR(1,:),'r-X','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR2(1,:),'b--o','LineWidth',2)
hold on
grid on
title('BER of SISO')
xlabel('SNR (dB)')
ylabel('BER')
legend('4QAM ZF matlab','4QAM ZF C')


