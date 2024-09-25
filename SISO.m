% SISO system
clear
clc
close all

Tx =1;         %�ǰe�ݭӼ�
Rx =1;         %�����ݭӼ�

%�i�H���ۤv�]
data_num = 150000;         % data�q(�`�NBER�n�]��10^(-3)!!)
SNR_in_dB = 0:5:40;        % �ۤv�]�q���T�j�p

SER_SNR_ZF=zeros(3,length(SNR_in_dB));
BER_SNR_ZF=zeros(3,length(SNR_in_dB));
SER_SNR_LMMSE=zeros(3,length(SNR_in_dB));
BER_SNR_LMMSE=zeros(3,length(SNR_in_dB));

%% choose QAM= 4/16/64;
for v=1:3
    QAM = 4^v;
    x=(0:QAM-1);              % 0��M���ƦC
    y=qammod(x,QAM);
    Eavg = mean(abs(y).^2);     %average power
    NF = 1/sqrt(Eavg);           % normal factor
    q_bit = log2(QAM);        % �@��symbol�i�H�ǴX��bit
    N = Tx;

    for  a=1:length(SNR_in_dB)
        SNR = 10^(SNR_in_dB(a)/10);
        No = 10^(-SNR_in_dB(a)/10);
        Es = 1;
        
        SER_ZF = 0;     
        BER_ZF = 0;                         % ��error rate�n�@����
        SER_LMMSE = 0;  
        BER_LMMSE = 0;
        
        parfor seperate = 1:data_num           %
            data = randi([0 QAM-1]);          % �H������0~3 for 4QAM
            bin_data = dec2bin(data,q_bit);      % �N 0~3 �ର '00'~'11'
            X = qammod(data, QAM,'gray') * NF;      % 0~3 to complex (Modulation); remember to normalize
            H = randn(Rx, Tx) / sqrt(2) + 1i * randn(Rx, Tx) / sqrt(2);     % randn����channel(�`�N���W�ƪ����D)
            n = sqrt(No/2) * (randn(Rx, 1) + 1i * randn(Rx, 1));   % randn����noise variance=No
            
            Y = H * X + n;
            
            %% type 1 = ZFs
            X_hat_ZF  =inv(H)*Y;
            data_hat_ZF = qamdemod( X_hat_ZF/NF , QAM, 'gray');     %complex to 0~3  ; remember to inverse-normalize
            bin_data_hat = dec2bin(data_hat_ZF, q_bit);          %0~3 to '00'~'11'
            % ��SER/BER
            BER_ZF =BER_ZF + sum(bin_data_hat ~= bin_data) ;
            if data_hat_ZF~=data 
                SER_ZF = SER_ZF + 1;
            end
            
            %% type 2 = LMMSE
            X_hat_LMMSE =  inv(H' * H + (No * eye(Tx))) * (H' * Y);
            data_hat_LMMSE = qamdemod(X_hat_LMMSE/NF , QAM, 'gray');        % complex to 0~3 ; remember to inverse-normalize
            bin_data_hat_LMMSE =  dec2bin(data_hat_LMMSE, q_bit);    %0~3 to '00'~'11'
            % SER/BER
            BER_LMMSE = BER_LMMSE + sum(bin_data_hat_LMMSE ~= bin_data);
            if data_hat_LMMSE~=data 
                SER_LMMSE = SER_LMMSE + 1;
            end
        end
    
%         ����SNR���n��SER/BER�s�b�x�}��
        SER_SNR_ZF(v,a) = SER_ZF/data_num;
        BER_SNR_ZF(v,a) = BER_ZF/(data_num*q_bit);
        
        SER_SNR_LMMSE(v,a) = SER_LMMSE/data_num;
        BER_SNR_LMMSE(v,a) = BER_LMMSE/(data_num*q_bit);
        
    end
end
% 
% ��JSNR_in_dB�MSER
figure(1)
semilogy(SNR_in_dB,SER_SNR_ZF(1,:),'r-X','LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_ZF(2,:),'r-diamond','LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_ZF(3,:),'r-O','LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(1,:),'b-X','LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(2,:),'b-diamond','LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(3,:),'b-O','LineWidth',2)
hold on
grid on
title('SER of SISO')
xlabel('SNR (dB)')
ylabel('SER')
legend('4QAM ZF','16QAM ZF','64QAM ZF','4QAM LMMSE','16QAM LMMSE','64QAM LMMSE')

%��JSNR_in_dB�MBER
figure(2)
semilogy(SNR_in_dB,BER_SNR_ZF(1,:),'r-X','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_ZF(2,:),'r-diamond','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_ZF(3,:),'r-O','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(1,:),'b-X','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(2,:),'b-diamond','LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(3,:),'b-O','LineWidth',2)
hold on
grid on
title('BER of SISO')
xlabel('SNR (dB)')
ylabel('BER')
legend('4QAM ZF','16QAM ZF','64QAM ZF','4QAM LMMSE','16QAM LMMSE','64QAM LMMSE')

