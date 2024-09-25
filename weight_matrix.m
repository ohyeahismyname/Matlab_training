function W= weight_matrix(windows,SNR_weight)

N=2048;%IFFT point
%  windows=6;%window
%  SNR_weight=40;

R_H_H=zeros(822,822);
R_H_HS=zeros(1644,822);
DMRS_num=2:2:1644; %計數
real_num=1:1:1644;
%% R_H_H
for k=1:822
    for p=1:822
        if DMRS_num(k)==DMRS_num(p)
            R_H_H(p,k)=1;
        else
%          if DMRS_num(k)~=real_num(p)
            R_H_H(p,k)=(1-exp(-1i*2*pi*windows*(DMRS_num(p)-DMRS_num(k))/N))/(1i*2*pi*windows*(DMRS_num(p)-DMRS_num(k))/N);
        end
    end
end

%% R_H_HS
for k=1:1:822
    for p=1:1644
        if DMRS_num(k)==real_num(p)
            R_H_HS(p,k)=1;
        else
            R_H_HS(p,k)=(1-exp(-1*1i*2*pi*windows*(real_num(p)-DMRS_num(k))/N))/(1i*2*pi*windows*(real_num(p)-DMRS_num(k))/N);
        end
    end
end


%% dB轉回能量 SNR=10*log(Ps/Pn)=>Ps/Pn=10^(SNR/10)
SNR_Power=10^(SNR_weight/10);
W=R_H_HS*inv(R_H_H+eye(822,822)/SNR_Power);



