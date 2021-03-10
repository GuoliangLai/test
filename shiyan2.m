clc
clear 
close all
%ASK(PAM)
M=[2,4,8];
r1=.2*log2(M);
x1=qfuncinv(0.5*0.0001)%算出Q函数对应X的值
SNR1=10*log10(x1*x1./(6*log2(M)./(M.*M-1)));
subplot(121);
plot(SNR1, r1,'-*');
ylabel('(R/W)/(bit/s)/Hz');xlabel('比特SNR(dB)')
hold on;
%PSK
M=[2,4,8,16];
SNR2=zeros(1,4);
r2=zeros(1,4);
SNR2(1:2)=SNR1(1);
x2=qfuncinv(0.5*0.000001)
 for n=1:4
     r2(n)=log2(M(n));
    if (n==1 || n==2)
      SNR2(n)=SNR1(1);
    else
     SNR2(n)=10*log10(x2*x2/(2*pi*pi*log2(M(n))/(M(n)*M(n))));
    end
 end
plot(SNR2, r2,'-+');
%DPSK
hold on;
SNR3=zeros(1,3);
r3=r2(1:3);
for n=1:3
    if (n==1)
      SNR3(n)=SNR2(1)+0.8;
    else
     SNR3(n)=SNR2(n)+2.3;
    end
 end
plot(SNR3, r3,'-o');
 
%QAM
hold on;
M1=[4,8,16,64];
x4=qfuncinv(0.0000001)
SNR4=10.*log10(x4*x4./(3*log2(M1)./(M1-1)));
r4=log2(M1);
plot(SNR4, r4,'-x');
title('无记忆调制-信研2008赖国良')
legend('(ASK)PAM','PSK','DPSK','QAM');

%正交信号相干检测
M2=[8,16,32,64];
r5=2.*log2(M2)./M2(n);
SNR5=[6 6.5 7 8.2];
subplot(122);
plot(SNR5,r5,'-s');
title('正交信号相干检测-信研2008赖国良')
ylabel('(R/W)/(bit/s)/Hz');xlabel('比特SNR(dB)')