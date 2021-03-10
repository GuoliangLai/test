function[]=RateDF(Pa,d,S)
format long 
d=input('ʧ�����d=');
Pa=input('������ʷֲ� Pa=');
r=input('������Դ��r=');
s=input('�����Դ��s=');
S=input('��ʽ����S=');
times=input('��������times=');
[r,s]=size(d);
if(length(find(Pa<=0))~=0)
    error('Not a prob.vector,shoud be positive component!');
end
if(abs(sum(Pa)-1)>10e-10)
    error('Not a prob.vector,component do not add up to 1!')
end
if(r~=length(Pa))
    error('The parameters do not match!');
end
pba=[];
RS=[];
DS=[];
m=1;
for z= 1: times
    Pba(1:r,1:s,1)=1/s*ones(r,s);
    for j=1:s
        Pb(j,1)=0;
        for i=1:r
            Pb(j,1)=Pb(j,1)+Pa(i)*Pba(i,j,1);
        end
    end
    for i=1:r
        temp(i)=0;
        for j=1:s
            temp(i)=temp(i)+Pb(j,1)*exp(S(m)*d(i,j));
        end
    end
    for i=1:r
        for j=1:s
            Pba(i,j,2)=(Pb(j,1)*exp(S(m)*d(i,j)))/temp(i);
        end
        D(1)=0;
        for i=1:r
            for j=1:s
                D(1)=D(1)+Pa(i)*Pba(i,j,1)*d(i,j);
            end
        end
        R(1)=0;
        for i=1:r
            for j=1:s
                if(Pba(i,j,1)~=0)
                    R(1)=R(1)+Pa(i)*Pba(i,j,1)*log2(Pba(i,j,1)/Pb(j,1));
                end
            end
        end
        n=2;
        while(1)
            for j=1:s
                Pb(j,n)=0;
                for i=1:r
                    Pb(j,n)=Pb(j,n)+Pa(i)*Pba(i,j,n);
                end
            end
            for i=1:r
                temp(i)=0;
                for j=1:s
                   % disp('SM:');disp(S(m));
                    temp(i)=temp(i)+Pb(j,n)*exp(S(m)*d(i,j));
                end
            end
            for i=1:r
                for j=1:s
                    if(temp(i)~=0)
                        Pba(i,j,n+1)=(Pb(j,n)*exp(S(m)*d(i,j)))/temp(i);
                    end
                end
            end
             D(n)=0;
             for i=1:r
                 for j=1:s
                     D(n)=D(n)+Pa(i)*Pba(i,j,n)*d(i,j);
                 end
             end
             R(n)=0;
             for i=1:r
                 for j=1:s
                     if(Pba(i,j,n)~=0)
                       R(n)=R(n)+Pa(i)*Pba(i,j,n)*log2(Pba(i,j,n)/Pb(j,n));
                     end
                 end
             end
             %disp('E1:');disp(abs(R(n)-R(n-1)));
             %disp('E2:');disp(abs(D(n)-D(n-1)));
             if(abs(R(n)-R(n-1))<=10^(-7))
                 if(abs(D(n)-D(n-1))<=10^(-7))
                     break;
                 end
             end
             n=n+1;
         end
         S(m+1)=S(m)+0.5;
         if(abs(R(n)<10^(-7)))
             
         end
         pba=[Pba(:,:,:)];
         RS=[RS R(n)];
         DS=[DS D(n)];
         m=m+1;
    end
end
     [k,l,q]=size(pba);
     Pba=pba(:,:,q);
     Rmin=min(RS);
     Dmax=max(DS);
     Smax=S(m-1);
disp('������ȷ������������£�');
disp('��С��Ϣ��Rmin:');disp(Rmin);
disp('���Dmax:');disp(Dmax);
disp('���ת�Ƹ��ʷֲ�Pba:');disp(Pba);
disp('�����ʽ����Smax:');disp(Smax);
plot(DS,RS)
xlabel('�����ʧ���D')
ylabel('��Ϣ��ʧ�溯��R(D)')
title('��Ϣ��ʧ�溯��R(D)������ͼ������2008������')