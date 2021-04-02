function [acc]=Sa_Matching(acc,dt,ksi,T,target)
% This function performs a spectral adjustment.

acc=reshape(acc,1,[]);
gamma_acc=0.50; % relaxation factor
w=2*pi./T;
a=acc;
t=0:dt:length(a)*dt-dt;

%% Sa misfit
[rsp,~,tpeak,tindex]=CalcRsp(a,dt,ksi,T);

for i=1:length(T)
    if rsp(i)>=0
        target(i)=abs(target(i));
    else
        target(i)=-abs(target(i));
    end
    misfit(i)=target(i)-rsp(i);
end
%% b computation
% Cosine-tapered wavelet function
dTj=zeros(1,length(T));
wavt=zeros(length(T),length(a));
for i=1:length(T)
    
    wj=w(i)*sqrt(1-ksi^2);
    tshift=tpeak(i);
    dTj(i)=atan((sqrt(1-ksi^2))/(ksi))/wj;
    tj=t-tshift+dTj(i);
    
    a1w=1.25;
    a2w=0.25;
    f1=1.0;
    f2=4.0;
    if 1/T(i)<1
        alpha=a1w*(1/T(i));
    elseif (f1<1/T(i)) && (1/T(i)<f2)
        alpha=(a1w+(a2w-a1w)*((1/T(i))-f1)/(f2-f1))*(1/T(i));
    else
        alpha=a2w*(1/T(i));
    end
    
    wavt(i,:)=(cos(wj.*tj).*exp(-abs(tj).*alpha)); %adjustment function
    
end

%apply highpass filter to adjustment function
for i=1:length(T)
    fc=1/6;
    nOrder=4;
    N = 2^nextpow2(length(wavt(i,:)));
    Y = fft(wavt(i,:),N);
    df=1/(N*dt);
    
    freq=0:df:N/2*df-df;
    ff=freq./fc;
    H=zeros(1,length(N));
    for j=2:length(ff)
        
        H(j)= sqrt((ff(j).^(2*nOrder)))./sqrt(((1+ff(j).^(2*nOrder))));
        Y(j)= Y(j).*H(j);
        Y(length(Y)-j+2)=Y(length(Y)-j+2).*H(j);
    end
    
    Y1= (real(ifft(Y,N)))';
    Y1=Y1(1:length(wavt));
    wavt(i,:)=Y1;
end

% C-matrix
for i=1:length(T) %each period
    for k=1:length(T)
        [~,disp(i,k),~]=CalcRsp(wavt(i,:),dt,ksi,T(k));
        C_mat(i,k)=(disp{i,k}(1,tindex(i))).*w(i)^2;
    end
end

%%%off-diagonal reduction of matrix
C1_mat=C_mat;
diagonal_coef=0.30;
C1_mat(setdiff(1:numel(C1_mat),1:length(C1_mat)+1:numel(C1_mat)))=diagonal_coef.*C1_mat(setdiff(1:numel(C1_mat),1:length(C1_mat)+1:numel(C1_mat)));

%%%Tikhonov regularization
lambda = 0.10*max(abs((diag(C1_mat))));
tikh_m = lambda*eye(size(C1_mat,1),size(C1_mat,2));
tikh_mT= tikh_m;
tikh_mProd = tikh_mT*tikh_m;

C1_mat_T = C1_mat';
C1_mat_Prod = C1_mat_T*C1_mat;
A1=cat(3,C1_mat_Prod,tikh_mProd);
D1=sum(A1,3);

b=(inv(D1))*C1_mat'*misfit';

for i=1:length(T)
    adj(i,:)=b(i).*wavt(i,:); % adjustment
end
adj_final=sum(adj,1);

acc=acc+gamma_acc.*adj_final; % adjusted acceleration

% % Apply lowpass filter to adjusted acceleration
sampleRate= 1/dt;
fc= 1/0.05; % corner frequency
filterOrder =4;
[b1,a1] = butter(filterOrder,fc/(sampleRate/2),'low');
acc=filtfilt(b1,a1,acc);
