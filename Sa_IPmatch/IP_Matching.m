function [acc]=IP_Matching(acc,dt,T,targetIP)
% This function performs an IP adjustment.

gamma_vel=0.10; %relaxation factor
w=2*pi./T;
b2=ones(length(T),1);

t=0:dt:length(acc)*dt-dt;
Vel_matched=cumsum(acc*981)*dt;

[~,tIPin,~,~]=Calc_IP(Vel_matched,dt,T);
wav_IP=zeros(length(T),length(Vel_matched));

% Ricker wavelet
for i=1:length(T)
    tshift=tIPin(i);
    tj=t-tshift; % centered on IP time
    wav_IP(i,:)=(1-0.5*(w(i).^2.*tj.^2)).*exp(-0.25*(w(i).^2*(tj).^2));
end

% Filtered wavelet
for i=1:length(T)
    fcutlow=1./(3*T);
    fcuthigh=1./(0.2*T);
    nOrder=4;
    wav_IP_filt(i,:)=bandpass_filter(wav_IP(i,:),dt,fcutlow(i),fcuthigh(i),nOrder);

end

% Adjusted Velocity

for k=1:length(T)
    dVel(k,:)=b2(k).*wav_IP(k,:);
    dVel_tot=sum(dVel,1);
end

Vel_matched=Vel_matched+gamma_vel*dVel_tot; 

% IP of the adjusted record
[IP,tIP,tIP_start,tIP_end]=Calc_IP(Vel_matched,dt,T);

% Filter the adjusted velocity
fcutlow=1./(3*T);
fcuthigh=1./(0.2*T);
nOrder=4;
Vel_filt=bandpass_filter(Vel_matched,dt,fcutlow,fcuthigh,nOrder);

% IP misfit after the adjustment
for m=1:length(T)
    misfit_IP(m)=targetIP(m)-IP(m);
    logmisfit_IP(m)=log(targetIP(m))-log(IP(m));
    
    [~,idx(m)]=min(abs(t-tIP(m)));
    [~,idx1(m)]=min(abs(t-tIP_start(m)));
    [~,idx2(m)]=min(abs(t-tIP_end(m)));
    
end  

%%%%C-Matrix %%%%
for i=1:length(T) % each period
    for k=1:length(T) % each time of the wavelet   
        dIP{i,k}=(2*Vel_filt(i,:).*wav_IP_filt(k,:));
        C2_mat{i,k}=((sum(dIP{i,k}(idx1(i):idx2(i)))).*dt)/(0.5*T(i));

    end
end

C2_mat=cell2mat((C2_mat));
diagonal_coef=0.30; 
C2_mat(setdiff(1:numel(C2_mat),1:length(C2_mat)+1:numel(C2_mat)))=diagonal_coef.*C2_mat(setdiff(1:numel(C2_mat),1:length(C2_mat)+1:numel(C2_mat)));

%%% Tikhonov regularization
lambda=0.10*max(abs((diag(C2_mat))));
tikh_m = lambda*eye(size(C2_mat,1),size(C2_mat,2));
tikh_mT= tikh_m';
tikh_mProd = tikh_mT*tikh_m;
        
C2_mat_T = C2_mat';
C2_mat_Prod = C2_mat_T*C2_mat;
A=cat(3,C2_mat_Prod,tikh_mProd);
D=sum(A,3);
b2=(inv(D))*C2_mat'*misfit_IP'; 

for i=1:length(T)
tshift=tIP(i);
tj=t-tshift;
wav_IP(i,:)=(1-0.5*(w(i).^2.*tj.^2)).*exp(-0.25*(w(i).^2*(tj).^2)); 
end

for k=1:length(T)
    dVel(k,:)=b2(k).*wav_IP(k,:);
    dVel_tot=sum(dVel,1);
end

Vel_matched=Vel_matched+gamma_vel.*dVel_tot;

% Adjusted acceleration
acc=diff(Vel_matched)./dt;
acc=[0,acc];
acc=acc*0.01/9.81; % convert cm/s^2 to g-unit.

