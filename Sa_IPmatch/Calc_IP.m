function [IP_rec,tIP,tIP_start,tIP_end]=Calc_IP(Vel,dt,T)
% This function calculates the IP of the record.

% Input Variables
% Vel       : velocity time history (in cm/s)
% dt        : time step of the record (s) 
% T         : Period 

% Output Variables
% IP_rec    : Instantaneous Power of a record (in cm^2/s^2)
% tIP       : IP center time 
% tIP_start : IP start time 
% tIP_end   : IP end time

t=0:dt:length(Vel)*dt-dt;
fcutlow=1./(3*T);
fcuthigh=1./(0.2*T);
nOrder=4;

Vel_filt=bandpass_filter(Vel,dt,fcutlow,fcuthigh,nOrder);

for i=1:length(T)
    
twin=round(0.5*T./dt);

for k=1:length(t)-twin(i)+1   
    energy_Vel(i,k)=(sum(Vel_filt(i,k:k+twin(i)-1).^2)).*dt/(0.5*T(i));
end

[IP_rec(i),idx(i)]=max(energy_Vel(i,:));

tIP(i)=(idx(i)-1).*dt+0.25*T(i); 
tIP_start(i)=tIP(i)-0.25*T(i); 
tIP_end(i)=tIP(i)+0.25*T(i); 
end
