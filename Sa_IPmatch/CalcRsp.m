function [rsp,disp,tpeak,tindex]=CalcRsp(acc,dt,ksi,T)

% This function computes the dynamic response of a single-degree-of-freedom system 
% using Newmark's Method.

% Input Variables
% acc   : acceleration time history (in g)
% dt    : Time step of the record
% ksi   : damping ratio
% T     : period 

% Output Variables
% rsp   : true value of peak pseudo-acceleration
% disp  : displacement time history
% tpeak : time of peak response
% tindex: index of the tpeak

acc=reshape(acc,1,[]);
acc=acc(all(~isnan(acc),2),:);
m=1; % mass
wn=2*pi./T; 
a=acc;
t=0:dt:length(a)*dt-dt;

% Newmark coefficients
gamma=0.5;
beta=1/4; % average acceleration

u=zeros(length(T),length(t));
tpeak=zeros(1,length(T));
disp=cell(1,length(T));
rsp=zeros(1,length(T));
tindex=zeros(1,length(T));

% acceleration response
for i=1:length(T)
udot=0;
ks=m.*wn(i)^2;
c=2*m*wn(i)*ksi;
uddot=(a(1)-c.*udot-ks.*u(i,1))/m;

a1=1/(beta*dt.^2).*m+gamma/(beta*dt).*c;
a2=1/(beta*dt).*m+(gamma/beta-1).*c;
a3=(1/2/beta-1).*m+dt*(gamma/2/beta-1).*c;
kh=ks+a1;

for j=1:length(t)-1
    ph=a(j+1)+a1*u(i,j)+a2*udot+a3*uddot;
    u(i,j+1)=ph/kh;
    udoti=udot;
    udot=gamma/(beta*dt).*(u(i,j+1)-u(i,j))+(1-gamma/beta).*udot+dt*(1-gamma/2/beta)*uddot;
    uddot=1/(beta*dt^2).*(u(i,j+1)-u(i,j))-1/(beta*dt).*udoti-(1/2/beta-1)*uddot;
end

PA=abs(u(i,:))*wn(i)^2; % absolute pseudo-acceleration
R=max(PA); % peak of absolute pseudo-acceleration
tindex(i)=find(PA==R);
tpeak(i)=(tindex(i)-1).*dt; 
disp{i}=u(i,:);
rsp(i)=u(i,tindex(i)).*wn(i)^2;
end
