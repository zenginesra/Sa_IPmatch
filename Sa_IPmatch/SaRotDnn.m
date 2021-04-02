function [SaH1,SaH2,RotD50,RotD100]=SaRotDnn(accH1,accH2,dt,ksi,T)

% Output Variables
% RotD50: median single-component horizontal ground motion across all
% non-redundant azimuths
% RotD100: maximum rotated component
% SaH1: pseudo acceleration spectrum (H1 component)
% SaH2: pseudo acceleration spectrum (H2 component)

[rsp1,disp1,~,~]=CalcRsp(accH1,dt,ksi,T);
[rsp2,disp2,~,~]=CalcRsp(accH2,dt,ksi,T);

SaH1=abs(rsp1); 
SaH2=abs(rsp2); 

Sa=zeros(180,length(T));

for j=1:length(T)
    for i = 1:length(Sa) 
        
        x = i * pi/180;
        u_dir{j}= disp1{1,j}.*cos(x) + disp2{1,j}.*sin(x);
        Sa(i,j) = max(abs(u_dir{1,j})).*(2*pi/T(j))^2;
                
    end
end

RotD50 = median(Sa);
RotD100 = max(Sa);
