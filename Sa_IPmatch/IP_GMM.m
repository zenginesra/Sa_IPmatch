function [IP,sqrtIP,sigma_sqrtIP]=IP_GMM(M,Rrup,T,targetSa)
% Conditional Ground-Motion Model for Instantaneous Power
% Input Variables
% M           : moment magnitude
% Rrup        : closest distance to rupture (km)
% T           : period (sec)
% targetSa    : target spectral acceleration (g)

% Output Variables
% IP          : Median Instantaneous Power (in cm^2/s^2)
% sqrtIP      : Median square root of IP.
% sigma_sqrtIP: standard deviation of ln(sqrtIP)

% Referenced manuscript:
% Zengin E. and Abrahamson N.(2020) Conditional Ground-Motion Model for
% Damaging Characteristics of Near-Fault Ground Motion Based on
% Instantaneous Power, Bulletin of the Seismological Society of America,110(6), 2828-2842.

T_period=[0.05 0.075 0.10 0.15 0.20 0.25 0.30 0.40 0.50 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10];


%% Coefficients

a1=[2.10 2.43 2.52 2.70 2.86 3.02 3.17 3.44 3.68 4.13 4.42 4.85 5.09 5.58 5.74 5.97 6.35 6.61];
a2=[1.01 0.98 0.91 0.87 0.85 0.85 0.86 0.87 0.88 0.89 0.87 0.85 0.82 0.82 0.82 0.82 0.82 0.82];
a3=[-0.105 -0.0115 0.0985 0.1935 0.2355 0.2515 0.2555 0.238 0.209 0.1455 0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125];
a4=[-0.10 -0.10 -0.10 -0.121 -0.136 -0.147 -0.157 -0.172 -0.183 -0.204 -0.219 -0.240 -0.254 -0.275 -0.275 -0.275 -0.275 -0.275];
phi1=[0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.31 0.34 0.35 0.35 0.35 0.35 0.35];
phi2=[0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28 0.28];
tau1=[0.20 0.20 0.20 0.20 0.20 0.20 0.20 0.20 0.20 0.20 0.20 0.26 0.32 0.39 0.41 0.41 0.41 0.41];
tau2=[0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15];
                  %1  %2  %3  %4   %5   %6  %7     %8   %9
coeff_all= [T_period',a1',a2',a3',a4',phi1',phi2',tau1',tau2'];
   
% Interpolation of coefficients
for j=1:length(T)
    
for i=2:size(coeff_all,2)
    coeff(i)=interp1(log(coeff_all(:,1)),coeff_all(:,i),log(T(j)),'linear','extrap');
end

    if M < 5.25
        phi=coeff(6);
        
    elseif M >=5.25 && M <= 6.5
        
        phi=coeff(6)+ (M-5.25)*(coeff(7)-coeff(6))/1.25;
        
    else
        phi=coeff(7);
    end
    
    if M < 5.25
        tau=coeff(8);
        
    elseif M >=5.25 && M <= 6.5
        
        tau=coeff(8)+ ((M-5.25)*(coeff(9)-coeff(8)))/1.25;
        
    else
        tau=coeff(9);
    end
    
    sigma=sqrt(tau^2+phi^2);
    
    %% sqrt_IP & IP

    ln_sqrtIP(j) = coeff(2)+ coeff(3)*log(targetSa(j))+coeff(4)*(M-6)+coeff(5)*log(Rrup+5*exp(0.4*(M-6)));
    
    sqrtIP(j)= exp(ln_sqrtIP(j)); 
   
    sigma_sqrtIP(j)=sigma;
   
    IP(j)=(sqrtIP(j))^2; 
    
end
