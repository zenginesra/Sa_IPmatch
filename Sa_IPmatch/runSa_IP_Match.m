close all; clc;
clear all;
%%-------------------------------------------------------------------------
% This program includes a matching algorithm that enables near-fault records to be adjusted so that they match 
% both a target response spectrum and Instantaneous Power(IP) spectrum. The procedure employs a time-domain
% adjustment approach for the IP.

% Esra Zengin, Norman A. Abrahamson
% <esrazengin@gmail.com>
% Last Updated: 20 March 2021

% Referenced Manuscripts: 
% Zengin E., and Abrahamson, N.A. (2021) A Procedure for Matching the Near-Fault
% Ground Motions based on Spectral Accelerations and Instantaneous Power,
% Earthquake Spectra, in press.

% Abrahamson, N.A. (1992) Non-stationary spectral matching, Seismological
% Research Letters, 63(1), 30.

% INPUT VARIABLES:
% accH1              : H1 component of the acceleration (in g)
% accH2              : H2 component of the acceleration (in g)
% dt_original        : original time step of the record
% M                  : moment magnitude of the target scenario earthquake
% Rrup               : rupture distance of the target scenario earthquake(km)
% targetSpec         : target spectral accelerations (in g)
% Tall               : target periods (s)
% T_lower            : lower period for matching 
% T_upper            : upper period for matching 
% targetPGA          : target peak ground acceleration (in g)
% ksi                : damping ratio
% matched_component  : component used in the matching process

% OUTPUT VARIABLE
% acc_final          : adjusted acceleration time history (in g)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% USER INPUT
% load H1 and H2 components of the acceleration
accH1=load('RSN1114_KOBE_H1.txt');
accH2=load('RSN1114_KOBE_H2.txt');
dt_original=0.01;

% load target response spectrum and periods
specData=load('Target_spectrum.txt');
Tall=specData(:,1);
targetSpec=specData(:,end);

% User Inputs
targetPGA=0.86;
matched_component='H1'; % record component(H1 or H2)
T_upper=5.0; 
T_lower=0.05;
M=7;
Rrup=3; 
%%%
if dt_original>0.01    
accH1=interp(accH1,2);
accH2=interp(accH2,2);
dt=dt_original/2;
else
    dt=dt_original;
end
%%%
% Default parameters
tolerance_Sa=0.05; % tolerance level for Sa misfit
tolerance_IP=0.10; % tolerance level for IP misfit
total_iter=40; % the maximum number of iterations
ksi=0.05;
showPlot=1; % (display=1, no display=0)

%%% User Inputs End Here %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resample target Sa spectrum
PerTgt=logspace(-1.4,log10(T_upper),100); % can be changed
targetall=interp1((Tall),(targetSpec),(PerTgt),'linear','extrap');

PerMatch=PerTgt((PerTgt>T_lower)&(PerTgt<T_upper));
targetMatch=targetall((PerTgt>T_lower)&(PerTgt<T_upper));

%% Target IP spectrum computation

targetSa=targetMatch;
[IP,sqrtIP,sigma_sqrtIP]=IP_GMM(M,Rrup,PerMatch,targetSa); 
targetIP=IP';

%% Initializing
[accH1, accH2]=Initialize(accH1,accH2,dt);
t=(0:dt:length(accH1)*dt-dt)';

% Computation of the original record spectrum
[rsp1,~,~]=CalcRsp(accH1,dt,ksi,PerMatch);
[rsp2,~,~]=CalcRsp(accH2,dt,ksi,PerMatch);

% Initial scaling of records
[acc_mod]=PGA_Scale(accH1,accH2,targetPGA);

% Scaled record spectrum
[rsp1_mod,tpeak1,tindex1]=CalcRsp(acc_mod.H1,dt,ksi,PerMatch);
[rsp2_mod,tpeak2,tindex2]=CalcRsp(acc_mod.H2,dt,ksi,PerMatch);

%% Component-specific target SA
[SaH1,SaH2,RotD50,RotD100]=SaRotDnn(acc_mod.H1,acc_mod.H2,dt,ksi,PerMatch);

% Component Variability Ratio (CVR)
SaH1_ratio=SaH1./RotD50;
SaH2_ratio=SaH2./RotD50;

targetSa_H1=targetMatch.*SaH1_ratio;
targetSa_H2=targetMatch.*SaH2_ratio;
%% Component-specific target IP spectrum

% Velocity time histories (cm/s)
VelH1=cumsum(acc_mod.H1*981)*dt; 
VelH2=cumsum(acc_mod.H2*981)*dt;

[IP_H1,~,~,~]=Calc_IP(VelH1,dt,PerMatch);
[IP_H2,~,~,~]=Calc_IP(VelH2,dt,PerMatch);

% Component Variability Ratio
IP_GM=sqrt(IP_H1.*IP_H2); % geometric-mean IP
IP_H1ratio=IP_H1./IP_GM;
IP_H2ratio=IP_H2./IP_GM;

targetIP_H1=targetIP.*IP_H1ratio';
targetIP_H2=targetIP.*IP_H2ratio';

if (strcmp(matched_component,'H1'))  
    acc=acc_mod.H1;
    targetSa_comp=targetSa_H1;
    targetIP_comp=targetIP_H1;
else
    acc=acc_mod.H2;
    targetSa_comp=targetSa_H2;
    targetIP_comp=targetIP_H2;
end

%% Subgroups
acc=reshape(acc,1,[]);
N_sample=length(PerMatch);
n_elements=20; % number of elements in the matrix
k=ceil(N_sample/n_elements); % total number of subgroup

%%% Sa subgroups
target_Saperiod=cell(length(k));
target_Sa=cell(length(k));
for i=1:k   
target_Saperiod{i}=PerMatch(i:k:end);
target_Sa{i}=targetSa_comp(i:k:end); 
end

%%%IP subgroups
[~,index] = min(abs(PerMatch-0.40)); % minimum period for IP matching

target_IPperiod=cell(length(k));
target_IP=cell(length(k));
for i=1:k    
target_IPperiod{i}=PerMatch(index+i:k:end);
target_IP{i}=targetIP_comp(index+i:k:end); 
end
%% Iterations

for iter=1:total_iter
  
    for i=1:k
    
        [acc]=Sa_Matching(acc,dt,ksi,target_Saperiod{1,i},target_Sa{1,i});

        [acc]=IP_Matching(acc,dt,target_IPperiod{1,i},target_IP{1,i});
        
    end
  
    % Check the misfits
    Vel_matched=cumsum(acc*981).*dt;

    [IP_final,~,~,~]=Calc_IP(Vel_matched,dt,PerMatch(index:end));
    
    logmisfit_IP=log((targetIP_comp(index:end)))-log((IP_final))';
     
    [rsp_matched,~,~,~]=CalcRsp(acc,dt,ksi,PerMatch);
    
    target1=targetSa_comp;
 
    for i=1:length(PerMatch)
        if rsp_matched(i)>=0
            target1(i)=abs(target1(i));
        else
            target1(i)=-abs(target1(i));
        end
        misfit(i)=target1(i)-rsp_matched(i);
        logmisfit_Sa(i)=log((target1(i)))-log((rsp_matched(i)));
        
    end
    
    Sa_misfit=(mean(abs(logmisfit_Sa)));
    IP_misfit=(mean(abs(logmisfit_IP)));
    
    fprintf('Iter %2.0f\n Sa_misfit %.3f\n IP_misfit %.3f\n',iter,Sa_misfit,IP_misfit);
    
%  Termination & saving results

if (mean(abs(logmisfit_Sa))<=tolerance_Sa) && (mean(abs(logmisfit_IP))<=tolerance_IP)
 
    % Baseline correction and tapering
   
    Vel_matched=cumsum(acc*981).*dt;
    Disp_matched=(cumsum(Vel_matched)).*dt;
    y=Disp_matched';% Apply polynomial fit to displacement
    polydeg =5; %degree of polynomial
    skip_pow = [1,0]; % Skip the constant and x term.
    opts = fitoptions( 'Method', 'LinearLeastSquares'); 
    opts.Lower = -inf(1, polydeg + 1);
    opts.Upper = inf(1, polydeg + 1);
    %Skip coefficients have a range from 0 to 0.
    opts.Lower(polydeg + 1 - skip_pow) = 0;
    opts.Upper(polydeg + 1 - skip_pow) = 0;
    [fitresult, gof] = fit( t, y, ['poly', num2str(polydeg)] , opts);
    
    Ydati=fitresult.p1*t.^5+fitresult.p2*t.^4+fitresult.p3*t.^3+fitresult.p4*t.^2;
    Ydati_d=diff(Ydati)./dt; 
    Ydati_d=[0;Ydati_d];
    
    Ydati_dd=diff(Ydati_d)./dt;
    Ydati_dd=[0;Ydati_dd];
    Ydati_dd=Ydati_dd*0.01/9.81; % should be subtracted from acceleration
    
    % Apply tapering (cosine bell)
    te=6; % taper percentage
    npts=length(acc);
    n1=round(npts-(npts*te)/100);
    
    for i=1:length(acc)
        
        if i<=n1
            w(i)=1;
            wPrime(i)=0;
            wPrime2(i)=0;
            
        else
            
            w(i)=0.5*(cos(pi*(i-n1)/(npts-n1))+1);
            
            wPrime(i)=-0.5*pi/((npts-n1)*dt)*sin(pi*(i-n1)/(npts-n1));
            
            wPrime2(i)=-0.5*(pi/((npts-n1)*dt))^2*cos(pi*(i-n1)/(npts-n1));
            
        end
        
        bPrime2=Ydati_dd;
        bPrime=Ydati_d;
        b=Ydati;
        
        % Corrected accelerogram
        acc_final(i)=(acc(i)-bPrime2(i)')*w(i)*981+...
            2.*(Vel_matched(i)-bPrime(i)')*wPrime(i)+(Disp_matched(i)-b(i))*wPrime2(i);
        
    end
    
    acc_final=acc_final*0.01/9.81; %(g)
    Vel_final=cumsum(acc_final*981).*dt; %(cm/s)
    Disp_final=cumsum(Vel_final)*dt; %(cm)
    
    % Save the matched acceleration (g)
    filename='Output_acc.txt';
    fid=fopen(filename,'w+');
    misfit_output='Sa mean misfit is %7.4f and IP mean misfit is %7.4f\n';
    fprintf(fid,misfit_output,Sa_misfit,IP_misfit);
    dlmwrite(filename,[t(:),acc_final(:)],'Delimiter', '\t','-append');
    fclose(fid)
    
    % Show Plots
    if (showPlot)
    figure (1)
    if (strcmp(matched_component,'H1'))
        acc_original=accH1;
        Vel_original=cumsum(acc_original*981)*dt;
        Disp_original=cumsum(Vel_original)*dt;
    else
        acc_original=accH2; 
        Vel_original=cumsum(acc_original*981)*dt;
        Disp_original=cumsum(Vel_original)*dt;
    end
    
    subplot(3,2,1)
    plot(t,acc_original); hold on;
    ylabel('Acceleration (g)');
    xlabel('Time (s)');
    legend('Original');
    
    subplot(3,2,2)
    plot(t,acc_final,'r');
    xlabel('Time (s)');
    ylabel('Acceleration (g)');
    legend('Adjusted');
    hold on;
    
    subplot(3,2,3)
    plot(t,Vel_original); hold on;
    xlabel('Time (s)');
    ylabel('Velocity (cm/s)');
    
    subplot(3,2,4)
    plot(t,Vel_final,'r');
    hold on;
    xlabel('Time (s)');
    ylabel('Velocity (cm/s)');
    
    subplot(3,2,5)
    plot(t,Disp_original,'b');
    xlabel('Time (s)');
    ylabel('Displacement (cm)');
    hold on;
    
    subplot(3,2,6)
    plot(t,Disp_final,'r');
    hold on;
    xlabel('Time (s)');
    ylabel('Displacement (cm)');
    
    % Plot Response Spectra and IP Spectra 
    figure (2)
    subplot(1,2,1)
    loglog(PerMatch,abs(rsp_matched),'-rd','Markersize',1.4);
    hold on;
    loglog(PerMatch,targetSa_comp,'--b')
    hold on;
    loglog(PerMatch,abs(rsp1_mod),'-.k')
    hold on
    loglog(PerMatch,targetMatch,'g','linewidth',1.2)
    hold on
    xlabel('Period(s)','Fontsize',12)
    hold on
    ylabel('Sa (g)','Fontsize',12)
    set(gca,'fontsize',12)
    legend({'Matched Record','Component Target','Linearly Scaled','Target RotD50'},'Fontsize',11,'Location','southwest');
    legend boxoff
    xlim([0.01 10])
    ylim([0.1 10])
    
    subplot(1,2,2)
    
    loglog(PerMatch(index:(end)),(IP_final),'-rd','Markersize',1.4)
    hold on;
    loglog(PerMatch(index:end),(targetIP_comp(index:end)),'--b')
    hold on
    loglog(PerMatch(index:end),(IP_H1(index:end)),'-.k');
    hold on;
    loglog(PerMatch(index:end),targetIP(index:end),'g','linewidth',1.2)
    hold on
    legend({'Matched Record','Component Target','Linearly Scaled','Target IP (Geomean)'},'Fontsize',11,'Location','southwest')
    legend boxoff
    xlabel('Period(s)','Fontsize',12)
    hold on
    ylabel('IP(cm^2/s^2)','Fontsize',12)
    ylim([100 10E4])
    set(gca,'fontsize',12)
    end
    
     return;
    
elseif iter==total_iter
    warning('iteration limit reached'); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%