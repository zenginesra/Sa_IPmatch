
function [acc_mod] = PGA_Scale(accH1,accH2,targetPGA)
% This function scales record to the target peak ground acceleration (PGA).

pgaH1=max(abs(accH1));
pgaH2=max(abs(accH2));

pga_geo=sqrt(pgaH1*pgaH2);

pgaH1_ratio=pgaH1/pga_geo;
pgaH2_ratio=pgaH2/pga_geo;

targetH1=targetPGA*pgaH1_ratio;
targetH2=targetPGA*pgaH2_ratio;

scaleH1=targetH1/pgaH1;
scaleH2=targetH2/pgaH2;

acc_mod.H1=accH1.*scaleH1;
acc_mod.H2=accH2.*scaleH2;


