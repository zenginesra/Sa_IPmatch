function [accH1, accH2] = Initialize(accH1,accH2,dt)

% delete the NANs in the array
accH1(all(isnan(accH1),2),:)=[];
accH2(all(isnan(accH2),2),:)=[];

% equalizing vector lengths
lengthData = min(length(accH1), length(accH2));
accH1=accH1(1:lengthData);
accH2=accH2(1:lengthData);

%zero padding at the end (3 sec)
accH1=[accH1;zeros(floor(3/dt),1)];
accH2=[accH2;zeros(floor(3/dt),1)];