function [a,R,gk,ek,gna,ena,gl,el,ni,mi,hi,c,vstart,b]=set_bio_params()
a=.001;     %radius
R=10;       %resistance
gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=153;    %sodium ion conductance
ena=112;    %115; %sodium reversal potential
ni=0.5; 
mi=0.4;
hi=0.2;
gl=0.3;     %leak conductance
el=0.6;     %leak potential
c=0.09;
vstart = 55;

%compute diff coefficient
b=a/(2*R*c);
end

