function [Uout,x,t]=rctSolve(NN,MM,HH,U,nX,frc,nT,x0, xf, t0, tf)
% c = for spilt time steps
% f = function
% U = initial condition
% nX = number of spatial steps
% nT = number of time steps
% x0, xf = initial and last spatial point
% t0, tf = initial and last time point

gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=120;    %sodium ion conductance
ena=112;    %115; %sodium reversal potential
gl=0.3;     %leak conductance
el=10.6;    %leak potential

c=.009;
k=abs(tf-t0)/nT;         % time step size
x=linspace(x0,xf,nX+1);    % space vector
t=linspace(t0,tf,nT+1);    % time vector

nR=length(x);            % size of spatial vector
nTime=nT+1;            % size of time vector

Uout=zeros(4*nR,nTime);  % initialize solution output
Uout(1:nR,1)=U;           % first column is initial input U
Uout(nR+1:2*nR,1)=NN;
Uout(2*nR+1:3*nR,1)=MM;
Uout(3*nR+1:4*nR,1)=HH;

%---For reaction solve-------%
f1 = @(nn,mm,hh,v) (-1/c).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*(v-ena)+gl.*(v-el)); %(-1/c).*
f2 = @(nn,v) an(v).*(1-nn)-bn(v).*nn;
f3 = @(mm,v) am(v).*(1-mm)-bm(v).*mm;
f4 = @(hh,v) ah(v).*(1-hh)-bh(v).*hh;

% f = @(V) [f1(V(1*nR+2:2*nR-1),V(2*nR+2:3*nR-1),V(3*nR+2:4*nR-1),V(2:nR-1));...
%           f2(V(1*nR+1:2*nR),V(1:nR));...
%           f3(V(2*nR+1:3*nR),V(1:nR));...
%           f4(V(3*nR+1:4*nR),V(1:nR))];
%-----------------------------%

for j=1:nTime-1
    Uout(2:nR-1,j+1)=Uout(2:nR-1,j)+frc.*k.*f1(Uout(nR+2:2*nR-1,j),...
                                               Uout(2*nR+2:3*nR-1,j),...
                                               Uout(3*nR+2:4*nR-1,j),...
                                               Uout(2:nR-1,j));
                                           
    Uout(nR+1:2*nR,j+1)=Uout(nR+1:2*nR,j)+frc.*k.*f2(Uout(nR+1:2*nR,j),Uout(1:nR,j));
    Uout(2*nR+1:3*nR,j+1)=Uout(2*nR+1:3*nR,j)+frc.*k.*f3(Uout(2*nR+1:3*nR,j),Uout(1:nR,j));
    Uout(3*nR+1:4*nR,j+1)=Uout(3*nR+1:4*nR,j)+frc.*k.*f4(Uout(3*nR+1:4*nR,j),Uout(1:nR,j));
end

end

function out=an(vin)
out=(0.01).*(10-vin)./(exp((10-vin)./10)-1);
end

function out=bn(vin)
out=(0.125).*exp(-vin./80);
end

function out=am(vin)
out=(0.1).*(25-vin)./(exp((25-vin)./10)-1);
end

function out=bm(vin)
out=4.*exp(-vin./18);
end

function out=ah(vin)
out=(0.07).*exp(-vin./20);
end

function out=bh(vin)
out=1./(exp((30-vin)./10)+1);
end

