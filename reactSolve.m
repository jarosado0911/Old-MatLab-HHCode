function [U,NN,MM,HH,x,t]=reactSolve(Vi,Ni,Mi,Hi,nX,nT,x0,xf,t0,tf,mthd, pltOn)
% t0 is the start-time
% tf is the end-time
% nT is the number of time steps
% nX is the number of spatial steps

% default
if (nargin < 2)
    if nargin==1
        [nX,nT, x0,xf, t0,tf,mthd,pltOn]=deal(500,25000,0,1,0,20,1,1);
        fprintf('Using Forward Euler.\n');
    else
        [nX,nT, x0,xf, t0,tf,mthd,pltOn]=deal(500,7500,0,1,0,20,0,1);
        fprintf('Using Runge-Kutta-4.\n');
    end
    ni=0.5; 
    mi=0.4;
    hi=0.2;
    nR=nX+1;
    Vi=zeros(nR,1);
    Vi(1,1)=55;
    Ni=ni*ones(nR,1);
    Mi=mi*ones(nR,1);
    Hi=hi*ones(nR,1);
end

k=abs(tf-t0)/nT;           % time step size
h=abs(xf-x0)/nX;           % equidistant spatial discretization width
x=linspace(x0,xf,nX+1);    % space vector
t=linspace(t0,tf,nT+1);    % time vector

%fprintf('Number of spatial steps = %i\n',nX);
%fprintf('Spatial dx = %d\n',h);
%fprintf('Number of time steps = %i\n',nT);
%fprintf('Temporal dt = %d\n',k);
%----------------------------------------%
% Parameters for AP
gk=36;      % potassium ion conductance
ek=-12;     % potassium reversal potential
gna=120;    % sodium ion conductance
ena=112;    % 115; %sodium reversal potential
gl=0.3;     % leak conductance
el=10.6;    % leak potential
c=0.009; %c=0.9e-2;   % capacitance
%----------------------------------------%

% solution
nR=length(x);
U=zeros(nR,length(t));
NN=zeros(nR,length(t));
MM=zeros(nR,length(t));
HH=zeros(nR,length(t));

f1 = @(nn,mm,hh,v) (-1/c).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*(v-ena)+gl.*(v-el));
f2 = @(nn,v) an(v).*(1-nn)-bn(v).*nn;
f3 = @(mm,v) am(v).*(1-mm)-bm(v).*mm;
f4 = @(hh,v) ah(v).*(1-hh)-bh(v).*hh;

f = @(V) [f1(V(1*nR+1:2*nR),V(2*nR+1:3*nR),V(3*nR+1:4*nR),V(1:nR));...
          f2(V(1*nR+1:2*nR),V(1:nR));...
          f3(V(2*nR+1:3*nR),V(1:nR));...
          f4(V(3*nR+1:4*nR),V(1:nR))];

%Initialize
U(:,1)=Vi;
NN(:,1)=Ni;
MM(:,1)=Mi;
HH(:,1)=Hi;

V=[U;NN;MM;HH];
%use RK4
    for j=1:length(t)-1
        if mthd ==1
            V(:,j+1)=V(:,j)+k.*f(V(:,j));
        else
            Y2=V(:,j)+(0.5).*k.*f(V(:,j));
            Y3=V(:,j)+(0.5).*k.*f(Y2);
            Y4=V(:,j)+k.*f(Y3);
            V(:,j+1)=V(:,j)+(1/6).*k.*(f(V(:,j))+2.*f(Y2)+2.*f(Y3)+f(Y4));
        end
    end    
U=V(1:nR,:);
NN=V(nR+1:2*nR,:);
MM=V(2*nR+1:3*nR,:);
HH=V(3*nR+1:4*nR,:);
    
    if pltOn==1
        figure(1)
        [X,T]=meshgrid(x,t);
        X=X';
        T=T';
        for jj=1:4
            subplot(2,2,jj)
            surf(X,T,V((jj-1)*nR+1:jj*nR,:));
            shading interp;
            colormap hot;
            colorbar
            xlabel('x pos');
            ylabel('time');
            zlabel('intensity');
        end
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