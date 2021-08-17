function [U,x,t]=simpleReact(nX,nT,X,T,mthd)
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

if nargin ==0
    [nX,nT, X, T,mthd]=deal(100,9000,1,20,1);
end

global k
k=T/nT; % time step size
%h=X/nX; % equidistant spatial discretization width
x=linspace(0,X,nX+1);    % space vector
t=linspace(0,T,nT+1);    % time vector

%----------------------------------------%
%Parameters for AP
%a=.01;      %radius
%R=40;       %resistance
gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=120;    %sodium ion conductance
ena=112;    %115; %sodium reversal potential
ni=0.5; 
mi=0.4;
hi=0.2;
gl=0.3;     %leak conductance
el=10.6;    %leak potential

%compute coefficients
c=0.9e-2;
%----------------------------------------%

% solution
V=zeros(nX+1,length(t)); 
N=zeros(nX+1,length(t));
M=zeros(nX+1,length(t));
H=zeros(nX+1,length(t));

U=[V;N;M;H];
nR=length(x);

f1 = @(nn,mm,hh,v) (-1/c).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*(v-ena)+gl.*(v-el));
f2 = @(nn,v) an(v).*(1-nn)-bn(v).*nn;
f3 = @(mm,v) am(v).*(1-mm)-bm(v).*mm;
f4 = @(hh,v) ah(v).*(1-hh)-bh(v).*hh;

f = @(V) [f1(V(1*nR+1:2*nR),V(2*nR+1:3*nR),V(3*nR+1:4*nR),V(1:nR));...
                        f2(V(1*nR+1:2*nR),V(1:nR));...
                        f3(V(2*nR+1:3*nR),V(1:nR));...
                        f4(V(3*nR+1:4*nR),V(1:nR))];
%Initialize
U(1,:)=-55;
U(1*nR+1:2*nR,1)=ni*ones(nR,1);
U(2*nR+1:3*nR,1)=mi*ones(nR,1);
U(3*nR+1:4*nR,1)=hi*ones(nR,1);

if mthd==1 %use RK4
    for j=1:length(t)-1
        %Y1=U(:,j);
        Y2=U(:,j)+(0.5).*k.*f(U(:,j));
        Y3=U(:,j)+(0.5).*k.*f(Y2);
        Y4=U(:,j)+k.*f(Y3);
        U(:,j+1)=U(:,j)+(1/6).*k.*(f(U(:,j))+2.*f(Y2)+2.*f(Y3)+f(Y4));
    end
else %use forward Euler
    for j=1:length(t)-1
        U(:,j+1)=U(:,j)+k.*f(U(:,j));
    end 
    
end

figure(1)
[X,T]=meshgrid(x,t);
X=X';
T=T';
subplot(2,2,1)
surf(X,T,U(1:nR,:));
shading interp;
colormap hot;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');

subplot(2,2,2)
surf(X,T,U(nR+1:2*nR,:));
shading interp;
colormap hot;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');

subplot(2,2,3)
surf(X,T,U(2*nR+1:3*nR,:));
shading interp;
colormap hot;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');

subplot(2,2,4)
surf(X,T,U(3*nR+1:4*nR,:));
shading interp;
colormap hot;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');
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


