function [U,x,t]=simpleDiffusionB(nX,nT,X,T)
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

%Test case: [U,x,t]=simpleDiffusion(100,1000,1,30);

%Default test case
if nargin ==0
    [nX,nT, X, T]=deal(100,100000,1,20);
end

%-------------------------------Parameters--------------------------------%
gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=120;    %sodium ion conductance
ena=112;    %115; %sodium reversal potential
ni=0.5; 
mi=0.4;
hi=0.2;
gl=0.3;     %leak conductance
el=10.6;    %leak potential

a=.01;      %radius
R=40;       %resistance
c=0.9e-2;   %membrane capacitance
alpha=a/(2*R*c); %diffusion coefficient
%-------------------------------------------------------------------------%

global k
k=T/nT; % time step size
h=X/nX; % equidistant spatial discretization width
x=linspace(0,X,nX+1);    % space vector
t=linspace(0,T,nT+1);    % time vector

%-------------------------------Output------------------------------------%
fprintf('Number of spatial steps = %i\n',nX);
fprintf('Spatial dx = %d\n',h);
fprintf('Number of time steps = %i\n',nT);
fprintf('Temporal dt = %d\n',k);
%-------------------------------------------------------------------------%


%solution
V=zeros(nX+1,length(t)); 
N=zeros(nX+1,length(t));
M=zeros(nX+1,length(t));
H=zeros(nX+1,length(t));
U=[V;N;M;H];    %output solution

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
%set boundary conditions at one end of line;
vstart=-55;
U(1,1:ceil(nT/4))=-vstart;
U(1*nR+1:2*nR,1)=ni*ones(nR,1);
U(2*nR+1:3*nR,1)=mi*ones(nR,1);
U(3*nR+1:4*nR,1)=hi*ones(nR,1);

%coefficient for Crank-Nicolson Solve
r=alpha*k/(2*h^2); 
%Matrices for Crank-Nicolson Solve
myI=eye(nX+1);
A=conv2(myI,[-r (1+2*r) -r],'same');
B=conv2(myI,[r (1-2*r) r],'same');
size(myI)

for j=1:length(t)-1
        
        U(:,j+1)=U(:,j)+k.*f(U(:,j));
        Utmp=U(1:nR,j+1);
        U(1:nR,j+1)=A\(B*Utmp+r*myI(:,1)*(g(x(1),t(j),X,T)+g(x(1),t(j+1),X,T))+r*myI(:,nX+1)*(g(x(nX+1),t(j),X,T)+g(x(nX+1),t(j+1),X,T)));
        
        %Utmp=U(:,j+1);
        %U(:,j+1)=Utmp+k.*f(Utmp);
        
        %Utmp=U(:,j+1);
        %Y2=Utmp+(0.5).*k.*f(Utmp);
        %Y3=Utmp+(0.5).*k.*f(Y2);
        %Y4=Utmp+k.*f(Y3);
        %U(:,j+1)=Utmp+(1/6).*k.*(f(Utmp)+2.*f(Y2)+2.*f(Y3)+f(Y4));

end


%ymax=max(max(U(1:nR,:)));
%ymin=min(min(U(1:nR,:)));
%figure(1)
%for j=1:nT+1
% plot(x,U(1:nR,j));
% xlabel('x pos');
% ylabel('u');
% title(sprintf('t = %f seconds',t(j)));
% ylim([ymin,ymax]);
% xlim([0,X]);
% drawnow
%end

figure(2)
[X,T]=meshgrid(x,t);
X=X';
T=T';
surf(X,T,U(1:nR,:));
shading interp;
colormap hot;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');
end

function v=g(X,t,L,T)
vstart=-55;
if (X==0) && t<=T
    %set boundary conditions at one end of line;
    v=-vstart;
else
    v=0;
end

%if X==0 && t<T/2
%    v=100e-3;
%elseif X==0 && t>=T/2
%    v=20e-3;
%else
%    v=0;
%end

%if t>=T/2 && X==0
%    v=-100e-3;
%end

%if t>=T/2 && X==L
%    v=000e-3;
%end
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