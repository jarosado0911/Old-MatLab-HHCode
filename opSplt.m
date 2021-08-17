function [U,x,t]=opSplt(nX,nT,x0,xf,t0,tf)

%default
if nargin==0
    [nX, nT,x0, xf, t0, tf]=deal(200,100000,0,.5,0,7);
end

k=abs(tf-t0)/nT;         % time step size
h=abs(xf-x0)/nX;         % equidistant spatial discretization width
x=linspace(x0,xf,nX+1);    % space vector
t=linspace(t0,tf,nT+1);    % time vector

fprintf('Number of spatial steps = %i\n',nX);
fprintf('Spatial dx = %d\n',h);
fprintf('Number of time steps = %i\n',nT);
fprintf('Temporal dt = %d\n',k);

nR=length(x);            % size of spatial vector
nTime=nT+1;            % size of time vector

%----------------------------------Parameters-----------------------------%
%Parameters for AP
a=.01;      % radius
R=40;       % resistance

gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=120;    %sodium ion conductance
ena=112;    %115; %sodium reversal potential
ni=0.5; 
mi=0.4;
hi=0.2;
gl=0.3;     %leak conductance
el=10.6;    %leak potential

% compute coefficients
c=.009;%0.005; %c=0.9e-2; %0.009 % capacitance
b=a/(2*R*c); %denom has c    % diffusion coefficient 
%-------------------------------------------------------------------------%

%---For reaction solve-------%
f1 = @(nn,mm,hh,v) -(gk.*(nn.^4).*(v-ek)+gna.*(mm.^3).*hh.*(v-ena)+gl.*(v-el)).*(1/c);
f2 = @(nn,v) an(v).*(1-nn)-bn(v).*nn;
f3 = @(mm,v) am(v).*(1-mm)-bm(v).*mm;
f4 = @(hh,v) ah(v).*(1-hh)-bh(v).*hh;
%-----------------------------%
U=zeros(4*nR,nTime);   % initialize solution

% Initialize one end of rod with fixed voltage
% This is the first row entries of U
vstart=55;
U(1,:)=vstart;

% Now initialize m,n,h entries
U(nR+1,:)=ni;      %initialize n whole row over time
U(2*nR+1,:)=mi;    %initialize m whole row over time
U(3*nR+1,:)=hi;    %initialize h whole row over time

% Now operator splitting
for j=1:nTime-1
    
    %Now Solve Reaction Equation using forward euler
         U(2:nR-1,j+1)=U(2:nR-1,j)+0.5*k*f1(U(nR+2:2*nR-1,j),U(2*nR+2:3*nR-1,j),U(3*nR+2:4*nR-1,j),U(2:nR-1,j));
      U(nR+1:2*nR,j+1)=U(nR+1:2*nR,j)+0.5*k*f2(U(nR+1:2*nR,j),U(1:nR,j));
    U(2*nR+1:3*nR,j+1)=U(2*nR+1:3*nR,j)+0.5*k*f3(U(2*nR+1:3*nR,j),U(1:nR,j));
    U(3*nR+1:4*nR,j+1)=U(3*nR+1:4*nR,j)+0.5*k*f4(U(3*nR+1:4*nR,j),U(1:nR,j));
   
    % Diffusion Solve
   [Uout,~,~]=diffSolve(U(1:nR,j+1),U(1,j+1),U(nR,j+1),b,nX,3,x0,xf,t(j),t(j+1));
   tmpU=Uout(:,end);
   
    %Now Solve Reaction Equation using forward euler
          U(2:nR-1,j+1)=tmpU(2:nR-1)+0.5*k*f1(U(nR+2:2*nR-1,j),U(2*nR+2:3*nR-1,j),U(3*nR+2:4*nR-1,j),tmpU(2:nR-1));
       U(nR+1:2*nR,j+1)=U(nR+1:2*nR,j)+0.5*k*f2(U(nR+1:2*nR,j),tmpU);
     U(2*nR+1:3*nR,j+1)=U(2*nR+1:3*nR,j)+0.5*k*f3(U(2*nR+1:3*nR,j),tmpU);
     U(3*nR+1:4*nR,j+1)=U(3*nR+1:4*nR,j)+0.5*k*f4(U(3*nR+1:4*nR,j),tmpU);
end

figure(1)
[X,T]=meshgrid(x,t);
X=X';
T=T';
surf(X,T,U(1:nR,:));
view(2);
shading interp;
colormap jet;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');

% figure(2)
% surf(X,T,U(nR+1:2*nR,:));
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');
% 
% figure(3)
% surf(X,T,U(2*nR+1:3*nR,:));
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');
% 
% figure(4)
% surf(X,T,U(3*nR+1:4*nR,:));
% shading interp;
% colormap jet;
% colorbar
% xlabel('x pos');
% ylabel('time');
% zlabel('intensity');

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
