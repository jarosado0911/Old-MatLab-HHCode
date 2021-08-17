function [U,x,t]=simpleDiffusion(nX,nT,X,T)
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

%Test case: [U,x,t]=simpleDiffusion(100,1000,1,30);
if nargin ==0
    [nX,nT, X, T]=deal(100,1000,1,30);
end

a=0.015;
global k
k=T/nT; % time step size
h=X/nX; % equidistant spatial discretization width
x=linspace(0,X,nX+1);    % space vector
t=linspace(0,T,nT+1);    % time vector
U=zeros(nX+1,length(t)); % solution

r=a*k/(2*h^2);

A=conv2(eye(nX+1),[-r (1+2*r) -r],'same');
B=conv2(eye(nX+1),[r (1-2*r) r],'same');

myI=eye(nX+1);
for j=1:length(t)-1
U(:,j+1)=A\(B*U(:,j)+r*myI(:,1)*(g(x(1),t(j),X,T)+g(x(1),t(j+1),X,T))+r*myI(:,nX+1)*(g(x(nX+1),t(j),X,T)+g(x(nX+1),t(j+1),X,T)));
%U(:,j+1)=A\(B*U(:,j)+r*myI(:,1)*(g(x(1),t(j),X,T))+r*myI(:,nX+1)*(g(x(nX+1),t(j),X,T)));
end

ymax=max(max(U));
ymin=min(min(U));
figure(1)
for j=1:nT+1
 plot(x,U(:,j));
 xlabel('x pos');
 ylabel('u');
 title(sprintf('t = %f seconds',t(j)));
 ylim([ymin,ymax]);
 xlim([0,X]);
 drawnow
end

figure(2)
[X,T]=meshgrid(x,t);
X=X';
T=T';
surf(X,T,U);
shading interp;
colormap hot;
colorbar
xlabel('x pos');
ylabel('time');
zlabel('intensity');
end

function v=g(X,t,L,T)
%global k
%scal=0.1;

if (X==0 ||X==L) && t<=T/4
    %v=50e-3*sin(t/15)*exp(-scal*k*t);
    v=70e-3;
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

