function BE_method(nX,nT,X,T)
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

x=linspace(0,X,nX);   % space vector
t=linspace(0,T,nT);   % time vector
k=t(2)-t(1);          % time step size
h=x(2)-x(1);          % equidistant spatial discretization width

%-------------------------------------------------------------------------%
%Biological Parameters for AP
a=.001;     %radius
R=100;       %resistance
gk=36;      %potassium ion conductance
ek=-12;     %potassium reversal potential
gna=120;    %sodium ion conductance
ena=112;    %115; %sodium reversal potential
ni=0.5; 
mi=0.4;
hi=0.2;
gl=0.3;     %leak conductance
el=0.6;     %leak potential


%compute coefficients
b=a/(2*R);
c=0.09;
vstart = 55;
%-------------------------------------------------------------------------%
%Initialize Solutions
u=zeros(length(x),1); u(1)=vstart; 
nn=zeros(length(x),1); nn(1)=ni;
mm=zeros(length(x),1); mm(1)=mi;
hh=zeros(length(x),1); hh(1)=hi;

%define stencil matrix
sten = [1 -2 1];
sysMat = spdiags(ones(nX,1)*sten,-1:1,nX,nX)*b/h^2; %sparse diagonal

%reaction term
f = @(nn,mm,hh,v) (-1).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*...
    (v-ena)+gl.*(v-el));

figure(1)
for i=1:nT
    fprintf(sprintf("Time step %i, t = %f\n",i, i*k))
    u(1)=vstart;
    
    A =eye(length(x))-(k/c).*sysMat;
    u = linsolve(A,u);
    u = u+(k/c).*f(nn,mm,hh,u);
    
    %Forward Euler for time dependent ODEs on n,m,h
    nn=nn+k*(an(u).*(1-nn)-bn(u).*nn);
    mm=mm+k*(am(u).*(1-mm)-bm(u).*mm);
    hh=hh+k*(ah(u).*(1-hh)-bh(u).*hh);
    subplot(2,1,1)
    plot(x,u)
    ylim([-20 120])
    xlim([x(1) x(end)])
    
    subplot(2,1,2)
    scatter(x,0.*x,160,u,'filled','s');
    xlim([x(1) x(end)])
    axis off
    colormap jet
    caxis([-15 100])
    %colorbar
    drawnow
    
end

end

%% Gating functions
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
