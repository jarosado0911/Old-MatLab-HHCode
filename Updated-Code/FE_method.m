function FE_method(nX,nT,X,T)
%% 
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

if nargin == 0
    [nX,nT, X, T]=deal(200,16000,1,60);
end

x=linspace(0,X,nX);   % space vector
t=linspace(0,T,nT);   % time vector
k=t(2)-t(1);          % time step size
h=x(2)-x(1);          % equidistant spatial discretization width

%% Biological Parameters for AP
[a,R,gk,ek,gna,ena,gl,el,ni,mi,hi,c,vstart]=set_params();
%compute diff coefficient
b=a/(2*R*c);

%% Initialize Solutions
u=zeros(length(x),1); u(1)=vstart; 
nn=zeros(length(x),1); nn(1)=ni;
mm=zeros(length(x),1); mm(1)=mi;
hh=zeros(length(x),1); hh(1)=hi;

%% get stencil matrix
[~,A,~]=get_stencils(b,h,k,length(x),0);

%% reaction term
f = @(nn,mm,hh,v) (-1/c).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*...
    (v-ena)+gl.*(v-el));

%% for ODE solves
fn = @(n,v) an(v).*(1-n)-bn(v).*n;
fm = @(m,v) am(v).*(1-m)-bm(v).*m;
fh = @(h,v) ah(v).*(1-h)-bh(v).*h;

%% for loops for solving
figure(1)

movie_name = 'FE_method';
vidfile = VideoWriter(sprintf('%s.mp4',movie_name),'MPEG-4');
open(vidfile);

    for i=1:nT
        u(1)=vstart; u(end)=0.0;

        %% Forward Euler for Diffusion and Reaction
        u = FE(u,A*u+f(nn,mm,hh,u),k);
        
        %% Forward Euler for time dependent ODEs on n,m,h
        nn = FE(nn,fn(nn,u),k);
        mm = FE(mm,fm(mm,u),k);
        hh = FE(hh,fh(hh,u),k);

        %% plotting
        clf()
        if mod(i,40) == 0
            hold on
            plot(x,u,'k','Linewidth',1.5)
            ylim([-20 120])
            xlim([x(1) x(end)])
            title(sprintf('time = %0.2f [ms]',i*k))
            
            scatter(x,0.*x-15,80,u,'filled','s');
            xlim([x(1) x(end)])
            axis off
            colormap jet
            caxis([-15 100])
            colorbar('southoutside')
            drawnow
            thisFrame = getframe(gcf);
            writeVideo(vidfile, thisFrame);
        end
    end
    
    close(vidfile);
end

%% Forward Euler Step Function
function out = FE(X1,X2,K)
    out = X1 + K.*X2;
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
