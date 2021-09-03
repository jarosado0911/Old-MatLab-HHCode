function [u,mm,nn,hh,x,t]=branchAP_Matrix()

%% assign simulation parameters
if nargin == 0
    [nX,nT, X, T,~]=set_sim_params();
    end_loc = floor(nX*2/3)-1;
    brc_ngbr_loc = end_loc+1;
    brc_loc = floor(nX/3);
    src_loc=1;
end

%%
x=linspace(0,X,nX);    % space vector
t=linspace(0,T,nT);    % time vector
k = t(2)-t(1);         % time step size
h = x(1)-x(2);         % spatial step size

%% Biological Parameters for AP
[~,~,gk,ek,gna,ena,gl,el,ni,mi,hi,c,vstart,b]=set_bio_params();

%% Initialize Solutions
u=zeros(length(x),1); u(src_loc)=vstart; 
nn=zeros(length(x),1); nn(1)=ni;
mm=zeros(length(x),1); mm(1)=mi;
hh=zeros(length(x),1); hh(1)=hi;

%% get stencil matrix
[~,A,~]=get_stencils(b,h,k,length(x),0);

% this takes care of the entries at the one branch point
entry1 = A(brc_loc,brc_loc-1);
A(brc_loc,brc_ngbr_loc) = entry1;
A(brc_loc,brc_loc) = A(brc_loc,brc_loc)*3/2;

A(brc_ngbr_loc,end_loc)=0;
A(brc_ngbr_loc,brc_loc)=entry1;

figure(1)
set(gcf, 'Position',  [100, 100, 2500, 600]);
subplot(1,2,1)
spy(A)
title('Sparsity Pattern')

%% reaction term
f = @(nn,mm,hh,v) (-1/c).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*...
    (v-ena)+gl.*(v-el));

%% for ODE solves
fn = @(n,v) gating_functions.an(v).*(1-n)-gating_functions.bn(v).*n;
fm = @(m,v) gating_functions.am(v).*(1-m)-gating_functions.bm(v).*m;
fh = @(h,v) gating_functions.ah(v).*(1-h)-gating_functions.bh(v).*h;

%% For plotting the y-branch geometry in 2d
xpts = [1:1:end_loc];
ypts = xpts.*0;

xpts_brch = [1:1:nX - end_loc]/sqrt(2)+brc_loc;
ypts_brch = [1:1:length(xpts_brch)]/2;

xpts = [xpts,xpts_brch];
ypts = [ypts,ypts_brch];

s2 = subplot(1,2,2);
%% Solver loop
for i=1:nT
    %% Forward Euler for Diffusion and Reaction (Lie-trotter opspilt)
    u = sim_functions.FE(u,A*u,k);
    u = sim_functions.FE(u,f(nn,mm,hh,u),k);
    
    %% I use dirichlet conditions here
    u(src_loc)=vstart; u(end_loc)=0; u(end)=0;

    %% Forward Euler for ODEs
    nn = sim_functions.FE(nn,fn(nn,u),k);
    mm = sim_functions.FE(mm,fm(mm,u),k); 
    hh = sim_functions.FE(hh,fh(hh,u),k);
    
    if mod(i,80) == 0
        cla(s2)
        scatter(xpts,ypts,80,u,'filled','o')
        xlim([min(xpts),max(xpts)])
        ylim([min(ypts),max(ypts)])
        title(sprintf('time = %0.2f [ms]',i*k))
        axis off
        colormap jet
        caxis([-15 100])
        colorbar('southoutside')
        drawnow
    end
end

end

