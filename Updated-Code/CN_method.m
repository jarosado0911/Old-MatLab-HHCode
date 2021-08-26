function CN_method(nX,nT,X,T)
%% 
% nX = number of spatial steps along a rod
% nT = number of time steps
% X = Length of rod
% T = end time

if nargin == 0
    [nX,nT, X, T,mthd]=set_sim_params();
end

x=linspace(0,X,nX);   % space vector
t=linspace(0,T,nT);   % time vector
k=t(2)-t(1);          % time step size
h=x(2)-x(1);          % equidistant spatial discretization width

%% Biological Parameters for AP
[~,~,gk,ek,gna,ena,gl,el,ni,mi,hi,c,vstart,b]=set_bio_params();

%% Initialize Solutions
u=zeros(length(x),1); u(1)=vstart; 
nn=zeros(length(x),1); nn(1)=ni;
mm=zeros(length(x),1); mm(1)=mi;
hh=zeros(length(x),1); hh(1)=hi;

%% define stencil matrices
[~,A,B]=get_stencils(b,h,k,length(x),2);

%% reaction term
f = @(nn,mm,hh,v) (-1/c).*(gk.*nn.^4.*(v-ek)+gna.*mm.^3.*hh.*...
    (v-ena)+gl.*(v-el));

%% for ODE solves
fn = @(n,v) gating_functions.an(v).*(1-n)-gating_functions.bn(v).*n;
fm = @(m,v) gating_functions.am(v).*(1-m)-gating_functions.bm(v).*m;
fh = @(h,v) gating_functions.ah(v).*(1-h)-gating_functions.bh(v).*h;

%% for loops for solving
figure(1)

movie_name = 'CN_method';
vidfile = VideoWriter(sprintf('%s.mp4',movie_name),'MPEG-4');
open(vidfile);

    for i=1:nT
        %% Operator splitting: first solve diffusion using BE, then solve ODEs
        u = linsolve(A,B*u);
                
        %% Forward Euler for time dependent ODEs on n,m,h
        if mthd == 0
            u = sim_functions.FE(u,f(nn,mm,hh,u),k);
            u(1)=vstart; u(end)=0.0;
        
            nn = sim_functions.FE(nn,fn(nn,u),k);
            mm = sim_functions.FE(mm,fm(mm,u),k); 
            hh = sim_functions.FE(hh,fh(hh,u),k);
        %% RK4 for time dependent ODEs on n,m,h
        else
            u = sim_functions.RK4_react_diff(nn,mm,hh,u,k,f,0);
            u(1)=vstart; u(end)=0.0;
            
            nn = sim_functions.RK4(nn,u,k,fn); 
            mm = sim_functions.RK4(mm,u,k,fm); 
            hh = sim_functions.RK4(hh,u,k,fh);
        end
        %% plotting
        if mod(i,40) == 0
            sim_functions.make_plot(x,u,i*k);
            thisFrame = getframe(gcf);
            writeVideo(vidfile, thisFrame);
        end
    end
    
    close(vidfile);
end