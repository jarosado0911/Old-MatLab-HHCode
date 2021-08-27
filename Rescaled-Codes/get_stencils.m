function [sysMat,A,B]=get_stencils(c,h,k,sz,mthd)
%%
% c is the diffusion coefficient
% h is the spatial width
% k is the time step size
% sz is the size of the 1D spatial x vector
% mthd is the temporal numerical method:
% 0 = Using Forward Euler time stepping
% 1 = Using Backward Euler time stepping
% 2 = Using Crank-Nicolson time stepping

%% initialize basic stencil for spatial discretization, here we are using
% center differencing
sten = [1 -2 1];
sysMat = spdiags(ones(sz,1)*sten,-1:1,sz,sz)*c/h^2; %sparse diagonal

%% FE mthd
if mthd ==0
   A = sysMat; B = eye(sz); 
end

%% BE mthd
if mthd ==1
    A =eye(sz)-k.*sysMat; B=eye(sz);
end

%% CN mthd
if mthd ==2
    A =eye(sz)-(k/2).*sysMat;
    B =eye(sz)+(k/2).*sysMat;
end


end

