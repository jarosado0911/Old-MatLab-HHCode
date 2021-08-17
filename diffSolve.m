function [Uout,x,t]=diffSolve(U,bc1,bc2,a,nX,nT,x0, xf, t0, tf)
% Uses Forward Euler and Center Difference
% bc1, bc2 = are the boundary conditions
% a = is the diffusion constant
% nX = number of spatial steps
% nT = number of time steps
% x0, xf = initial and last spatial point
% t0, tf = initial and last time point

k=abs(tf-t0)/nT;         % time step size
h=abs(xf-x0)/nX;         % equidistant spatial discretization width
x=linspace(x0,xf,nX+1);    % space vector
t=linspace(t0,tf,nT+1);    % time vector

nR=length(x);            % size of spatial vector
nTime=nT+1;            % size of time vector

Uout=zeros(nR,nTime);  % initialize solution output
Uout(:,1)=U;           % first column is initial input U

% initialize first and last rows with b.c.
Uout(1,:)=bc1;         
Uout(end,:)=bc2;

% for tridiagonal matrix
r=a*k/(h^2);
% tridiagonal matri
K1D=spdiags(ones(nR-2,1)*[r 1-2*r r],-1:1,nR-2,nR-2);
% used for booking keeping of dirichelet conditions
myI=eye(nR-2);
for j=1:nTime-1
    Uout(2:nR-1,j+1)=K1D*Uout(2:nR-1,j)+r*(myI(:,1)*Uout(1,j+1)+myI(:,end)*Uout(end,j+1));
end

end

