close all
clc
clear
format short

addpath(genpath(pwd))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% | |
% |  \
% |_|\|
%
% k = Pi/3; 
% X = {Cos[t], Sin[t], s} ; 
% r = 1/2;
% X2 = {r Cos[u], r Sin[u], v} ;
% M = {    { Cos[k], 0, -Sin[k]},
%          { 0,      1,       0},
%          { Sin[k], 0,  Cos[k]}  };
% M*X2'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------------------parameters-----------------
R = 1;    % radius of big cylinder
k = pi/3; % rotation angle
r = 2/3;  % radius of thin cylinder
Ni = 2;   % number of points from angle -alpha to alpha
%ub = 5;   % z=ub for top circle 
%lb = -5;  % z=lb for bottom circle and side circle
addN = 60; % add to Y shape

%%

alpha = asin(r/R);

tt = linspace(-alpha,alpha,Ni+1);
t  = [tt(1),2/3*tt(1)+1/3*tt(2),tt(2:end-1),1/3*tt(end-1)+2/3*tt(end),tt(end)];
Ni = Ni+2;
u1 = asin(R*sin(t)/r);
v1 = -(R*cos(t)-r*cos(k)*cos(u1))/sin(k);
s1 = v1*cos(k)+r*cos(u1)*sin(k);

%t = linspace(-alpha,alpha,Ni+1);
u2 = pi-asin(R*sin(t)/r);
v2 = -(R*cos(t)-r*cos(k)*cos(u2))/sin(k);
s2 = v2*cos(k)+r*cos(u2)*sin(k);

Xi = [cos(t)',sin(t)',s1'];  % upper side of the interface curve 
dh = 0.3;
T1 = max(Xi(:,3))+dh;
ub = T1;
X_top0 = [cos(t)',sin(t)',(t*0+ub)'];

X1i = [cos(t)',sin(t)',s2']; % down side of the interface curve 
T2 = min(X1i(:,3));
lb = T2-dh;

v = lb;
Yi = [(-v*sin(k)+r*cos(k)*cos(u1))',(r*sin(u1))',(v*cos(k)+r*cos(u1)*sin(k))'];  %side circle up
Y1i = [(-v*sin(k)+r*cos(k)*cos(u2))',(r*sin(u2))',(v*cos(k)+r*cos(u2)*sin(k))']; %side circle down
%% 
midv = Xi(1,3);
t = linspace(alpha,pi,floor(Ni/2*(pi-alpha)/alpha)+1);
X_top1 = [cos(t)',sin(t)',(t*0+ub)'];
X_mid1 = [cos(t)',sin(t)',(t*0+midv)'];
%%
t = linspace(-pi,-alpha,floor(Ni/2*(pi-alpha)/alpha)+1);
X_top2 = [cos(t)',sin(t)',(t*0+ub)'];
X_mid2 = [cos(t)',sin(t)',(t*0+midv)'];

%% Vl is the list of vertices 
ind = 0;
% ind=ind+1;
% Vl(ind,:) = [ind,0,0,ub]; % top center 1
% ind=ind+1;
% Vl(ind,:) = [ind,0,0,lb]; % bottom center 2
% ind=ind+1;
% Vl(ind,:) = [ind,-v*sin(k),0,v*cos(k)]; % side center 3
Vl(1,:) = [1,0,0,0];
marker3D(1,:) = [0,0]; %[index of Vertices, marker]
swc1D = zeros(addN*3+4,8); %[index of Vertices,type,x,y,z,radius, ParentID,marker]
%% I_side is the interface curve on side (points marker is 0)
swc1D(addN+2,:) =[addN+2,2,0,0,0,R,addN+1,0];

I_side = zeros(2*Ni,4);
for i=1:Ni
    ind = ind+1;
    I_side(i,:) = [ind,Xi(i,1),Xi(i,2),Xi(i,3)]';
    Vl(ind,:) = [ind,Xi(i,1),Xi(i,2),Xi(i,3)]';
    marker3D(ind,:) = [ind,0];
end
for i=1:Ni
    ind = ind+1;
    I_side(Ni+i,:) = [ind,X1i(Ni+2-i,1),X1i(Ni+2-i,2),X1i(Ni+2-i,3)]';
    Vl(ind,:) = I_side(Ni+i,:);
    marker3D(ind,:) = [ind,0];
end
%% I_top is the top of interface curve connected with mid circle (points on mid circle is 0)
it1 = size(Xi,1);
it2 = size(X_mid1,1)-2;
it3 = size(X_mid2,1)-1;
Nit = it1+it2+it3;
I_top = zeros(Nit,4);
 for i=1:it1
     I_top(i,:) = I_side(i,:);
 end
 for i=1:it2
      ind = ind+1;
      I_top(it1+i,:) = [ind,X_mid1(i+1,1),X_mid1(i+1,2),X_mid1(i+1,3)]';
      Vl(ind,:)= I_top(it1+i,:);
      marker3D(ind,:) = [ind,0];
 end
 for i=1:it3
      ind = ind+1;
      I_top(it1+it2+i,:) = [ind,X_mid2(i,1),X_mid2(i,2),X_mid2(i,3)]';
      Vl(ind,:) = I_top(it1+it2+i,:);
      marker3D(ind,:) = [ind,0];
 end
%% I_down is the lower side of interface curve connected with mid circle 
I_down = zeros(Nit,4);
I_down(1,:) =I_side(1,:);
I_down(2:it1-1,:)=I_side(end:-1:it1+1,:);
I_down(it1,:) =I_side(it1,:);
I_down(it1+1:Nit,:) = I_top(it1+1:end,:);
%% swc file

%% C_top is the vertical top circle (points on top circle with marker 1)
swc1D(addN+1,:) =[addN+1,2,0,0,ub,R,addN,1];
C_top = zeros(Nit,4);
 for i=1:it1
     ind=ind+1;
     C_top(i,:) = [ind,X_top0(i,1),X_top0(i,2),X_top0(i,3)]';
     Vl(ind,:)  = C_top(i,:); 
     marker3D(ind,:) = [ind,1];
 end
 for i=1:it2
     ind = ind+1;
     C_top(it1+i,:) = [ind,X_top1(i+1,1),X_top1(i+1,2),X_top1(i+1,3)]';
     Vl(ind,:)  = C_top(it1+i,:); 
     marker3D(ind,:) = [ind,1];
 end
 for i=1:it3
     ind = ind+1;
     C_top(it1+it2+i,:) = [ind,X_top2(i,1),X_top2(i,2),X_top2(i,3)]';
     Vl(ind,:)  = C_top(it1+it2+i,:); 
     marker3D(ind,:) = [ind,1];
 end
 
 
%% C_bot is the vertical bottom circle (points on bottom circle with marker -1)
swc1D(addN+3,:) =[addN+3,2,0,0,lb,R,addN+2,-1];
C_bot = zeros(Nit,4);
 for i=1:Nit
     ind=ind+1;
     C_bot(i,:) = [ind,C_top(i,2),C_top(i,3),C_top(i,4)-ub+lb]';
     Vl(ind,:)  = C_bot(i,:); 
     marker3D(ind,:) = [ind,-1];
 end
 
%% C_side is the side circle (on the small cylinder) marker -10000
itm = 2*(addN+1)+1+1;
swc1D(itm,:) =[itm,2,-v*sin(k),0,v*cos(k),r,addN+2,-10000];
C_side = zeros(2*Ni,4);
for i=1:Ni
    ind = ind+1;
    C_side(i,:) = [ind,Yi(i,1),Yi(i,2),Yi(i,3)]';
    Vl(ind,:)  = C_side(i,:);
    marker3D(ind,:) = [ind,-10000];
end
for i=1:Ni
    ind = ind+1;
    C_side(Ni+i,:) = [ind,Y1i(Ni+2-i,1),Y1i(Ni+2-i,2),Y1i(Ni+2-i,3)]';
    Vl(ind,:)  = C_side(Ni+i,:);
    marker3D(ind,:) = [ind,-10000];
end

%% Plot the dots for Y shape geometry

% scatter3(I_side(:,2),I_side(:,3),I_side(:,4),'filled')
% hold on
% scatter3(C_side(:,2),C_side(:,3),C_side(:,4),'filled')
% hold on
% s=size(I_side,1);
% for i=1:s
%     pts = [I_side(i,2:end); C_side(i,2:end)];
%     plot3(pts(:,1), pts(:,2), pts(:,3))
%     hold on
% end
% s=size(I_top,1);
% for i=1:s
%     pts = [I_top(i,2:end); C_top(i,2:end)];
%     plot3(pts(:,1), pts(:,2), pts(:,3))
%     hold on
% end
% s=size(I_top,1);
% for i=1:s
%     pts = [I_down(i,2:end); C_bot(i,2:end)];
%     plot3(pts(:,1), pts(:,2), pts(:,3))
%     hold on
% end
% scatter3(X_top0(:,1),X_top0(:,2),X_top0(:,3),'filled')
% hold on
% scatter3(X_top1(:,1),X_top1(:,2),X_top1(:,3),'filled')
% hold on
% scatter3(X_top2(:,1),X_top2(:,2),X_top2(:,3),'filled')
% hold on
% scatter3(C_bot(:,2),C_bot(:,3),C_bot(:,4),'filled')
% hold on
% scatter3(X_mid1(:,1),X_mid1(:,2),X_mid1(:,3),'filled')
% hold on
% scatter3(X_mid2(:,1),X_mid2(:,2),X_mid2(:,3),'filled')
% grid on
% daspect([1 1 1])
% xlabel('X')
% ylabel('Y')
% zlabel('Z')


%% Nodes first
%% top cylinder
%collect all triangles
jt = 0;
Triset(1,:) = [0,0,0];
Nt = size(C_top,1);
for j=1:Nt
%     fprintf(fileID,'%s %d %d %d\n','f',C_top(j,1),I_top(j,1),C_top(mod(j,Nt)+1,1));
%     fprintf(fileID,'%s %d %d %d\n','f',I_top(j,1),I_top(mod(j,Nt)+1,1),C_top(mod(j,Nt)+1,1));
    jt=jt+1;
    Triset(jt,:) = [C_top(j,1),I_top(j,1),C_top(mod(j,Nt)+1,1)];
    jt=jt+1;
    Triset(jt,:) = [I_top(j,1),I_top(mod(j,Nt)+1,1),C_top(mod(j,Nt)+1,1)];
end

%% bottom cylinder
Nt = size(C_bot,1);
for j=1:Nt
%     fprintf(fileID,'%s %d %d %d\n','f',I_down(j,1),C_bot(j,1),I_down(mod(j,Nt)+1,1));
%     fprintf(fileID,'%s %d %d %d\n','f',C_bot(j,1),C_bot(mod(j,Nt)+1,1),I_down(mod(j,Nt)+1,1));
    jt=jt+1;
    Triset(jt,:) = [I_down(j,1),C_bot(j,1),I_down(mod(j,Nt)+1,1)];
    jt=jt+1;
    Triset(jt,:) = [C_bot(j,1),C_bot(mod(j,Nt)+1,1),I_down(mod(j,Nt)+1,1)];
end

%% side cylinder
Nt = size(C_side,1);
for j=1:Nt
%     fprintf(fileID,'%s %d %d %d\n','f',I_side(j,1),C_side(j,1),I_side(mod(j,Nt)+1,1));
%     fprintf(fileID,'%s %d %d %d\n','f',C_side(j,1),C_side(mod(j,Nt)+1,1),I_side(mod(j,Nt)+1,1));
    jt=jt+1;
    Triset(jt,:) = [I_side(j,1),C_side(j,1),I_side(mod(j,Nt)+1,1)];
    jt=jt+1;
    Triset(jt,:) = [C_side(j,1),C_side(mod(j,Nt)+1,1),I_side(mod(j,Nt)+1,1)];
end

%%
% %% add to top
  addh = 1;
%%
  addtop1 = C_top;
  addtop2 = zeros(Nit,4);
  for k1=1:addN
      itm = addN+1-k1;
      swc1D(itm,:) =[itm,2,0,0,C_top(i,4)+k1*addh,R,itm-1,k1+1];
  for i=1:Nit
     ind = ind+1;
     addtop2(i,:) = [ind,C_top(i,2:3),C_top(i,4)+k1*addh]';
     Vl(ind,:)  = addtop2(i,:); 
     marker3D(ind,:) = [ind,1+k1];
  end
  Nt = Nit;
  for j=1:Nt
    %fprintf(fileID,'%s %d %d %d\n','f',addtop2(j,1),addtop1(j,1),addtop2(mod(j,Nt)+1,1));
    %fprintf(fileID,'%s %d %d %d\n','f',addtop1(j,1),addtop1(mod(j,Nt)+1,1),addtop2(mod(j,Nt)+1,1));
    jt=jt+1;
    Triset(jt,:) = [addtop2(j,1),addtop1(j,1),addtop2(mod(j,Nt)+1,1)];
    jt=jt+1;
    Triset(jt,:) = [addtop1(j,1),addtop1(mod(j,Nt)+1,1),addtop2(mod(j,Nt)+1,1)];
  end
  addtop1 = addtop2;
  end
%% add to bottom   
  addtop1 = C_bot;
  addtop2 = zeros(Nit,4);
  for k1=1:addN
      itm = addN+3+k1;
      swc1D(itm,:) =[itm,2,0,0,C_bot(i,4)-k1*addh,R,itm-1,-1-k1];
  for i=1:Nit
     ind = ind+1;
     addtop2(i,:) = [ind,C_bot(i,2:3),C_bot(i,4)-k1*addh]';
     Vl(ind,:)  = addtop2(i,:); 
     marker3D(ind,:) = [ind,-1-k1];
  end
  Nt = Nit;
  for j=1:Nt
    %fprintf(fileID,'%s %d %d %d\n','f',addtop2(j,1),addtop1(j,1),addtop2(mod(j,Nt)+1,1));
    %fprintf(fileID,'%s %d %d %d\n','f',addtop1(j,1),addtop1(mod(j,Nt)+1,1),addtop2(mod(j,Nt)+1,1));
    jt=jt+1;
    Triset(jt,:) = [addtop1(j,1),addtop2(j,1),addtop1(mod(j,Nt)+1,1)];
    jt=jt+1;
    Triset(jt,:) = [addtop2(j,1),addtop2(mod(j,Nt)+1,1),addtop1(mod(j,Nt)+1,1)];
  end
  addtop1 = addtop2;
  end
  
%% add to side   
  Nt = size(C_side,1);
  addtop1 = C_side;
  addtop2 = zeros(Nt,4);
  for k1=1:addN
      v = v-addh;
      itm = 2*addN+4+k1;
      swc1D(itm,:) =[itm,2,-v*sin(k),0,v*cos(k),r,itm-1,-10000-k1*10000];
      %%
      Ys = [(-v*sin(k)+r*cos(k)*cos(u1))',(r*sin(u1))',(v*cos(k)+r*cos(u1)*sin(k))'];  %side circle up
      Y1s = [(-v*sin(k)+r*cos(k)*cos(u2))',(r*sin(u2))',(v*cos(k)+r*cos(u2)*sin(k))']; %side circle down
      C=[Ys(1:end-1,:);Y1s(end:-1:2,:)];
  for i=1:Nt
     ind = ind+1;
     addtop2(i,:) = [ind,C(i,:)]';
     Vl(ind,:)  = addtop2(i,:); 
     marker3D(ind,:) = [ind,-10000-k1*10000];
  end
  for j=1:Nt
    %fprintf(fileID,'%s %d %d %d\n','f',addtop2(j,1),addtop1(j,1),addtop2(mod(j,Nt)+1,1));
    %fprintf(fileID,'%s %d %d %d\n','f',addtop1(j,1),addtop1(mod(j,Nt)+1,1),addtop2(mod(j,Nt)+1,1));
    jt=jt+1;
    Triset(jt,:) = [addtop1(j,1),addtop2(j,1),addtop1(mod(j,Nt)+1,1)];
    jt=jt+1;
    Triset(jt,:) = [addtop2(j,1),addtop2(mod(j,Nt)+1,1),addtop1(mod(j,Nt)+1,1)];
  end
  addtop1 = addtop2;
  end
  
  
%% Print *.obj file

%cd '/home/qg11b/Desktop/Weak_galerkin_elliptic_obstacle/Mesh2d_v24/Mesh2d v24/test';
name1 = 'Y_shape.obj';
fileID = fopen(name1,'w');  
for i =1:size(Vl,1)
    fprintf(fileID,'%s %12.16f %12.16f %12.16f\n','v',Vl(i,2),Vl(i,3),Vl(i,4)); % vertices
end

for i =1:size(Triset,1)
    fprintf(fileID,'%s %12.16f %12.16f %12.16f\n','f',Triset(i,1),Triset(i,2),Triset(i,3)); % faces
end

fclose(fileID);

%% Print *.swc file and *.txt

%cd '/home/qg11b/Desktop/Weak_galerkin_elliptic_obstacle/Mesh2d_v24/Mesh2d v24/test';
swc1D(1,7) = -1;
name2 = 'Y_shape.swc';
name2txt = 'Y_shape.txt';
fileID = fopen(name2,'w');  
fileID2 = fopen(name2txt,'w');
for i =1:size(swc1D,1)
    fprintf(fileID,'%d %d %12.4f %12.4f %12.4f %12.4f %d %d',...
        swc1D(i,1:7));
    fprintf(fileID,'\n');
        %swc1D(i,1),swc1D(i,2),swc1D(i,3),swc1D(i,4),swc1D(i,5),swc1D(i,6),swc1D(i,7)); % vertices
    fprintf(fileID2,'%d %d %12.4f %12.4f %12.4f %12.4f %d %d',...
        swc1D(i,1:7));
    fprintf(fileID2,'\n');
        
end
fclose(fileID);
fclose(fileID2);

%% Print 1D3D map file

%cd '/home/qg11b/Desktop/Weak_galerkin_elliptic_obstacle/Mesh2d_v24/Mesh2d v24/test';
swc1D(1,7) = -1;
name2 = 'Y_shape_1D3D.map';
name2txt = 'Y_shape_1D3D.txt';
fileID = fopen(name2,'w');
fileID2 = fopen(name2txt,'w');

fprintf(fileID,'%s \n','#index of vertices in swc file, marker'); % vertices, marker
for i =1:size(swc1D,1)
    fprintf(fileID,'%d %d',swc1D(i,1),swc1D(i,8)); % vertices, marker
    fprintf(fileID,'\n');
    
    % for txt version of file
    fprintf(fileID2,'%d %d',swc1D(i,1),swc1D(i,8)); % vertices, marker
    fprintf(fileID2,'\n');
end
fclose(fileID2);

fprintf(fileID,'%s \n','#index of vertices in obj file, marker'); % vertices, marker
for i =1:size(marker3D,1)
    fprintf(fileID,'%d %d',marker3D(i,:)); % vertices, marker
    fprintf(fileID,'\n');
end
fclose(fileID);
%% 
% 
% % comm2 = ['~/.local/bin/tetgen -pq',  ' cylinB1.poly' ];
% % system(comm2);

