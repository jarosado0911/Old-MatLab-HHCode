function [sysMat,A,B]=get_stencils(c,k,sz,mthd)

sten = [1 -2 1];
sysMat = spdiags(ones(sz,1)*sten,-1:1,sz,sz)*c; %sparse diagonal

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

