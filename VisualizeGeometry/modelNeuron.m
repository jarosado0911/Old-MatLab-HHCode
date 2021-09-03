function combinedobject = modelNeuron()

if nargin == 0
    filename = 'grids/geometry6.swc';
end

[~,id,pid,coord,r,~]=readSWC(filename);
[x,y,z]=sphere;

figure(1)
set(gcf, 'Position',  [100, 100, 2000, 900]);
myax = axes();
view(3);
axis equal;
hold on 

mysRad = r(1);
h(1)=surf(x*mysRad+coord(id(1),1),y*mysRad+coord(id(1),2),z*mysRad+coord(id(1),3));
for i=2:length(id)
    fprintf(sprintf('piece %i\n',i))
    myRad=r(i);
    surf(x*myRad+coord(id(i),1),y*myRad+coord(id(i),2),z*myRad+coord(id(i),3),'edgecolor','none');
    h(i)=Frustum(coord(id(i),:),coord(pid(i),:),r(id(i)),r(pid(i)),10,'b',0);
end
combinedobject = hgtransform('parent',myax);
set(h,'parent',combinedobject);
set(gcf,'color',[1 1 1]*0.5);
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
axis equal
axis vis3d
axis off
light
shading interp
lighting gouraud
end