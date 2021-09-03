function Frustum=Frustum(X1,X2,r1,r2,n,cyl_color,lines)
frustumLength = norm(X1-X2);
t = linspace(0,frustumLength,n);

mslope= (r2-r1)/(t(end)-t(1));
theLine = mslope*(t-t(1))+r1;
[X,Y,Z]=cylinder(theLine,20);
Z=frustumLength*Z;
Frustum=mesh(X,Y,Z);
xlabel('x')
ylabel('y')
zlabel('z')

% Defining Unit vector along the X-direction
unit_Vx=[0 0 1];
% Calulating the angle between the x direction and the required direction
% of cylinder through dot product
angle_X1X2=acos( dot( unit_Vx,(X2-X1) )/( norm(unit_Vx)*norm(X2-X1)) )*180/pi;
% Finding the axis of rotation (single rotation) to roate the cylinder in
% X-direction to the required arbitrary direction through cross product
axis_rot=cross(unit_Vx,(X2-X1) );
% Rotating the plotted cylinder and the end plate circles to the required
% angles
if angle_X1X2~=0 % Rotation is not needed if required direction is along X
    rotate(Frustum,axis_rot,angle_X1X2,[0 0 0])
    
end
set(Frustum,'XData',get(Frustum,'XData')+X1(1))
set(Frustum,'YData',get(Frustum,'YData')+X1(2))
set(Frustum,'ZData',get(Frustum,'ZData')+X1(3))
% Setting the color to the cylinder and the end plates
set(Frustum,'FaceColor',cyl_color)
% If lines are not needed making it disapear
if lines==0
    set(Frustum,'EdgeAlpha',0)
end

end

