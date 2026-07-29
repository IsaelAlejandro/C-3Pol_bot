function[major_axis, minor_axis]=ellipse3D(phi,theta,psi,rba)
phi=phi*pi/180;
theta=theta*pi/180;
psi=psi*pi/180;

ax=cos(phi)*sin(theta);
ay=sin(phi)*sin(theta);
az=cos(theta);
major_axis = [ax, ay, az]; major_axis=major_axis./norm(major_axis);  % Major axis direction and length
if  abs(az)<1
minor_axis0=[-ay, ax,0]; minor_axis0=minor_axis0./norm(minor_axis0);
else
minor_axis0=[0, 1,0]; 
end
K=[0,-az, ay; az, 0, -ax; -ay, ax, 0];
I=[1, 0, 0; 0, 1, 0; 0, 0, 1];
R=I+sin(psi)*K+(1-cos(psi))*K*K;

minor_axis =(R*minor_axis0')'; minor_axis = rba*minor_axis./norm(minor_axis);% Minor axis direction and length %
normab=sqrt(norm(minor_axis)^2+norm(major_axis)^2);
minor_axis=minor_axis/normab;
major_axis=major_axis/normab;

end