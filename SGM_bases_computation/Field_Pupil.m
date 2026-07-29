function[Ex,Ey]=Field_Pupil(r0,p0,Gt0,k1,kz1,f,x,y,phi_u,th_u,r_u,acf2)
[n,m]=size(r0); N_dipoles=n;
L=length(x);
Ex=zeros(L,L); Ey=zeros(L,L);
for ii=1:N_dipoles    
px=p0(ii,1);    
py=p0(ii,2);
pz=p0(ii,3);
x0=r0(ii,1);    
y0=r0(ii,2);
z0=r0(ii,3);
Gt=(exp(-1i.*(k1.*((x/f).*(x0)+(y/f).*(y0))-z0.*kz1))).*Gt0;
%%%%%%%%%  Far field of each dipole %%%%%%%%%%%
E_far(:,:,1)=(px*Gt(:,:,1)+py*Gt(:,:,4)+pz*Gt(:,:,7));
E_far(:,:,2)=(px*Gt(:,:,2)+py*Gt(:,:,5)+pz*Gt(:,:,8));
E_far(:,:,3)=(px*Gt(:,:,3)+py*Gt(:,:,6)+pz*Gt(:,:,9));
%%%%%%%%% Projection of the far field into phi and theta directions %%%%%%%
Eip0=E_far.*phi_u;
Eit0=E_far.*th_u;
Eip=zeros(L,L);
Eit=zeros(L,L);
for jj=1:3
Eip=Eip+Eip0(:,:,jj);
Eit=Eit+Eit0(:,:,jj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Transformation of the theta compo into -rho vector and phi into phi.
Epf=Eip.*phi_u;
Etf=Eit.*(-r_u); %we are propagating the field in the -z direction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Multiplication by the aperture and the term cf2 for energy conservation
Exp0=acf2.*(Epf(:,:,1)+Etf(:,:,1));
Eyp0=acf2.*(Epf(:,:,2)+Etf(:,:,2));
%%%%%%%%%%%%%%Field at pupil and plane %%%%%%%%%%%%%%%

Ex=Ex+Exp0;
Ey=Ey+Eyp0;
end
