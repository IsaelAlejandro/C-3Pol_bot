function[Gt0,M,k1,k2,kz1,kz2,x,y,phi_u,th_u,r_u,acf2,rho,a,a_xy,a_ad,Rmax]=Green_Parameters(L,dxf,NA,n1,n2,f,l,dfocus,reduction_factor)
M=round(2*NA*L*dxf/(l));    %Pupil diameter in pixels 240 is the magnification  %IT SHOULD BE MULTIPLIDED BY N2

if mod(M,2)==1
    M=M-1;
end
%%%%%%%%%%%%%%%%%%%%%% Spatial & angular coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rmax=f*NA/n2; dx=Rmax/(M/2); % pixel size
X=(1:L)-L/2-.5; X=dx*X; [x,y]=meshgrid(X,-X); %Integration grid in spatial coordinates 
z=-sqrt(f.^2-x.^2-y.^2); %z coordinate
a_xy=zeros(L,L);
a_ad=zeros(L,L);
a=zeros(L,L);
for ii=1:L
    for jj=1:L
%    if sqrt((x(ii,jj)/dx).^2+(y(ii,jj)/dx).^2)<=M/2
   if sqrt((x(ii,jj)/dx).^2+(y(ii,jj)/dx).^2)<reduction_factor*M/2 %scatter .82
   a(ii,jj)=1;
    a_xy(ii,jj)=1;
      a_ad(ii,jj)=1;
   end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% auxiliar coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=sqrt(x.^2+y.^2); 
sph=y./rho;
cph=x./rho;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% angular coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k0=2*pi/l; k1=n1*k0; k2=n2*k0; %wave numbers in vaccun, top medium, bottom medium
k_t=k2.*rho./f;
kz1=sqrt(k1.^2-k_t.^2);
kz2=sqrt(k2.^2-k_t.^2);
%%%%%%%%%%%%% Fresnel coefficients (interface) %%%%%%%%%%%%%%%%
m1=1; e1=n1^2; m2=1; e2=n2^2;  %permeability and permitivity in both mediums 
ts=2*m2*kz1./(m2.*kz1+m1.*kz2);
tp=(sqrt((m2.*e1)/(m1.*e2))).*(2.*e2.*kz1./(e2.*kz1+e1.*kz2));
%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
rho2=rho.^2;
Phi_1=(tp.*(n2/n1).*(k2.*z./f))./(kz1);
Phi_2=-tp.*(n2/n1);
Phi_3=(ts.*(k2.*z./f))./(kz1);
cth=z./sqrt(rho2+z.^2);
sth=sqrt(1.^2-cth.^2);
cf2=-1i.*sqrt(n2)*(cth).^(-.5); %conservation energy
xy=x.*y;
xz=x.*z;
zy=y.*z;
xx=x.^2;
yy=y.^2;
zz=z.^2;
P2_r2_f2=(Phi_2)./((rho2).*(f^2));
P3_r2=(Phi_3)./((rho2));
%%%%%%%%%%%%%%%%%%%%     Green tensor components        %%%%%%%%%%%%%%%%%%%%
Gt0(:,:,1)=((xx).*(zz)).*P2_r2_f2+((yy).*P3_r2);
Gt0(:,:,2)=((xy).*(zz)).*P2_r2_f2-((xy).*P3_r2);
Gt0(:,:,3)=-((xz).*(Phi_2))./(f.^2);
Gt0(:,:,4)=((xy).*(zz)).*P2_r2_f2-((xy).*P3_r2);
Gt0(:,:,5)=((yy).*(zz)).*P2_r2_f2+((xx).*P3_r2);
Gt0(:,:,6)=-((zy).*(Phi_2))./(f.^2);
Gt0(:,:,7)=-((xz).*(Phi_1))./(f.^2);
Gt0(:,:,8)=-((zy).*(Phi_1))./(f.^2);
Gt0(:,:,9)=(1-((zz)./(f.^2))).*(Phi_1);
m0=.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gt0=(exp(1i*kz2.*dfocus)./(4*pi*f)).*(exp(1i*k2*f)).*Gt0.*exp(-1i*(m0*2*pi*x/(dx*L)-m0*2*pi*y/(dx*L)));
% Gt0=(exp(1i*kz2.*dfocus)./(4*pi*f)).*(exp(1i*k2*f)).*Gt0;

%%%%%%%%%%%  Unitary vectors phi, rho (cylindrical basis) and theta (spherical basis) %%%%%%%%%%%%%%%%%%%
phi_u(:,:,1)=-sph;
phi_u(:,:,2)=cph;
phi_u(:,:,3)=zeros(L,L);
th_u(:,:,1)=cph.*cth;
th_u(:,:,2)=sph.*cth;
th_u(:,:,3)=-sth;
r_u(:,:,1)=cph;
r_u(:,:,2)=sph;
r_u(:,:,3)=zeros(L,L);
acf2=a.*cf2;

end