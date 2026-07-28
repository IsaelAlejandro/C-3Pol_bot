clear all; close all; clc; format long g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates a library of DSF bases stored as a 450x10 matrix,
% as described in the paper "3D Stokes Polarimetric Imaging at the Nanoscales."
%
% Each element of the library depends on the spatial coordinates of the
% dipole relative to the coverslip.
%
% Each circular polarization channel is represented by a 15x15 pixel array,
% where the center pixel corresponds to the coordinate (x, y) = (0, 0)
% in pixel units.
%
% The axial position is controlled by the parameters dz and df, in microns
% while the lateral position is controlled by rx (in pixel units).
%
% The case dz=df=0,rx = 0 corresponds to a single DSF basis centered at (x, y) = (0, 0).
%
% The case dz=df=0, rx = [-1, 0, 1] generates a set of 9 DSF bases centered at the
% following pixel coordinates:
% (-1,1),  (0,1),  (1,1),
% (-1,0),  (0,0),  (1,0),
% (-1,-1), (0,-1), (1,-1).
% The nine lateral positions are assigned indices from 1 to 9 according
% to the following convention:
%
%   1  2  3    <->   (-1, 1)  (0, 1)  (1, 1)
%   4  5  6    <->   (-1, 0)  (0, 0)  (1, 0)
%   7  8  9    <->   (-1,-1)  (0,-1)  (1,-1)

% We can also generate a set of bases across a volume. For example, by
% setting df = [-0.050, 0.0, 0.050] (in microns) and rx = [-1 0 1] (in pixel units), the script
% generates a set of 27 basis elements. The index assigned to each element
% follows the same ordering described previously within each plane of
% constant df. The numbering then continues sequentially from the plane
% with the lowest df value to the plane with the highest df value.

% The script saves the following files:
%
% I0Rm.mat
%   A (15 × 15 × 10 × length(df) × Total) array, where
%       Total = length(df) * length(rx)^2,
%   representing the total number of basis elements.
%
%   The array contains the 9 right-handed circular polarization basis
%   components, each stored as a 15 × 15 image. The 10th component
%   corresponds to the background.
%
% I0Lm.mat
%   A (15 × 15 × 10 × length(df) ×Total) array, where
%       Total = length(df) * length(rx)^2,
%   representing the total number of basis elements.
%
%   The array contains the 9 left-handed circular polarization basis
%   components, each stored as a 15 × 15 image. The 10th component
%   corresponds to the background.
%
% MBasea.mat
%   Contains the combined basis arrays:
%       Base = [I0Rm I0Lm].
%
% PinvBaseb.mat
%   A 10 × 450 × Total array containing the pseudoinverse matrices
%   associated with each basis, including the constant background
%   component.
%
% MBaseb.mat
%   A 10 × 450 × Total array containing the basis matrices associated
%   with each basis, including the constant background component.
%
% PinvBase.mat
%   A 9 × 450 × Total array containing the pseudoinverse matrices
%   associated with each basis, excluding the constant background
%   component.
%
% MBase.mat
%   A 9 × 450 × Total array containing the basis matrices associated
%   with each basis, excluding the constant background component.

L=2^9;   %Grid length in pixels
pixel_size=0.0527242; %camera pixel size in microns

NA =1.49 ; %Objective numerical aperture
l=0.5320;  %wavelength  in microns    
n1=1.518;  %index of refraction of dipole medium
n2=1.518; % index of refraction of glass

reduction_factor=1;%%reduction factor of the numercial aperture at the pupil plane, if there is not reduction, this parameter is 1
dz=0; %distance between dipole and coverslip in microns.
df=[-0.05,0.0,0.05]; %distance between coverslip and focal plane in microns.  it can be a vector
%% Spatial coordinates%%%%%%%%%%%%%%
rx=[-1,0,1];  %it can be a vector, In pixel units.




totalind=0;
for plane=1:length(df)
ind=1;
%%%%%%%%%%%%%%%%% Initial Parameters %%%%%%%%%%%%%%%%%
factorr=4;
d=40*factorr; %actual dimension 2*d+1;
totalN=2*d+1;
dfocus=df(plane); %dfocus<0, coverslip closer to MO, dfocus>0 coverllip above NFP.
dx=15*factorr;
dxf=pixel_size/factorr; 
rx=dxf*rx*factorr;
ry=-rx;


for iy=1:length(rx)
    for ix=1:length(rx)
      tic
r0=[rx(ix),ry(iy),dz;rx(ix),ry(iy),dz;rx(ix),ry(iy),dz];

p0=[1,0,0;0,1,0;0,0,1];
ap_xy=1;
ap_ad=1;
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=1000;   %focal lenght in microns
[Gt0,M,k1,k2,kz1,kz2,x,y,phi_u,th_u,r_u,acf2,rho,a,a_xy,a_ad,Rmax]=Green_Parameters(L,dxf,NA,n1,n2,f,l,dfocus,reduction_factor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:3
[Ekx,Eky]=Field_Pupil(r0(ii,:),p0(ii,:),Gt0,k1,kz1,f,x,y,phi_u,th_u,r_u,acf2);
Ex10=a.*Ekx;  Ey10=a.*Eky;
phi=mod(atan2(y,x),2*pi); u=sqrt((x./Rmax).^2+(y./Rmax).^2);
%%%%%%%%%%%%%%%%%%%%%
ax=-1i*k1.*((x/f)); ay=-1i*k1.*((y/f)); az=1i*kz1;
% %%%%%%%%%%%%%%%%%%%%% PSF computation %%%%%%%%%%%%%%%
Eimx(:,:,ii)=(fftshift(ifft2(ifftshift(Ex10))));
Eimy(:,:,ii)=(fftshift(ifft2(ifftshift(Ey10))));
%%%%%%%%%%%%%%%%%%%%%%%%PSF derivatives computation %%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[I0RRm,I0LLm,EE]=base_stokes_derivatives_better2(Eimx,Eimy,L,d,ind,factorr);
I0Rm00(:,:,:,plane)=I0RRm;
I0Lm00(:,:,:,plane)=I0LLm;

   
        totalind=totalind+1; 
        newcenterx=(totalN+1)/2; 
        newcentery=(totalN+1)/2; 
        for elementss=1:10
            if elementss==1

        I0Rm(:,:,elementss,plane,totalind)=I0Rm00(:,:,elementss,plane);
        I0Lm(:,:,elementss,plane,totalind)=I0Lm00(:,:,elementss,plane);
        I0_0=[ I0Rm(:,:,elementss,plane,totalind)  I0Lm(:,:,elementss,plane,totalind)];
        I00=sum(I0_0(:));
         I0Rm(:,:,elementss,plane,totalind)=I0Rm(:,:,elementss,plane,totalind)./I00;
        I0Lm(:,:,elementss,plane,totalind)=I0Lm(:,:,elementss,plane,totalind)./I00;
            end
        I0Rm(:,:,elementss,plane,totalind)=I0Rm00(:,:,elementss,plane)./I00;
        I0Lm(:,:,elementss,plane,totalind)=I0Lm00(:,:,elementss,plane)./I00;
      


        end
toc
    end
 end

end



save('I0Rm.mat','I0Rm', '-v7.3');
save('I0Lm.mat','I0Lm', '-v7.3');
base=[I0Rm I0Lm];
save('MBasea.mat','base', '-v7.3');
base=[I0Rm I0Lm];
[n1,n2,n3,n4,n5]=size(base);

k=0;
for j=1:n5
if mod(j,n5/n4)==1  
   k=k+1;
end

if n5==1 && n4==1
   k=k+1;
end

for i=1:9
    M11=base(:,:,i,k,j);
    MK(:,i)=M11(:);

end

MBase(:,:,j)=MK;
M1=MK';
M2=(MK')*MK;
M2inv=inv(M2);
PinvBase(:,:,j)=(M2inv)*M1;

end
save('PinvBase.mat','PinvBase', '-v7.3');
save('MBase.mat','MBase', '-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k=0;
for j=1:n5
if mod(j,n5/n4)==1  
   k=k+1;
end

if n5==1 && n4==1
   k=k+1;
end

for i=1:10
    M11=base(:,:,i,k,j);
    MK(:,i)=M11(:);

end

MBaseb(:,:,j)=MK;
M1=MK';
M2=(MK')*MK;
M2inv=inv(M2);
PinvBaseb(:,:,j)=(M2inv)*M1;

end
save('PinvBaseb.mat','PinvBaseb', '-v7.3');
save('MBaseb.mat','MBaseb', '-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imagesc(MBaseb(:,:,13)); colormap(jet)
