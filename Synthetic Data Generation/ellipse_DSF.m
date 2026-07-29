clc; clear; close all; format long

% This code generates N_rep realizations of the Dipole Spread Function (DSF)
% for a given rigid dipole configuration defined by the parameters ξ, η, ψ, and r_ba.
%
% The simulation accounts for a specified number of detected signal photons (NP)
% and a mean background photon count per pixel (bn), both modeled using
% Poisson statistics.
%
% The parametrization of the dipole elliptical state in terms of angles is described in the
% Supplementary Material of the paper:
% "3D Stokes Polarimetric Imaging at the Nanoscales."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Mbasea.mat %load the SGM basis

N_rep=100; %realizations
NP=7500; %number of signal photons
bn=10;  %mean background photon count per pixel
%Ellipse parameters: 
r_ba=0; %b/a, ratio between the minor/major axis length, the sign indicates the handennes.
xi=45;  %in degrees
eta=45; %in degrees
psi=0;  %in degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Degree of polarization: 
%L2=L3=0, corresponds to the case of a fully polarized state. 
L2=0;
L3=L2;
L1=1-L2-L3; 


for j0=1:length(NP)
    k=0;
    ll=0;
   
    for j1=1:length(r_ba)
        for j2=1:length(eta)
            for j3=1:length(psi)           
                for j4=1:length(xi)
                    for j6=1:length(L2)
                        ll=1+ll;
                        
                        
                        [major_axis, minor_axis]=ellipse3D(xi(j4),eta(j2),psi(j3),r_ba(j1)); % Minor and major axis calculation
                        
                        [Mdep]=second_moment_matrix_iso(major_axis,minor_axis,L1(j6),L2(j6),L3(j6));
                        m2m(:,:,ll)=Mdep; %second moment matrix
                        s(ll,:)=mm2_to_s(m2m(:,:,ll)); %SGM parameters from the second moment matrix
                        %Ideal DSF calculation
                        DSF0=base(:,:,1); 
                        for ii=1:8
                        DSF0=DSF0+(2/sqrt(3)).*s(ll,ii+1).*base(:,:,ii+1);
                        end
                        DSF0=DSF0*NP(j0)/sum(DSF0(:))+bn;
                        for j5=1:N_rep
                           
                            k=k+1;
                      
                            I_Cf(:,:,k)=poissrnd(DSF0);                     
                
                       end
                   end         
               end
           end
        end
   end
  
    name2=['I_Cf' num2str(j0) '.mat']; %Dipole circular DSFs
    save(name2,'I_Cf')
    

end

plot_ellipse3D(major_axis,minor_axis); % Ellipse visualization.
title('Elliptical Dipole');

figure
imagesc(DSF0)
title('Ideal Circ DSF');

figure
imagesc(I_Cf(:,:,1))
title('Synthetic circular Dipole Spread Function');



