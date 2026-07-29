clear all; close all; clc; format long


load I0Lm.mat 
load I0Rm.mat
load MBaseb.mat
load PinvBaseb.mat


for indx=1:1
    
    name=['I_Cf' num2str(indx) '.mat'];
    load(name)


    base=[I0Rm I0Lm];
    [n1,n2,n3,n4,n5]=size(base);
    [nr,nc]=size(base(:,:,1,1,1)); denoIk=zeros(nr,nc);
    CF=[];
    FF=[];
    FF1=[];
    n1n=0;
    indd=length(I_Cf);
    totalframe=0;
    Imff=I_Cf;    

parfor i=1:indd %parallel loop

           [xi1,eta1,mm1,s11,indmax,III,plane,am,bm,DegP,FFF,NTphotons]  =  inv_stokes(Imff(:,:,i),base,denoIk,PinvBaseb,MBaseb,n1,n2,n4,n5);
           backPh{i}=FFF;
           NPhot{i}=NTphotons;
           c_Imfs{i}=III;
           c_aam{i}=am;
           c_bbm{i}=bm;
           c_xi{i}=xi1;
           c_eta{i}=eta1;
           c_mm{i}=mm1;
           c_zcoor{i}=plane;
           c_posr2{i}=indmax;
           c_s{i}=s11;
           c_degpol{i}=DegP;


end

for ii = 1:indd
       i=ii;
        NPhoton(i)=NPhot{ii};
        BNPho(i)=backPh{ii};
        Imfs(:,:,i)=c_Imfs{ii};
        aam(:,i)=c_aam{ii};
        bbm(:,i)=c_bbm{ii};
        xi(:,i)=c_xi{ii};
        eta(:,i)=c_eta{ii};
        mm(:,:,i)= c_mm{ii};
        zcoor(:,i)=c_zcoor{ii};
        posr2(:,i)=c_posr2{ii};
        s(i,:)=c_s{ii};
        degpol(:,i)=c_degpol{ii};


    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1_NPhoton=NPhoton; name=['s1_NPhoton' num2str(indx) '.mat']; save(name,'s1_NPhoton'); %Signal Photons
s1_BNPho=BNPho; name=['s1_BNPho' num2str(indx) '.mat']; save(name,'s1_BNPho'); %Background Photons /pixels
s1_Imfs=Imfs; name=['Imfs' num2str(indx) '.mat']; save(name,'Imfs'); % Best fit
s1_aam=aam; name=['s1_aam' num2str(indx) '.mat']; save(name,'s1_aam'); %major axis
s1_bbm=bbm; name=['s1_bbm' num2str(indx) '.mat']; save(name,'s1_bbm');%minor axis
s1_xi=xi; name=['s1_xi' num2str(indx) '.mat']; save(name,'s1_xi'); %in plane angle \xi
s1_eta=eta; name=['s1_eta' num2str(indx) '.mat']; save(name,'s1_eta'); %off plane angle \eta
s1_mm=mm;   name=['s1_mm' num2str(indx) '.mat']; save(name,'s1_mm'); %second moment matrix
s1_zcoor=zcoor; name=['s1_zcoor' num2str(indx) '.mat']; save(name,'s1_zcoor'); %index associated with the axial coordinate of the bases 
s1_s=s; name=['s1_s' num2str(indx) '.mat']; save(name,'s1_s'); %normalized estimated stokes gell mann parameters
s1_degpol=degpol; name=['s1_degpol' num2str(indx) '.mat']; save(name,'s1_degpol'); %degree of polarization
s1_posr2=posr2; name=['s1_posr2' num2str(indx) '.mat']; save(name,'s1_posr2');% global index of the basis that better fits the data. 
end
% 
