function [xi,eta,mm,s,indmax,Ibest,planeee,am,bm,DegP,FFF0,NTphotons] = inv_stokes(file,base,denoIk,PinvBaseb,MBaseb,n1,n2,n4,n5)
Imff=file;


k=0;
for j=1:n5
    if mod(j,n5/n4)==1  
       k=k+1;
    end

     if n5==1 && n4==1  
       k=k+1;
    end
    bs00(:,j)=PinvBaseb(:,:,j)*Imff(:);

    [bs00(1:9,j),mm0(:,:,j)]= physical_Stokes(bs00(1:9,j));

    FFF(j)=bs00(10,j)*MBaseb(1,10,j);
 
    Imff0=Imff-FFF(j);
    NPh(j)=sum(sum(Imff0));
    B(j)=(sum(sum(Imff))-NPh(j))/(n1*n2);
    bs0(1:9,j)=bs00(1:9,j);
    Imff2=Imff0;
    [Imfss(:,:,j)]=simPSF_arbitrary_base(bs0(:,j)',base(:,:,:,k,j),denoIk);
    Imfsv= Imfss(:,:,j).^1;
    ccc(j)=(1)*((Imff2(:)'*Imfsv(:))/(sqrt(Imff2(:)'*Imff2(:))*sqrt(Imfsv(:)'*Imfsv(:))));
end
[valmax,indmax]=max(ccc);
% Ibest=Imfss(:,:,indmax);
FFF0=B(indmax);
NTphotons=NPh(indmax);
for kk=1:n4
    if (kk-1)*(n5/n4)+1<=indmax && indmax<=(n5/n4)*kk
       planeee=kk;
    end
end


s(1:9)=bs0(1:9,indmax);
s01=bs0(1:9,indmax);
mm=mm0(:,:,indmax);
Ibest=simPSF_arbitrary_base(s,base(:,:,:,planeee,indmax),denoIk);
Ibest=(Ibest);

s=(sqrt(3)/2)*s./s(1);
s(1)=1;

s=real(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,D] =eig(mm);
[i1,i2]=max([D(1,1) D(2,2) D(3,3)]);
normL=D(1,1)+D(2,2)+D(3,3);
l1=D(1,1)/normL;
l2=D(2,2)/normL;
l3=D(3,3)/normL;
DegP=norm(s(2:9));
A=sqrt(sum(conj(V(:,i2)).*V(:,i2)));
EE=sum(V(:,i2).*V(:,i2));
PHII=mod(0.5*atan2(imag(EE),real(EE)),2*pi);
f=V(:,i2)*exp(-1i*PHII)/(A);
am=real(f);
bm=imag(f);
if am(3,1)<0
    am=-am;
end
eta=acos(abs(am(3,1))/sqrt(am(1,1)^2+am(2,1)^2+am(3,1)^2))*180/pi;
xi=mod(atan2(am(2,1),am(1,1)),2*pi)*180/pi;
if xi<0
    xi=xi+360;
end
xi=mod(xi,360);




end