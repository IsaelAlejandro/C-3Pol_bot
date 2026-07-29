function [I0Rm,I0Lm,EfR]=base_stokes_derivatives_better2(Eimx,Eimy,L,d,ind,factorr)
d=(factorr)*((15-1)/2)+factorr/2;
if ind==1

EfR=(1/sqrt(2)).*(Eimx(L/2+2-d:L/2+1+d,L/2+2-d:L/2+1+d,:)-1i.*Eimy(L/2+2-d:L/2+1+d,L/2+2-d:L/2+1+d,:));
EfL=(1/sqrt(2)).*(Eimx(L/2+2-d:L/2+1+d,L/2+2-d:L/2+1+d,:)+1i.*Eimy(L/2+2-d:L/2+1+d,L/2+2-d:L/2+1+d,:));



else 
EfR=Eimx(L/2+1-d:L/2+1+d,L/2+1-d:L/2+1+d,:);
EfL=Eimy(L/2+1-d:L/2+1+d,L/2+1-d:L/2+1+d,:);



end

I0Rx=(EfR(:,:,1).*conj(EfR(:,:,1)));
I0Ry=(EfR(:,:,2).*conj(EfR(:,:,2)));
I0Rz=(EfR(:,:,3).*conj(EfR(:,:,3)));
I0Rpyz=(EfR(:,:,2)+EfR(:,:,3)).*conj(EfR(:,:,2)+EfR(:,:,3))/2;
I0Rpzx=(EfR(:,:,3)+EfR(:,:,1)).*conj(EfR(:,:,3)+EfR(:,:,1))/2;
I0Rpxy=(EfR(:,:,1)+EfR(:,:,2)).*conj(EfR(:,:,1)+EfR(:,:,2))/2;
I0Rmyz=(EfR(:,:,2)-EfR(:,:,3)).*conj(EfR(:,:,2)-EfR(:,:,3))/2;
I0Rmzx=(EfR(:,:,3)-EfR(:,:,1)).*conj(EfR(:,:,3)-EfR(:,:,1))/2;
I0Rmxy=(EfR(:,:,1)-EfR(:,:,2)).*conj(EfR(:,:,1)-EfR(:,:,2))/2;
I0Rryz=(EfR(:,:,2)-1i*EfR(:,:,3)).*conj(EfR(:,:,2)-1i*EfR(:,:,3))/2;
I0Rrzx=(EfR(:,:,3)-1i*EfR(:,:,1)).*conj(EfR(:,:,3)-1i*EfR(:,:,1))/2;
I0Rrxy=(EfR(:,:,1)-1i*EfR(:,:,2)).*conj(EfR(:,:,1)-1i*EfR(:,:,2))/2;
I0Rlyz=(EfR(:,:,2)+1i*EfR(:,:,3)).*conj(EfR(:,:,2)+1i*EfR(:,:,3))/2;
I0Rlzx=(EfR(:,:,3)+1i*EfR(:,:,1)).*conj(EfR(:,:,3)+1i*EfR(:,:,1))/2;
I0Rlxy=(EfR(:,:,1)+1i*EfR(:,:,2)).*conj(EfR(:,:,1)+1i*EfR(:,:,2))/2;

I0Lx=(EfL(:,:,1).*conj(EfL(:,:,1)));
I0Ly=(EfL(:,:,2).*conj(EfL(:,:,2)));
I0Lz=(EfL(:,:,3).*conj(EfL(:,:,3)));
I0Lpyz=(EfL(:,:,2)+EfL(:,:,3)).*conj(EfL(:,:,2)+EfL(:,:,3))/2;
I0Lpzx=(EfL(:,:,3)+EfL(:,:,1)).*conj(EfL(:,:,3)+EfL(:,:,1))/2;
I0Lpxy=(EfL(:,:,1)+EfL(:,:,2)).*conj(EfL(:,:,1)+EfL(:,:,2))/2;
I0Lmyz=(EfL(:,:,2)-EfL(:,:,3)).*conj(EfL(:,:,2)-EfL(:,:,3))/2;
I0Lmzx=(EfL(:,:,3)-EfL(:,:,1)).*conj(EfL(:,:,3)-EfL(:,:,1))/2;
I0Lmxy=(EfL(:,:,1)-EfL(:,:,2)).*conj(EfL(:,:,1)-EfL(:,:,2))/2;
I0Lryz=(EfL(:,:,2)-1i*EfL(:,:,3)).*conj(EfL(:,:,2)-1i*EfL(:,:,3))/2;
I0Lrzx=(EfL(:,:,3)-1i*EfL(:,:,1)).*conj(EfL(:,:,3)-1i*EfL(:,:,1))/2;
I0Lrxy=(EfL(:,:,1)-1i*EfL(:,:,2)).*conj(EfL(:,:,1)-1i*EfL(:,:,2))/2;
I0Llyz=(EfL(:,:,2)+1i*EfL(:,:,3)).*conj(EfL(:,:,2)+1i*EfL(:,:,3))/2;
I0Llzx=(EfL(:,:,3)+1i*EfL(:,:,1)).*conj(EfL(:,:,3)+1i*EfL(:,:,1))/2;
I0Llxy=(EfL(:,:,1)+1i*EfL(:,:,2)).*conj(EfL(:,:,1)+1i*EfL(:,:,2))/2;




for i = 1:15
    for j = 1:15
        % Extract the 4x4 block
        bR = I0Rx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRx(i, j) = sum(bR(:));
        bL = I0Lx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILx(i, j) = sum(bL(:));
        bR = I0Ry((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRy(i, j) = sum(bR(:));
        bL = I0Ly((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILy(i, j) =  sum(bL(:));
        bR= I0Rz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);   IRz(i, j)= sum(bR(:));
        bL = I0Lz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILz(i, j)= sum(bL(:));
        bR=I0Rpyz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRpyz(i, j)= sum(bR(:));
        bR=I0Rmyz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRmyz(i, j)=sum(bR(:));
        bL=I0Lpyz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILpyz(i, j)= sum(bL(:));
        bL=I0Lmyz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILmyz(i, j)= sum(bL(:));
        bR=I0Rpzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRpzx(i, j)= sum(bR(:));
        bR=I0Rmzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRmzx(i, j)= sum(bR(:));
        bL=I0Lpzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILpzx(i, j)= sum(bL(:));
        bL=I0Lmzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILmzx(i, j)= sum(bL(:));
        bR=I0Rpxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRpxy(i, j)= sum(bR(:));
        bR=I0Rmxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRmxy(i, j)= sum(bR(:));
        bL=I0Lpxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILpxy(i, j)= sum(bL(:));
        bL=I0Lmxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILmxy(i, j)= sum(bL(:));
        bR=I0Rryz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRryz(i, j)= sum(bR(:));
        bR=I0Rrzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRrzx(i, j)= sum(bR(:));
        bR=I0Rrxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRrxy(i, j)= sum(bR(:));
        bR=I0Rlyz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRlyz(i, j)= sum(bR(:));
        bR=I0Rlzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRlzx(i, j)= sum(bR(:));
        bR=I0Rlxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  IRlxy(i, j)= sum(bR(:));
        bL=I0Lryz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILryz(i, j)= sum(bL(:));
        bL=I0Lrzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILrzx(i, j)= sum(bL(:));
        bL=I0Lrxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILrxy(i, j)= sum(bL(:));
        bL=I0Llyz((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILlyz(i, j)= sum(bL(:));
        bL=I0Llzx((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILlzx(i, j)= sum(bL(:));
        bL=I0Llxy((i-1)*factorr+1:i*factorr, (j-1)*factorr+1:j*factorr);  ILlxy(i, j)= sum(bL(:));


    end
end


I0R(:,:,1)=(IRx+IRy+IRz)/3;
I0L(:,:,1)=(ILx+ILy+ILz)/3;
I0R(:,:,2)=(IRx-IRy)/2;
I0L(:,:,2)=(ILx-ILy)/2;
I0R(:,:,3)=(IRx+IRy-2*IRz)./(2*sqrt(3));
I0L(:,:,3)=(ILx+ILy-2*ILz)./(2*sqrt(3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%
I0R(:,:,4)=.5*(IRpyz-IRmyz);
I0L(:,:,4)=.5*(ILpyz-ILmyz);
I0R(:,:,5)=.5*(IRpzx-IRmzx);
I0L(:,:,5)=.5*(ILpzx-ILmzx);
I0R(:,:,6)=.5*(IRpxy-IRmxy);
I0L(:,:,6)=.5*(ILpxy-ILmxy);
I0R(:,:,7)=.5*(IRryz-IRlyz);
I0L(:,:,7)=.5*(ILryz-ILlyz);
I0R(:,:,8)=.5*(IRrzx-IRlzx);
I0L(:,:,8)=.5*(ILrzx-ILlzx);
I0R(:,:,9)=.5*(IRrxy-IRlxy);
I0L(:,:,9)=.5*(ILrxy-ILlxy);
I0R(:,:,10)=max(max(I0R(:,:,1)))*ones(15,15)/8;
I0L(:,:,10)=max(max(I0R(:,:,1)))*ones(15,15)/8;
I0=[I0R(:,:,1) I0L(:,:,1)];
I00=sum(I0(:));
I0Rm=I0R./I00;
I0Lm=I0L./I00;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dz

% 