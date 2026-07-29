function [s,mm] = physical_Stokes(bs00);

ss(1:9)=bs00(1:9);
s(1:9)=ss(1:9);
%Normalizada%
%no normalizada
mxx=(2/3)*s(1)+s(2)+s(3)/sqrt(3);
myy=(2/3)*s(1)-s(2)+s(3)/sqrt(3);
mzz=(2/3)*s(1)-2*s(3)/sqrt(3);
mxy=s(6)-1i*s(9); myx=conj(mxy);
mxz=s(5)+1i*s(8); mzx=conj(mxz);
myz=s(4)-1i*s(7); mzy=conj(myz);
mm=0.5*[mxx mxy mxz; myx myy myz; mzx mzy mzz];
% [V,D] =eig(mm);
[V,D0] =eig(mm);

% mm2=zeros(3,3);
D0(D0<0)=0;

mm=V*D0*V';
% mm=mm2;
s(1)=trace(mm);
s(2)=mm(1,1)-mm(2,2);
s(3)=(mm(1,1)+mm(2,2)-2*mm(3,3))/sqrt(3);
s(4)=2*real(mm(2,3));
s(5)=2*real(mm(3,1));
s(6)=2*real(mm(1,2));
s(7)=2*imag(mm(2,3));
s(8)=2*imag(mm(3,1));
s(9)=2*imag(mm(1,2));

% 
% s=(sqrt(3)/2)*s./s(1);
% s(1)=1;
s=real(s);
end