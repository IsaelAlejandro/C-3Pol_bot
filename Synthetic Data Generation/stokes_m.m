function[s]=stokes_m(Mdep)
[n1,n2,n3]=size(Mdep);
for i=1:n3
mm2=Mdep(:,:,i);
mm=mm2/trace(mm2);
s(i,1)=trace(mm);
s(i,2)=mm(1,1)-mm(2,2);
s(i,3)=(mm(1,1)+mm(2,2)-2*mm(3,3))/sqrt(3);
s(i,4)=2*real(mm(3,2));
s(i,5)=2*real(mm(1,3));
s(i,6)=2*real(mm(2,1));
s(i,7)=2*imag(mm(3,2));
s(i,8)=2*imag(mm(1,3));
s(i,9)=2*imag(mm(2,1));
s(i,:)=(sqrt(3)/2)*s(i,:)./s(i,1);
s(i,1)=1;
s(i,:)=real(s(i,:));
end

end