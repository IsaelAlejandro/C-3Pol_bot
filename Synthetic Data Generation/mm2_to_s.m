function[s]=mm2_to_s(mm2);
mm=mm2/trace(mm2);
s(1)=trace(mm);
s(2)=mm(1,1)-mm(2,2);
s(3)=(mm(1,1)+mm(2,2)-2*mm(3,3))/sqrt(3);
s(4)=2*real(mm(2,3));
s(5)=2*real(mm(3,1));
s(6)=2*real(mm(1,2));
s(7)=2*imag(mm(2,3));
s(8)=2*imag(mm(3,1));
s(9)=2*imag(mm(1,2));
s=(sqrt(3)/2)*s./s(1);
s(1)=1;

% s=real(s);
end