function[Mdep]=second_moment_matrix_iso(major_axis,minor_axis,a,b,c)
%.055-----degpol 0.8
%0.1--------0.7
%0.16-------0.6
%0.25--------0.5
%0.43-------0.4
L1=a;
L2=b;
L3=c;
E1=major_axis+1i*minor_axis;
M1=E1'*E1; M3=eye(3,3); 
Mdep=(L1-L2)*(M1)+L3*(M3); Mdep=Mdep/trace((Mdep));

l1=L1/(L1+L2+L3);
l2=L2/(L1+L2+L3);
l3=L3/(L1+L2+L3);
P=sqrt(((l1-l2)^2+(l1-l3)^2+(l2-l3)^2)/2);


end