function [P]=polarization_degree(a,b,c)
L1=a;
L2=b;
L3=c;
l1=L1./(L1+L2+L3);
l2=L2./(L1+L2+L3);
l3=L3./(L1+L2+L3);
P=sqrt(((l1-l2).^2+(l1-l3).^2+(l2-l3).^2)/2);
end