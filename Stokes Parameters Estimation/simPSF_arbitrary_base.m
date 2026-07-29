function[denoIk]=simPSF_arbitrary_base(bs,base,denoIk)

for k=1:9
    denoIk=denoIk+bs(1,k).*base(:,:,k);
end

end
