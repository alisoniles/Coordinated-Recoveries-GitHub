function [psi,pdi] = phasefcimf(imf1,imf2,imfsize)

%imfsize=size(imf1);
%imfsize=imfsize(2);
%imfsize=imfsize-1;

psi=zeros(imfsize,1);
pdi=zeros(imfsize,1);
for i=1:imfsize
    h1=hilbert(imf1(:,i));
    h2=hilbert(imf2(:,i));
    a1=angle(h1);
    a2=angle(h2);
    [psi(i),pdi(i)]=pscoherence(unwrap(a1)-unwrap(a2));
end

