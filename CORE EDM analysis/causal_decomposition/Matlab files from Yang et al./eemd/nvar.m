function nv=nvar(imf)

% calculate variance
v=var(imf);

% calculate normalized variance
nv=v(1:end-1)/sum(v(1:end-1));
%nv=v(1:end)/sum(v(1:end));
