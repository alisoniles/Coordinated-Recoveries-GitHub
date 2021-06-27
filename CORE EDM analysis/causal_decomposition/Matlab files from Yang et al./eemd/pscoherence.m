function [pci,pdi] = pscoherence(phasediff);

% transform phase difference to a complex plane with radius of one
z=exp((phasediff)*sqrt(-1));

% phase coherence
pcv=sum(z)/length(z);
pci=abs(pcv);

% phase shift (in degree)
pdi=angle(pcv)/pi*180;


