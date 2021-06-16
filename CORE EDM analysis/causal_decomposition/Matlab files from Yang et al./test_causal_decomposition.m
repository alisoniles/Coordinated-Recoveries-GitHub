
clearvars
close all


load ecosystem_data.mat

plot_pairedimfs(time_DIDIPARA,DIDI,PARA,0,1000,[]);

causal_matrix = causal_decomposition(DIDI,PARA,0,1000);



%Hilbert Transform
%To build a CFC measure, the first step is to find the phase or amplitude envelope of the signal. We can do this by % Generate sinusoid signal
dt = 0.001; % sampling interval [s]
t = dt:dt:1; % time axis [s]
f1 = 2.0; % freq of sinusoid [Hz]
phi0 = 0.0; % initial phase of sinusoid 
d = sin(2.0*pi*t*f1 + phi0);

% Compute analytic signal
dA = hilbert(d);
% The built-in "hilbert" function in MATLAB returns the analytic signal. We can plot the original signal, and the real and imaginary parts of the signal.
figure; clf() % clf = clear current figure window subplot(4,1,1) % 4x1 plot, 1st plot
plot(d); ylabel('Data');
subplot(4,1,2) % 4x1 plot, 2nd plot 
plot(real(dA));
hold on; 
plot(imag(dA), 'r');
hold off;
ylabel('Real (blue), Imag (red)'); axis tight
% Note that the imaginary part of the analytic signal is the real part of the analytic signal shifted by 90 degrees (pi/2)
% Compute phase of analytic signal
phi = angle(dA);
% Compute amplitude envelope
amp = abs(dA);
% Plot results
subplot(4,1,3); plot(phi); ylabel('Angle'); axis tight subplot(4,1,4); plot(amp); ylabel('Amplitude'); axis tight