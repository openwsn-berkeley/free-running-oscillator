% Matlab code to model the White and Flicker of frequency noise of a free running oscillator.
% Author: Dr. Osama Khan
% Email: oukhan@berkeley.edu
% Company: SWARM lab UC Berkeley
% Date: 1/19/2016
% Version: 1.0


Fc=2447e6; % Oscillator Frequency
Fs=6e9;    % Simulation Sampling Frequency
time=0:1/Fs:0.1e-3-(1/Fs);% Time vector
 
% Variable initialization
N=numel(time);
xPN=zeros(1,N);
alphaN=zeros(1,N);
phiN=zeros(1,N);
h=zeros(1,N);
beta=-3;
cvco=0.5e-15;% White Frequency noise coefficient
cvcofn=0.5e-10;% Flicker of Frequency noise coefficient
 
deltaf=Fs/N;   % FFT bin size
fundindex=Fc/deltaf+1;% FFT bin number for the fundamental tone
freq=-Fs/2:deltaf:Fs/2-deltaf; % FFT Frequency vector
 
% Computing Impulse response coefficients for Flicker of Frequency noise
h(1)=1;
for i=2:N
h(i)=(i-2-beta/2)*h(i-1)/(i-1);
end
 
w=randn(1,N);% Generating AWGN noise
% The convolution operation below can be speeded up by multiplying
% the Fourier transform in the frequency domain and taking the 
% inverse FFT. Care must be taken to perform the linear convolution instead
% of the circular convolution.
fconvt=conv(h,w);% Generating 1/f noise. 
 
 
alphaN(1)=sqrt(cvco)*sqrt((1/Fs))*randn(1,1);
phiN(1)=2*pi*Fc*time(1)+2*pi*Fc*alphaN(1);
xPN(1)=cos(2*pi*Fc*time(1)+2*pi*Fc*alphaN(1));
 
% Discrete time simulation of Oscillator with PN
for i=2:N
    alphaN(i)=alphaN(i-1)+sqrt(cvco)*sqrt((1/Fs))*randn(1,1);
end
 
phiN=2*pi*Fc*time(1:N)+2*pi*Fc*alphaN(1:N)+2*pi*Fc*sqrt(2*pi*cvcofn)*(1/Fs)*fconvt(1:N);
xPN=cos(phiN(1:N));
 
% Compute FFT
xfPN=(1/Fs).*fft(xPN);
xfPNpos=2*xfPN(1:N/2);
xPNptot=sum(abs(xfPNpos).^2);
pxfPN=10*log10(abs(xfPNpos).^2/xPNptot);  % dBc per FFT bin size in Hz.
 
% Plot the normalized RF spectral density in units of dBc per Hz
figure;
semilogx(freq(N/2+fundindex+1:end)-Fc,pxfPN(fundindex+1:end)-10*log10(deltaf),'r','LineWidth',2);
grid on;
hold on;
xlabel('Frequency (Hz)')
ylabel('RF spectral density (dBc/Hz)')
legend('Matlab DT PN model White + Flicker','AFS')
set(gca,'FontSize',16)
