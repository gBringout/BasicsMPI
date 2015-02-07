% This example demonstrated the basic of MPI signal generation using the
% example described in the Ph.D. Thesis from S. Biederer untitled:
% "Entwicklung eines Spektrometers zur Analyse superparamagnetischer Eisenoxid-Nanopartikel f�r Magnetic-Particle-Imaging"
% Gael Bringout - Oct 2014

close all
clear all

addpath(genpath(fullfile('.')))

    
mu0         = 4*pi*1e-7; % permeability of the air.
%particle
c           = 0.03*0.500*1000; % [mol(Fe)/m^3] Using Resovist, with the assumption that only 3% of the iron is usefull to the signal
v           = 10*10^-9; %[m^3] volume of the test sample (10 �l) (page 99 of the thesis)
d           = 30*10^-9; % [m] diameter of the particle. From the thesis. P98.
Ms          = 0.6/mu0; %[A/m] Magnetisation at saturation from Magnetit. From the thesis. P54
temperature = 310;
%signal generation and procesing
Fs          = 10*10^6;%2.5e7; % [Hz] sampling frequency
fStart      = 2; % [-] staring frequency components to be displayed
fStop       = 800; % [-] Last frequency components to be displayed. Please respect the shanon theorem :)
dt          = 1/(50*10^6); %delta used to calculate the derivative. Used to compensate the "small" sampling frequencies

% MPI signal generation
f           = 25000; % [Hz] Frequency of the signal
nbrPeriod   = 4; % nbr of signal period
time        = 0:1/Fs:nbrPeriod/f-1/Fs; % [s] time vector. Remove the last point as we start at zero.
N           = length(time); %Number of time points
H0          = 31831/2;% [A/m] Drive field amplitude (31831 A/m to generate 40 mT)
B0          = H0*mu0; % [T] - H = B/mu

% MPI scanner properties
% According to Knopp "introduction into MPI" pdf 35 equ 2.36
%we assume that the sensibility is constante in this voxel
% the sensibilte is there define as the field generated by a 1A current
s           = 1366*mu0; % [T/A] From the thesis. P98. Eq. 4.67

% Noise model
kB          = 1.380650424e-23; % [J/K] Boltzmann constant 
deltaF      = 10*10^6; % [Hz] according to Weizenecker 2007 - A simulation study...
Rp          = 0.0001*185 *10^-3; % [Ohm]  according to Weizenecker 2007 - A simulation study...
Scaling     = 1000; % To easily scale the maximal voltage induced by the noise
%% Particle Magnetization in function of the applied magnetic field

nbrPoint = 1000;
a=zeros(nbrPoint,1);
mu0 = 4*pi*1e-7; % permeability of the air.
Hpart=linspace(-31831,31831,nbrPoint);
Bpart = Hpart.*mu0;
Mpart=zeros(nbrPoint,1);
for i=1:nbrPoint,
  [Mpart(i),~,~,a(i)] = langevinParticle4( Bpart(i),0,0,abs(Bpart(i)), d,Ms,temperature,c,v);

end
figure
subplot(2,3,1)
hold off
plot(Bpart, Mpart, 'k');
%plot(a, Mpart, 'k');
%xlim([-max(Bpart) max(Bpart)])
xlabel('B / (T)')
ylabel('M / (A/m)')
title('Particle magnetisation curve');

%% Drive field Amplitude in function of the time
H=H0*sin(2*pi*f*time);
Hdt=H0*sin(2*pi*f*(time-dt));
B = H.*mu0;
Bdt = Hdt.*mu0;

subplot(2,3,4)
plot(time,B,'k');
xlabel('t / s')
ylabel('B / (T)')
title(sprintf('Drive field amplitude: %0.2g mT peak',B0*1000));

%% Particle Magnetization in function of the time
M=zeros(N,1);
Mdt=zeros(N,1);
for i=1:N,
  [M(i),~,~,~] = langevinParticle4( B(i),0,0,abs(B(i)), d,Ms,temperature,c,v);
  [Mdt(i),~,~,~] = langevinParticle4( Bdt(i),0,0,abs(Bdt(i)), d,Ms,temperature,c,v);
end

subplot(2,3,2)
plot(time, M, 'k');
xlabel('t / s')
ylabel('M / (A/m)')
title('Particle magnetisation');

Umax = Scaling*sqrt(4*kB*temperature*deltaF*Rp);% [V]
u=zeros(N,1);

for i=1:N,
	u(i) = s*(M(i)-Mdt(i))/dt + normrnd(0,Umax);
end

subplot(2,3,3)
plot(time, u, 'k');
xlabel('t / s')
ylabel('u / V')
title('Time signal');

%% Spectrum

subplot(2,3,6)
periodogram = fft(u)/size(u,1);
uhat = abs((2*periodogram(1:(end-1)/2+1)));
Fs = 1/(time(2)-time(1));
freq = Fs/2*linspace(0,1,length(u)/2+1);
plot(freq(fStart:fStop)/f,real(uhat(fStart:fStop)), 'o');
title('Spectrum - # of harmonic');
xlim([ 0 50])

% filter the results
%figure; semilogy(abs(periodogram));
i1 = find(freq == f);
i2 = find(freq == 2*f);
i3 = find(freq == 3*f);

AS = abs(fft(u));
AS(i1) = AS(i1)/400;
AS(size(u,1)-i1+2) = AS(i1);
AS(i2) = AS(i2)/200;
AS(size(u,1)-i2+2) = AS(i2);
AS(i3) = AS(i3);
AS(size(u,1)-i3+2) = AS(i3);

%ybis = ifft(AS);
%figure;plot(time,real(ybis));