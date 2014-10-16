clear all
close all

addpath(genpath('.'))
addpath(genpath('..\SphericalHarmonics\'))
addpath(genpath('..\ScannerDesign\'))

%% Load the fields and prepae the data
disp('Loading and pre-calculation of datas')
load('IdealFFL.mat');

%% Define the parameters
disp('calculate the frequency etc.')
%we asumme that we have a squarre voxel with size dx
calculation.dx = 1*10^-3/2; % [m] size of a voxelin the x direction
calculation.dy = 1*10^-3/2; % [m] size of a voxelin the y direction
calculation.dz = 1*10^-3/2; % [m] size of a voxelin the z direction

system.frequencyLineRotation = 100;% The corresponding line rotation frequency
system.frequencyQuadrupol = 2*system.frequencyLineRotation;
system.frequencyDrive = 25000; % Hz

system.nbrPeriods = 2; % To adapte the resolution of the spectra
system.samplingFrequency = 1e6; %[Hz] Sampling Frequency of the acquisition Board.
%Maximal and minimal frequencies in the spectrum
%system.fSMmin = 330*1000; %because of the filters
%system.fSMmax = 3e6; % to keep the data usage low
system.fSMmin = 57*1000; %because of the filters
system.fSMmax = 0.5e6; % to keep the data usage low

system.timeLength = system.nbrPeriods/(2*system.frequencyLineRotation);
system.numberOfTimePoints  = system.timeLength*system.samplingFrequency;
calculation.time = linspace(0,system.timeLength,system.numberOfTimePoints+1);
calculation.time = calculation.time(1:end-1);% Remove the last point as we start at zero.

calculation.dt= 1/25*10^-6;
calculation.deltaf = 1/system.timeLength;

calculation.numberOfFrequencies = system.timeLength*system.samplingFrequency/2.+1; % maximal number of frequencies
system.nbrPointPerDrivePeriod = system.samplingFrequency/system.frequencyDrive;
system.freq = system.samplingFrequency/2*linspace(0,1,calculation.numberOfFrequencies);

%general
calculation.mu0 = 4*pi*1e-7; % permeability of the air.

% Particle parameters
system.diameterPart = 30*10^-9;
system.MsatSinglePart = 0.6/calculation.mu0;
system.Tempe = 310;
system.concentrationPartiMax = 0.03*0.5*1000;% [mol/m^3] 3% of "good" particle,1000 is for l to m^�,/10 for a dilution
system.volumeSample = calculation.dx*calculation.dy*calculation.dz; % [m^3] volume of a voxel

% We use the noise model Weiznecker 2007
calculation.kB  = 1.380650424e-23;
noise.Rp = 185*10^-3;
noise.T = 310;
noise.deltaF = 10*10^6; % [Hz] according to Weizenecker 2007 - A simulation study...
noise.maxAmplitude = sqrt(4*calculation.kB*noise.T*noise.deltaF*noise.Rp)/1;
noise.maxAmplitudeSM = sqrt(4*calculation.kB*noise.T*noise.deltaF*noise.Rp)/10; % For the System matrix
system.SNRLimits = 6;
system.maxIterationReco = 10;

% defining the geometry of the Field of View
system.x = -0.0055:calculation.dx:0.005; %The resolution have to be limited in order to be able to reconstruct
system.y = -0.0065:calculation.dy:0.006;
system.z = 0:calculation.dz:0;

%Size of the matrix
system.sizeX = size(system.x,2);
system.sizeY = size(system.y,2);
system.sizeZ = size(system.z,2);


%% Calculate the fields
disp('Reconstruct the fields from SHC')
tic
Selection_Z.B = RebuildField7(Selection_Z.bc,Selection_Z.bs,Selection_Z.rhoReference,system.x,system.y,system.z,'sch');
Quadru_0.B = RebuildField7(Quadru_0.bc,Quadru_0.bs,Quadru_0.rhoReference,system.x,system.y,system.z,'sch');
Quadru_45.B = RebuildField7(Quadru_45.bc,Quadru_45.bs,Quadru_45.rhoReference,system.x,system.y,system.z,'sch');
Drive_X.B = RebuildField7(Drive_X.bc,Drive_X.bs,Drive_X.rhoReference,system.x,system.y,system.z,'sch');
Drive_Y.B  = RebuildField7(Drive_Y.bc,Drive_Y.bs,Drive_Y.rhoReference,system.x,system.y,system.z,'sch');

% filter the field value outside a diameter of 0.5 m
% (to simulate cylindrical coils)
for i=1:system.sizeX
    for j=1:system.sizeY
        if sqrt(system.x(i)^2+system.y(j)^2)>0.25
            Selection_Z.B(:,i,j)    = 10^-16;
            Quadru_0.B(:,i,j)       = 10^-16;
            Quadru_45.B(:,i,j)      = 10^-16;
            Drive_X.B(:,i,j)        = 10^-16;
            Drive_Y.B(:,i,j)        = 10^-16;
        end
    end
end
toc

% Make the phantom
% It has to have the same size as the fields
phantom.shape = [system.sizeX system.sizeY system.sizeZ]; % Number of pixel per direction
%phantom.shape = [1820 1860 1]; % Number of pixel per direction
phantom.size = [2*max(system.x) 2*max(system.y) 0.0001]; % size of thr FoV
%phantom.diameters = [0.002, 0.004, 0.006, 0.008]; % [m] diameter of the circles 
phantom.particleDiameter = system.diameterPart; %Best resolution with 50nm, Worst with 30nm
phantom.concentrationPartiMax =  system.concentrationPartiMax; % [mol/l]undiluted Resovist: 5.0e-4;
%phantom.shapeScaled = createResolutionPhantomGael(phantom.shape, phantom.size ,phantom.diameters);
phantom.shapeScaled = createResolutionPhantomGael4(phantom.shape, 8);
%phantom.shapeScaled = createResolutionPhantomGael3(phantom.shape);

%Make sure of the scaling of the phantom
phantom.shapeScaled = phantom.shapeScaled/max(phantom.shapeScaled(:));
%figure;imagesc(system.x,system.y,phantom.shapeScaled)

% the sensibilte is here define as the field generated by a 1A current
% this is why we have to save the bc and bs for a unit current
system.s1 = Drive_X.B;
calculation.s1x = system.s1(1,:);
calculation.s1y = system.s1(2,:);
calculation.s1z = system.s1(3,:);

system.s2 = Drive_Y.B;
calculation.s2x = system.s2(1,:);
calculation.s2y = system.s2(2,:);
calculation.s2z = system.s2(3,:);

%It is not important to have all the particle in the field of view
% as long as the field have been deletec outside the FoV :)
for i=1:system.sizeX
    for j=1:system.sizeY
        if sqrt(system.x(i)^2+system.y(j)^2)>0.25
            phantom.shapeScaled(i,j) = 0;
        end
    end
end

clear('i','j')

%%
% We have to calculate the derivative,
% We thus need to now the field after a dt time
% We then have the H matric where the field at time t is stored
% And the Hdt matrix where the field at time t+dt is stored

disp('Generating the time varying fields')
tic
c1 = cos(2*pi*system.frequencyQuadrupol*calculation.time(:)');
c2 = sin(2*pi*system.frequencyQuadrupol*calculation.time(:)');
c3 = cos(0.5*2*pi*system.frequencyQuadrupol*calculation.time(:)');
c4 = sin(0.5*2*pi*system.frequencyQuadrupol*calculation.time(:)');
c5 = sin(2*pi*system.frequencyDrive*calculation.time(:)');

c1_dt = cos(2*pi*system.frequencyQuadrupol*(calculation.time(:)'+calculation.dt));
c2_dt = sin(2*pi*system.frequencyQuadrupol*(calculation.time(:)'+calculation.dt));
c3_dt = cos(0.5*2*pi*system.frequencyQuadrupol*(calculation.time(:)'+calculation.dt));
c4_dt = sin(0.5*2*pi*system.frequencyQuadrupol*(calculation.time(:)'+calculation.dt));
c5_dt = sin(2*pi*system.frequencyDrive*(calculation.time(:)'+calculation.dt));

system.coefSelection_Z = ones(1,system.numberOfTimePoints);
system.coefQuadru_0 =   c1;
system.coefQuadru_45 =  c2;
system.coefDrive_X =    c5.*c3;
system.coefDrive_Y =    c5.*c4;

system.coefSelection_Z_dt = ones(1,system.numberOfTimePoints);
system.coefQuadru_0_dt =   c1_dt;
system.coefQuadru_45_dt =  c2_dt;
system.coefDrive_X_dt =    c5_dt.*c3_dt;
system.coefDrive_Y_dt =    c5_dt.*c4_dt;

Bx =system.coefSelection_Z'*   Selection_Z.current*Selection_Z.B(1,:)+...
    system.coefQuadru_0'   *	Quadru_0.current*Quadru_0.B(1,:)+...
    system.coefQuadru_45'  *	Quadru_45.current*Quadru_45.B(1,:)+...
    system.coefDrive_X'    *	Drive_X.current*Drive_X.B(1,:)+...
    system.coefDrive_Y'    *	Drive_Y.current*Drive_Y.B(1,:);

By =system.coefSelection_Z'*   Selection_Z.current*Selection_Z.B(2,:)+...
    system.coefQuadru_0'   *	Quadru_0.current*Quadru_0.B(2,:)+...
    system.coefQuadru_45'  *	Quadru_45.current*Quadru_45.B(2,:)+...
    system.coefDrive_X'    *	Drive_X.current*Drive_X.B(2,:)+...
    system.coefDrive_Y'    *	Drive_Y.current*Drive_Y.B(2,:);

Bz =system.coefSelection_Z'*   Selection_Z.current*Selection_Z.B(3,:)+...
    system.coefQuadru_0'   *	Quadru_0.current*Quadru_0.B(3,:)+...
    system.coefQuadru_45'  *	Quadru_45.current*Quadru_45.B(3,:)+...
    system.coefDrive_X'    *	Drive_X.current*Drive_X.B(3,:)+...
    system.coefDrive_Y'    *	Drive_Y.current*Drive_Y.B(3,:);

Babs = sqrt(Bx.^2+By.^2+Bz.^2);

Bx_dt = system.coefSelection_Z_dt'*   Selection_Z.current*Selection_Z.B(1,:)+...
        system.coefQuadru_0_dt'   *	Quadru_0.current*Quadru_0.B(1,:)+...
        system.coefQuadru_45_dt'  *	Quadru_45.current*Quadru_45.B(1,:)+...
        system.coefDrive_X_dt'    *	Drive_X.current*Drive_X.B(1,:)+...
        system.coefDrive_Y_dt'    *	Drive_Y.current*Drive_Y.B(1,:);

By_dt = system.coefSelection_Z_dt'*   Selection_Z.current*Selection_Z.B(2,:)+...
        system.coefQuadru_0_dt'   *	Quadru_0.current*Quadru_0.B(2,:)+...
        system.coefQuadru_45_dt'  *	Quadru_45.current*Quadru_45.B(2,:)+...
        system.coefDrive_X_dt'    *	Drive_X.current*Drive_X.B(2,:)+...
        system.coefDrive_Y_dt'    *	Drive_Y.current*Drive_Y.B(2,:);

Bz_dt = system.coefSelection_Z_dt'*   Selection_Z.current*Selection_Z.B(3,:)+...
        system.coefQuadru_0_dt'   *	Quadru_0.current*Quadru_0.B(3,:)+...
        system.coefQuadru_45_dt'  *	Quadru_45.current*Quadru_45.B(3,:)+...
        system.coefDrive_X_dt'    *	Drive_X.current*Drive_X.B(3,:)+...
        system.coefDrive_Y_dt'    *	Drive_Y.current*Drive_Y.B(3,:);

Babs_dt = sqrt(Bx_dt.^2+By_dt.^2+Bz_dt.^2);
    
clear('c1','c2','c3','c4','c5','c1_dt','c2_dt','c3_dt','c4_dt','c5_dt')
fprintf('Time taken %2.0f s.\n', toc)
%% Display the fields
% figure
% threshold = 3*10^-3;
% for i=1:system.nbrPointPerDrivePeriod/20:system.numberOfTimePoints
%     image = reshape(Babs(i,:),[system.sizeX,system.sizeY]);
%     imagesc(system.x,system.y,image);
%     xlabel('x axis /m')
%     ylabel('y axis /m')
%     axis square
%     caxis([-threshold threshold]);
%     pause(1/25)
%     %title(sprintf('Magnetic field density /T. Time %3.3f ms',calculation.time(i)*1000))
% end

%% Here we calculate the needed magnetization for the system matrix.
% it is in fact a phantom with a concentration of 1 in any position
% Later, we will just have to multiply this magnetization matrix by the
% real concentration of the phantom, to get the measurements
% Note that in order to optimize the needed RAM, we put the result into B
% and then rename it M

disp('Calculting the time varying magnetization')
tic
for i=1:system.numberOfTimePoints,
    [Bx(i,:),By(i,:),Bz(i,:)] = langevinParticle4(Bx(i,:),By(i,:),Bz(i,:),Babs(i,:),phantom.particleDiameter,system.MsatSinglePart,system.Tempe,phantom.concentrationPartiMax,system.volumeSample);
    [Bx_dt(i,:),By_dt(i,:),Bz_dt(i,:)] = langevinParticle4(Bx_dt(i,:),By_dt(i,:),Bz_dt(i,:),Babs_dt(i,:),phantom.particleDiameter,system.MsatSinglePart,system.Tempe,phantom.concentrationPartiMax,system.volumeSample);
end
Mx = Bx;
My = By;
Mz = Bz;
Mx_dt = Bx_dt;
My_dt = By_dt;
Mz_dt = Bz_dt;
clear('Bx','By','Bz','Babs','Bx_dt','By_dt','Bz_dt','Babs_dt')
fprintf('Time taken %2.0f s.\n', toc)

%% Calculate the system matrix
% Here we have the derivative of the Magnetization
% If we concidere that the voltage at the begining of the experiment is nul
% Then the first line are null
% and the derivative is calculated with M(i-1)-M(i)
disp('Calculating the system matrix')
tic
u1=zeros(system.numberOfTimePoints,system.sizeX*system.sizeY*system.sizeZ);
u2=zeros(system.numberOfTimePoints,system.sizeX*system.sizeY*system.sizeZ);
for i=1:system.numberOfTimePoints
    u1(i,:) = ((calculation.s1x.*(Mx_dt(i,:)-Mx(i,:))/calculation.dt)+(calculation.s1y.*(My_dt(i,:)-My(i,:))/calculation.dt)+(calculation.s1z.*(Mz_dt(i,:)-Mz(i,:))/calculation.dt)) + normrnd(0,noise.maxAmplitudeSM);
    u2(i,:) = ((calculation.s2x.*(Mx_dt(i,:)-Mx(i,:))/calculation.dt)+(calculation.s2y.*(My_dt(i,:)-My(i,:))/calculation.dt)+(calculation.s2z.*(Mz_dt(i,:)-Mz(i,:))/calculation.dt)) + normrnd(0,noise.maxAmplitudeSM);
end
%u1 = u1+u2;
results.SM1 = fft(u1)'; % we take the FFT
results.SM1 = results.SM1(:,1:calculation.numberOfFrequencies);% and use the one sided part
results.SM2 = fft(u2)'; % we take the FFT
results.SM2 = results.SM2(:,1:calculation.numberOfFrequencies);% and use the one sided part

%SM = SM(:,1:floor(system.numberOfTimePoints/2+1)); % and use the one sided part
% filter directly the frequencies
% finding the index
[~,Imin] = find(abs(system.freq-system.fSMmin) < calculation.deltaf,1,'first');
[~,Imax] = find(abs(system.freq-system.fSMmax) < calculation.deltaf,1,'first');
% We now remove the useless frequencies
results.SM1 = results.SM1(:,1:Imax);
results.SM2 = results.SM2(:,1:Imax);

%save('SM_FFL6_HR_361x363.mat','SM','-v7.3')

clear('Imax','Imin','u1','u2','i','j');
fprintf('Time taken %2.0f s.\n', toc)

%% Multiply by the phantom
disp('Applaying the phantom');
tic
for i=1:system.numberOfTimePoints,
    Mx(i,:) = Mx(i,:).*phantom.shapeScaled(:)';
    My(i,:) = My(i,:).*phantom.shapeScaled(:)';
    Mz(i,:) = Mz(i,:).*phantom.shapeScaled(:)';
    Mx_dt(i,:) = Mx_dt(i,:).*phantom.shapeScaled(:)';
    My_dt(i,:) = My_dt(i,:).*phantom.shapeScaled(:)';
    Mz_dt(i,:) = Mz_dt(i,:).*phantom.shapeScaled(:)';
end

clear('i');
fprintf('Time taken %2.0f s.\n', toc)

%% Calculating the signal
% Here we have the derivative of the Magnetization
% If we concidere that the voltage at the begining of the experiment is nul
% Then the first line are null
% and the derivative is calculated with M(i-1)-M(i)
disp('Calculating the signal')
tic

results.u1_2=zeros(system.numberOfTimePoints,1);
results.u2_2=zeros(system.numberOfTimePoints,1);
for i=1:system.numberOfTimePoints
    results.u1_2(i) = sum(((calculation.s1x.*(Mx_dt(i,:)-Mx(i,:))/calculation.dt)+(calculation.s1y.*(My_dt(i,:)-My(i,:))/calculation.dt)+(calculation.s1z.*(Mz_dt(i,:)-Mz(i,:))/calculation.dt)))+ normrnd(0,noise.maxAmplitude);
    results.u2_2(i) = sum(((calculation.s2x.*(Mx_dt(i,:)-Mx(i,:))/calculation.dt)+(calculation.s2y.*(My_dt(i,:)-My(i,:))/calculation.dt)+(calculation.s2z.*(Mz_dt(i,:)-Mz(i,:))/calculation.dt)))+ normrnd(0,noise.maxAmplitude);
end
clear('Mx','My','Mz','Mx_dt','My_dt','Mz_dt');

%results.u1_2 = results.u1_2+results.u2_2;
%  Extracting the spectrum of the particle signal
results.signal1FFT = fft(results.u1_2)'; %We store in the line each point, in the column the time point
results.signal2FFT = fft(results.u2_2)';
results.signal1FFT_oneSided = results.signal1FFT(1:calculation.numberOfFrequencies);% We now remove the folded part of the spectrum
results.signal2FFT_oneSided = results.signal2FFT(1:calculation.numberOfFrequencies);
results.signal1AbsFFT_oneSided = abs(results.signal1FFT_oneSided);
results.signal2AbsFFT_oneSided = abs(results.signal2FFT_oneSided);
calculation.S1 = 1; %as we are using a rectangular windows
results.signal1PowerSpectrum = 2*results.signal1AbsFFT_oneSided.^2/calculation.S1^2; % According to report GH_FFT from G. Heinzel
results.signal2PowerSpectrum = 2*results.signal2AbsFFT_oneSided.^2/calculation.S1^2; % According to report GH_FFT from G. Heinzel
results.signal1AmplitudeSpectrum = sqrt(results.signal1PowerSpectrum);
results.signal2AmplitudeSpectrum = sqrt(results.signal2PowerSpectrum);

% Calculation of the Power on the noise signal
noise.uNoise = zeros(system.numberOfTimePoints,1);
for i=1:system.numberOfTimePoints
	noise.uNoise(i) = normrnd(0,noise.maxAmplitude);
end
noise.noiseFFT = fft(noise.uNoise)';   
noise.noiseFFT_oneSided= noise.noiseFFT(1:calculation.numberOfFrequencies); % We now remove the folded part of the spectrum and multiply by 2 to keep the same energy
noise.noiseAbsFFT_oneSided = abs(noise.noiseFFT_oneSided); % Is this really an energy?
calculation.S1 = 1; %as we are using a rectangular windows
noise.noisePowerSpectrum = 2*noise.noiseAbsFFT_oneSided.^2/calculation.S1^2; % According to report GH_FFT from G. Heinzel
noise.noiseAmplitudeSpectrum = sqrt(noise.noisePowerSpectrum);
%calculation.S2 = 1; %as we are using a rectangular windows
%noiseAmplitudeSpectralDensity = 2*results.noiseAbsFFT_oneSided.^2/(calculation.S2^2*samplingFrequency);
noise.stdNoiseAmplitudeSpectrum = std(noise.noiseAmplitudeSpectrum);
%noise.meanNoiseAmplitudeSpectrum = mean(noise.noiseAmplitudeSpectrum);
%maxNoiseAmplitudeSpectrum = max(results.noiseAmplitudeSpectrum);

% Calculating the SNR
results.noiseLevel = noise.stdNoiseAmplitudeSpectrum;
results.nbrGoodFrequency1 = 1;
results.nbrGoodFrequency2 = 1;
results.signal1SNR = zeros(1,size(noise.noiseFFT_oneSided,2));
results.signal2SNR = zeros(1,size(noise.noiseFFT_oneSided,2));
%count the number of frequency
for i=2:calculation.numberOfFrequencies
    results.signal1SNR(i) = results.signal1AmplitudeSpectrum(i)/results.noiseLevel;
    results.signal2SNR(i) = results.signal2AmplitudeSpectrum(i)/results.noiseLevel;
    if  results.signal1SNR(i) > system.SNRLimits && system.freq(i) > system.fSMmin && system.freq(i) < system.fSMmax
        results.nbrGoodFrequency1 = results.nbrGoodFrequency1+1;
    end
    if  results.signal2SNR(i) > system.SNRLimits && system.freq(i) > system.fSMmin && system.freq(i) < system.fSMmax
        results.nbrGoodFrequency2 = results.nbrGoodFrequency2+1;
    end
end

% Filtering signals and SM
index = 1;
results.tSNR1 = zeros(1,results.nbrGoodFrequency1-1);
%goodFrequency = zeros(1,nbrGoodFrequency-1);
%tSignalAmplitudeSpectrum = zeros(1,nbrGoodFrequency-1);
results.tFreq1 = zeros(1,results.nbrGoodFrequency1-1);
results.tSM1 = zeros(size(results.SM1,1),results.nbrGoodFrequency1-1);
%tSignalPower = zeros(1,nbrGoodFrequency-1); % The signal with zero where it is not taken into account
results.tSignal1FFT_oneSided = zeros(1,results.nbrGoodFrequency1-1);
for i=2:size(results.signal1FFT_oneSided,2)
    if results.signal1SNR(i) > system.SNRLimits && system.freq(i) > system.fSMmin && system.freq(i) < system.fSMmax
        results.tSNR1(index) = results.signal1SNR(i); % used in a figure
        results.tSignal1FFT_oneSided(index) = results.signal1FFT_oneSided(i); % used in reco
        %tSignalAmplitudeSpectrum(index) = signalAmplitudeSpectrum(i);
        %tSignalPower(index) = signalPowerSpectrum(i);
        results.tSM1(:,index) = results.SM1(:,i); % used in reco
        results.tFreq1(index) = system.freq(i); % used in a figure
        %goodFrequency(index) = 1;
        index = index+1;
    end
end

% second signal
index = 1;
results.tSNR2 = zeros(1,results.nbrGoodFrequency2-1);
%goodFrequency = zeros(1,nbrGoodFrequency-1);
%tSignalAmplitudeSpectrum = zeros(1,nbrGoodFrequency-1);
results.tFreq2 = zeros(1,results.nbrGoodFrequency2-1);
results.tSM2 = zeros(size(results.SM2,1),results.nbrGoodFrequency2-1);
%tSignalPower = zeros(1,nbrGoodFrequency-1); % The signal with zero where it is not taken into account
results.tSignal2FFT_oneSided = zeros(1,results.nbrGoodFrequency2-1);
for i=2:calculation.numberOfFrequencies
    if results.signal2SNR(i) > system.SNRLimits && system.freq(i) > system.fSMmin && system.freq(i) < system.fSMmax
        results.tSNR2(index) = results.signal2SNR(i); % used in a figure
        results.tSignal2FFT_oneSided(index) = results.signal2FFT_oneSided(i); % used in reco
        %tSignalAmplitudeSpectrum(index) = signalAmplitudeSpectrum(i);
        %tSignalPower(index) = signalPowerSpectrum(i);
        results.tSM2(:,index) = results.SM2(:,i); % used in reco
        results.tFreq2(index) = system.freq(i); % used in a figure
        %goodFrequency(index) = 1;
        index = index+1;
    end
end

clear('index','i','j');
fprintf('Time taken %2.0f s.\n', toc)

%% Reconstrucing with partial signal (SNR thresholding)

disp('Reconstructing with partial signal')

tic
%[results.X,results.rho,results.eta] = artGael([results.tSM1]',[results.tSignal1FFT_oneSided]',50);
%[results.X,results.rho,results.eta] = artGael([results.tSM2]',[results.tSignal2FFT_oneSided]',50);
%[results.X,results.rho,results.eta] = artGael([results.tSM1+results.tSM2]',[results.tSignal1FFT_oneSided+results.tSignal1FFT_oneSided]',50);
[results.X,~,~] = artGael([results.tSM1,results.tSM2]',[results.tSignal1FFT_oneSided,results.tSignal2FFT_oneSided]',system.maxIterationReco);

results.errorEstimate = zeros(1,system.maxIterationReco);
for i=1:system.maxIterationReco
    res = reshape(results.X(:,i),[system.sizeX,system.sizeY]);
    results.errorEstimate(i) = sum(sum(sqrt((phantom.shapeScaled-res).^2)));
end

figure
for i=1:system.maxIterationReco
    subplot(2,2,1)
    imagesc(phantom.shapeScaled)
    axis square
    subplot(2,2,2)
    res = reshape(results.X(:,i),[system.sizeX,system.sizeY]);
    imagesc(system.x,system.y,real(res));
    colormap('gray')
    axis square
    title(sprintf('i=%i',i));
    pause(1/25)
end
subplot(2,2,[3 4])
plot(results.errorEstimate)

clear('res')
fprintf('Time taken %2.0f s.\n', toc)

%% Figure
disp('display the results')
figure('Name','Signal')

subplot(4,1,1)
hold all
plot(calculation.time*1000,system.coefDrive_X*Drive_X.current)
plot(calculation.time*1000,system.coefDrive_Y*Drive_Y.current)
xlabel('Time / ms')
ylabel('Current amplitude / A')
title('Drive fields current')

subplot(4,1,2)
plot(calculation.time*1000,system.coefQuadru_0*Quadru_0.current)
hold all
plot(calculation.time*1000,system.coefQuadru_45*Quadru_45.current)
xlabel('Time / ms')
ylabel('Current amplitude / A')
title('Quadrupoles current')

subplot(4,1,3)
plot(calculation.time*1000,results.u1_2+results.u2_2)
xlabel('Time / ms')
ylabel('Voltage amplitude / V')
title('Sum of the induced voltage in the x and y canal.')

subplot(4,1,4)
stem(system.freq/system.frequencyDrive,results.signal1SNR,'Marker','None')
set(gca,'yscal','log')
hold all
stem(results.tFreq1/system.frequencyDrive,results.tSNR1,'Marker','None')
set(gca,'yscal','log')
xlabel('# Harmonic')
ylabel('SNR')
title('SNR (based on the amplitude spectrum)')
xlim([0 21])
ylim([1 10^2]);%std(signalSNR)])

% Calculated the amplitude of the noise associted with the SM
noise.uNoise = zeros(system.numberOfTimePoints,1);
for i=1:system.numberOfTimePoints
	noise.uNoise(i) = normrnd(0,noise.maxAmplitudeSM);
end
noise.noiseFFT = fft(noise.uNoise)';   
noise.noiseFFT_oneSided= noise.noiseFFT(1:calculation.numberOfFrequencies); % We now remove the folded part of the spectrum and multiply by 2 to keep the same energy
noise.noiseAbsFFT_oneSided = abs(noise.noiseFFT_oneSided); % Is this really an energy?
calculation.S1 = 1; %as we are using a rectangular windows
noise.noisePowerSpectrum = 2*noise.noiseAbsFFT_oneSided.^2/calculation.S1^2; % According to report GH_FFT from G. Heinzel
noise.noiseAmplitudeSpectrum = sqrt(noise.noisePowerSpectrum);
calculation.S2 = 1; %as we are using a rectangular windows
noise.noiseAmplitudeSpectralDensity = 2*noise.noiseAbsFFT_oneSided.^2/(calculation.S2^2*system.samplingFrequency);
noise.stdNoiseAmplitudeSpectrum = std(noise.noiseAmplitudeSpectrum);
noise.meanNoiseAmplitudeSpectrum = mean(noise.noiseAmplitudeSpectrum);
noise.maxNoiseAmplitudeSpectrum = max(noise.noiseAmplitudeSpectrum);

results.wantedFreq = [3 4 5 6]*system.frequencyDrive;
figure('Name','SM')
for i=1:size(results.wantedFreq,2)
    [~,freqIndex] = find(system.freq-results.wantedFreq(i)==0);

    subplot(1,4,i)
    results.SMAbsFFT_oneSided = real(results.SM1(:,freqIndex)); % take the absolut value
    results.SMPowerSpectrum = 2*results.SMAbsFFT_oneSided.^2/calculation.S1^2; % According to report GH_FFT from G. Heinzel
    results.SMAmplitudeSpectrum = sqrt(results.SMPowerSpectrum);
    res = reshape(results.SMAmplitudeSpectrum(:),[system.sizeX,system.sizeY]);
    imagesc(system.x,system.y,log(res./noise.stdNoiseAmplitudeSpectrum));
    axis square
    caxis([1 5])
    title(sprintf('log(SM SNR) @ %i kHz',system.freq(freqIndex)/1000))
    xlabel('x axis / m');
    ylabel('y axis / m');
end
%colorbar; %display the color bar on the last picture

figure('Name','Reco')
[~,i] = min(results.errorEstimate);
%i=system.maxIterationReco;
res = reshape(results.X(:,i),[system.sizeX,system.sizeY]);
imagesc(system.x,system.y,real(res));
colormap('gray')
axis square
hold all
title(sprintf('Reconstructed image at iteration %i - SNR > %i',i,system.SNRLimits));
xlabel('x axis / m');
ylabel('y axis / m');

clear('res','freqIndex','i')