clear all
close all

%% 1. Define the access path to the required other packages.
% This are typically the sphericalHarmonics and  ScannerDesign packages.
disp('1. Define the access path to the required other packages.')

addpath(genpath(fullfile('.')))
addpath(genpath(fullfile('..','SphericalHarmonics')))
addpath(genpath(fullfile('..','ScannerDesign')))

%% 2. Loading the scanner model.
% It is made of at least the spherical harmonics coefficient, the radius in which they are define and the nominal current used for the coil.
disp('2. Loading the scanner''s coil codel')
load('IdealFFL.mat');

% Scanner properties
system.c11 = Quadru_0.bc(1).coefficient(2,2)/Quadru_0.rhoReference*Quadru_0.current; % "Gradient" on the line
system.aDx = Drive_X.bc(1).coefficient(1,1)*Drive_X.current;
system.aDy = Drive_Y.bc(2).coefficient(1,1)*Drive_Y.current;
fprintf('c11 %g T/m; aDx %g T; aDy %g T;\n',system.c11,system.aDx,system.aDx);
fprintf('Drive peak translation x:%g m; y:%g m\n',system.aDx/system.c11,system.aDy/system.c11);

%% 3. Define all the other system parameters.
% The spatial resolution of the phantom and the system matrix. The used frequencies, time vector, sampling frequencies, noise parameters etc. This is highly dependent on the scanner you want to model.
disp('3. Define all the other system parameters.')

% We asumme that we have squarre voxels with size dx*dy*dz
% The resolution used for the SM and the phantom have to be different
% This is done to avoid the 'Inverse Crime'. See "Introduction into MPI"
% Knopp & Buzug 2011 p140 �5.4.4
calculation.dxSM = 1*10^-3/2; % [m] size of a voxel in the x direction
calculation.dySM = 1*10^-3/2; % [m] size of a voxel in the y direction
calculation.dzSM = 1*10^-3/2; % [m] size of a voxel in the z direction
calculation.dxPH = calculation.dxSM/1.3; % [m] size of a voxel in the x direction
calculation.dyPH = calculation.dySM/1.3; % [m] size of a voxel in the y direction
calculation.dzPH = calculation.dzSM; % [m] size of a voxel in the z direction
calculation.radiusFoV = 0.09;

% frequencies & time signal
system.frequencyLineRotation = 100;% The corresponding line rotation frequency
system.frequencyQuadrupol = 2*system.frequencyLineRotation;
system.frequencyDrive = 25000; % Hz

system.nbrPeriods = 2; % To adapte the resolution of the spectra
system.samplingFrequency = 1e6; %[Hz] Sampling Frequency of the acquisition Board.
%Maximal and minimal frequencies in the spectrum
system.fSMmin = 45000; % to remove the first harmonics and related distortion around it
system.fSMmax = 500000; % to keep the memory usage low

system.timeLength = system.nbrPeriods/(2*system.frequencyLineRotation);
system.numberOfTimePoints  = system.timeLength*system.samplingFrequency; % This should be an even integer.
if mod(system.numberOfTimePoints,2) ~= 0
  error('system.numberOfTimePoints should be even for this script')
end
calculation.time = linspace(0,system.timeLength,system.numberOfTimePoints+1); % +1 to ge the right numer of point after the nextstep
calculation.time = calculation.time(1:end-1);% Remove the last point as we start at zero, and thus form a period.

calculation.dt= 1/50*10^-6;% time difference used to calculate the derivatives
calculation.deltaf = 1/system.timeLength;% resolution of our spectrum

calculation.numberOfFrequencies = system.timeLength*system.samplingFrequency/2+1; % maximal number of frequencies
system.freq = system.samplingFrequency/2*linspace(0,1,calculation.numberOfFrequencies); % frequencies of the spectrum
system.nbrPointPerDrivePeriod = (system.numberOfTimePoints-1)/system.nbrPeriods; % number of acquired point for a period of the drive signal

calculation.mu0 = 4*pi*1e-7; % permeability of the air.

% Particle parameters
system.particleDiameter = 30*10^-9;
system.MsatSinglePart = 0.6/calculation.mu0;
system.Tempe = 310;
system.concentrationPartiMax = 0.03*0.5*1000/10;% [mol/m^3] 3% of "good" particle,1000 is for l to m^�,/10 for a dilution
system.volumeSample = calculation.dxSM*calculation.dySM*calculation.dzSM; % [m^3] volume of a voxel

% We use the noise model Weiznecker 2007
noise.multiplier = 10;
calculation.kB  = 1.380650424e-23;
noise.Rp = 185*10^-6; % [Ohm]  according to Weizenecker 2007 - A simulation study...
noise.T = 310; % [K]
noise.deltaF = system.samplingFrequency/2; % [Hz]
noise.maxAmplitude = noise.multiplier*sqrt(4*calculation.kB*noise.T*noise.deltaF*noise.Rp);
noise.ASD = noise.multiplier*sqrt(4*calculation.kB*noise.T*noise.Rp);
noise.maxAmplitudeSM = noise.multiplier*sqrt(4*calculation.kB*noise.T*noise.deltaF*noise.Rp)/30; % To simulate the fact that we can average the system matrix , we reduce by a factor 30 the noise in the SM

% Reconstruction related parameters
system.SNRLimits = 8; % choosen alsmost arbitrarly
system.maxIterationReco = 20; % choosen alsmost arbitrarly

% defining the geometry of the Field of View
system.xSM = -0.0090:calculation.dxSM:0.0090; %The resolution have to be limited in order to be able to reconstruct
system.ySM = -0.0100:calculation.dySM:0.0100;
system.zSM = 0:calculation.dzSM:0;
system.xPH = -0.0100:calculation.dxPH:0.0100;
system.yPH = -0.0100:calculation.dyPH:0.0100;
system.zPH = 0:calculation.dzPH:0;

%Size of the matrix
system.sizeXSM = size(system.xSM,2);
system.sizeYSM = size(system.ySM,2);
system.sizeZSM = size(system.zSM,2);
system.sizeXPH = size(system.xPH,2);
system.sizeYPH = size(system.yPH,2);
system.sizeZPH = size(system.zPH,2);


%% 4. Calculate the field for the system matrix.
% From the spherical harmonics coefficient, the magnetic flux density is calculated. The receive coils fields are also defined here.
disp('4. Calculate the field for the system matrix.')
tic
Selection_Z.B = RebuildField7(Selection_Z.bc,Selection_Z.bs,Selection_Z.rhoReference,system.xSM,system.ySM,system.zSM,'sch');
Quadru_0.B = RebuildField7(Quadru_0.bc,Quadru_0.bs,Quadru_0.rhoReference,system.xSM,system.ySM,system.zSM,'sch');
Quadru_45.B = RebuildField7(Quadru_45.bc,Quadru_45.bs,Quadru_45.rhoReference,system.xSM,system.ySM,system.zSM,'sch');
Drive_X.B = RebuildField7(Drive_X.bc,Drive_X.bs,Drive_X.rhoReference,system.xSM,system.ySM,system.zSM,'sch');
Drive_Y.B  = RebuildField7(Drive_Y.bc,Drive_Y.bs,Drive_Y.rhoReference,system.xSM,system.ySM,system.zSM,'sch');

% the sensibility is here define as the field generated by a 1A current
% this is why we have to save the bc and bs for a unit current
system.s1 = Drive_X.B;
calculation.s1x = system.s1(1,:);
calculation.s1y = system.s1(2,:);
calculation.s1z = system.s1(3,:);

system.s2 = Drive_Y.B;
calculation.s2x = system.s2(1,:);
calculation.s2y = system.s2(2,:);
calculation.s2z = system.s2(3,:);

clear('i','j')
fprintf('Time taken %2.0f s.\n', toc)
%% 5. Define the time-varying amplitude applied on the different coils. Also named sequence.
% As we are going to make a derivative, two similare vectors are created for each amplitude.
% One at t and one at t+dt.

disp('5. Define the time-varying amplitude applied on the different coils.')
c1 = cos(2*pi*system.frequencyQuadrupol*calculation.time(:).');
c2 = sin(2*pi*system.frequencyQuadrupol*calculation.time(:).');
c3 = sin(0.5*2*pi*system.frequencyQuadrupol*calculation.time(:).');
c4 = -cos(0.5*2*pi*system.frequencyQuadrupol*calculation.time(:).');
c5 = sin(2*pi*system.frequencyDrive*calculation.time(:)');

c1_dt = cos(2*pi*system.frequencyQuadrupol*(calculation.time(:).'+calculation.dt));
c2_dt = sin(2*pi*system.frequencyQuadrupol*(calculation.time(:).'+calculation.dt));
c3_dt = sin(0.5*2*pi*system.frequencyQuadrupol*(calculation.time(:).'+calculation.dt));
c4_dt = -cos(0.5*2*pi*system.frequencyQuadrupol*(calculation.time(:).'+calculation.dt));
c5_dt = sin(2*pi*system.frequencyDrive*(calculation.time(:).'+calculation.dt));

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
clear('c1','c2','c3','c4','c5','c1_dt','c2_dt','c3_dt','c4_dt','c5_dt')

%% 6. Create the time-varying magnetic flux-density for the SM
disp('6. Create the time-varying magnetic flux-density for the SM')
tic
Bx =system.coefSelection_Z.'*   Selection_Z.current*Selection_Z.B(1,:)+...
    system.coefQuadru_0.'   *	Quadru_0.current*Quadru_0.B(1,:)+...
    system.coefQuadru_45.'  *	Quadru_45.current*Quadru_45.B(1,:)+...
    system.coefDrive_X.'    *	Drive_X.current*Drive_X.B(1,:)+...
    system.coefDrive_Y.'    *	Drive_Y.current*Drive_Y.B(1,:);

By =system.coefSelection_Z.'*   Selection_Z.current*Selection_Z.B(2,:)+...
    system.coefQuadru_0.'   *	Quadru_0.current*Quadru_0.B(2,:)+...
    system.coefQuadru_45.'  *	Quadru_45.current*Quadru_45.B(2,:)+...
    system.coefDrive_X.'    *	Drive_X.current*Drive_X.B(2,:)+...
    system.coefDrive_Y.'    *	Drive_Y.current*Drive_Y.B(2,:);

Bz =system.coefSelection_Z.'*   Selection_Z.current*Selection_Z.B(3,:)+...
    system.coefQuadru_0.'   *	Quadru_0.current*Quadru_0.B(3,:)+...
    system.coefQuadru_45.'  *	Quadru_45.current*Quadru_45.B(3,:)+...
    system.coefDrive_X.'    *	Drive_X.current*Drive_X.B(3,:)+...
    system.coefDrive_Y.'    *	Drive_Y.current*Drive_Y.B(3,:);

Babs = sqrt(Bx.^2+By.^2+Bz.^2);

Bx_dt = system.coefSelection_Z_dt.'*   Selection_Z.current*Selection_Z.B(1,:)+...
        system.coefQuadru_0_dt.'   *	Quadru_0.current*Quadru_0.B(1,:)+...
        system.coefQuadru_45_dt.'  *	Quadru_45.current*Quadru_45.B(1,:)+...
        system.coefDrive_X_dt.'    *	Drive_X.current*Drive_X.B(1,:)+...
        system.coefDrive_Y_dt.'    *	Drive_Y.current*Drive_Y.B(1,:);

By_dt = system.coefSelection_Z_dt.'*   Selection_Z.current*Selection_Z.B(2,:)+...
        system.coefQuadru_0_dt.'   *	Quadru_0.current*Quadru_0.B(2,:)+...
        system.coefQuadru_45_dt.'  *	Quadru_45.current*Quadru_45.B(2,:)+...
        system.coefDrive_X_dt.'    *	Drive_X.current*Drive_X.B(2,:)+...
        system.coefDrive_Y_dt.'    *	Drive_Y.current*Drive_Y.B(2,:);

Bz_dt = system.coefSelection_Z_dt.'*   Selection_Z.current*Selection_Z.B(3,:)+...
        system.coefQuadru_0_dt.'   *	Quadru_0.current*Quadru_0.B(3,:)+...
        system.coefQuadru_45_dt.'  *	Quadru_45.current*Quadru_45.B(3,:)+...
        system.coefDrive_X_dt.'    *	Drive_X.current*Drive_X.B(3,:)+...
        system.coefDrive_Y_dt.'    *	Drive_Y.current*Drive_Y.B(3,:);

Babs_dt = sqrt(Bx_dt.^2+By_dt.^2+Bz_dt.^2);
    
fprintf('Time taken %2.0f s.\n', toc)
%% 7. The magnetic flux-density can be displayed in a 2D plane for all time t,
% to check your MPI-signal generating volume's shape
% 
% figure
% threshold = [1 2 3]*10^-3;
% for i=1:13:system.numberOfTimePoints/2
%     image = reshape(Babs(i,:),[system.sizeXSM,system.sizeYSM])';
%     for j=1:size(threshold,2)
%         [C,h] = contour(system.xSM,system.ySM,image,[threshold(j) threshold(j)]);
%         set(h,'LineWidth',2);
%         hold all;
%     end
%     hold off
%     %imagesc(system.xSM,system.ySM,image);
%     xlabel('x axis /m')
%     ylabel('y axis /m')
%     axis square
%     set(gca,'YDir','normal');
%     legend('1 mT','2 mT','3 mT')
%     %caxis([-threshold threshold]);
%     title(sprintf('|B| /T. Time %3.3f ms',calculation.time(i)*1000))
%     pause(1/100)
% end

%% 8. Calculate the time-varying magnetic moment for each point in space for the SM.
% Often using the Langevin model. 
% Note that the magnetic flux density matrix are used to stored the 
% obtained values, in order to save memory. If you have enough memory on 
% your system, you can of-course save the magnetic moment in another matrix.

disp('8. Calculate the time-varying magnetic moment for each point in space for the SM.')
tic
for i=1:system.numberOfTimePoints,
    [Bx(i,:),By(i,:),Bz(i,:)] = langevinParticle4(Bx(i,:),By(i,:),Bz(i,:),Babs(i,:),system.particleDiameter,system.MsatSinglePart,system.Tempe,system.concentrationPartiMax,system.volumeSample);
    [Bx_dt(i,:),By_dt(i,:),Bz_dt(i,:)] = langevinParticle4(Bx_dt(i,:),By_dt(i,:),Bz_dt(i,:),Babs_dt(i,:),system.particleDiameter,system.MsatSinglePart,system.Tempe,system.concentrationPartiMax,system.volumeSample);
end
mx = Bx;
my = By;
mz = Bz;
mx_dt = Bx_dt;
my_dt = By_dt;
mz_dt = Bz_dt;
clear('Bx','By','Bz','Babs','Bx_dt','By_dt','Bz_dt','Babs_dt')
fprintf('Time taken %2.0f s.\n', toc)

%% 9. Calculate the induce votlage for the SM and the SM.
% Here the SM is calculated, reduced to his one-sided form and the energy is corrected.
disp('9. Calculate the induced voltage for the SM and the SM.')
tic
u1=zeros(system.numberOfTimePoints,system.sizeXSM*system.sizeYSM*system.sizeZSM);
u2=zeros(system.numberOfTimePoints,system.sizeXSM*system.sizeYSM*system.sizeZSM);

for i=1:system.numberOfTimePoints
    u1(i,:) = ((calculation.s1x.*(mx_dt(i,:)-mx(i,:))/calculation.dt)+(calculation.s1y.*(my_dt(i,:)-my(i,:))/calculation.dt)+(calculation.s1z.*(mz_dt(i,:)-mz(i,:))/calculation.dt)) + normrnd(0,noise.maxAmplitudeSM);
    u2(i,:) = ((calculation.s2x.*(mx_dt(i,:)-mx(i,:))/calculation.dt)+(calculation.s2y.*(my_dt(i,:)-my(i,:))/calculation.dt)+(calculation.s2z.*(mz_dt(i,:)-mz(i,:))/calculation.dt)) + normrnd(0,noise.maxAmplitudeSM);
end

clear('results','mx','my','mz','mx_dt','my_dt','mz_dt')
results.SM1 = fft(u1).'/system.numberOfTimePoints; % we take the FFT
results.SM1 = 2*results.SM1(:,1:calculation.numberOfFrequencies);% and use the one sided part
results.SM1(:,1) = results.SM1(:,1)/2;% Correct the energy
results.SM1(:,end) = results.SM1(:,end)/2;% Correct the energy

results.SM2 = fft(u2).'/system.numberOfTimePoints; % we take the FFT
results.SM2 = 2*results.SM2(:,1:calculation.numberOfFrequencies);% and use the one sided part
results.SM2(:,1) = results.SM2(:,1)/2;% Correct the energy
results.SM2(:,end) = results.SM2(:,end)/2;% Correct the energy

%save('SM_FFP.mat','SM','-v7.3') %Watch out, the SM are often quite big.

clear('Imax','Imin','u1','u2','i','j');
fprintf('Time taken %2.0f s.\n', toc);


figure('Name','SM - absolute value of the first frequency components')
freqWithEnergie = [248,250,252,254,499,501,503,748,750,752,754];
maxSM = zeros(1,size(freqWithEnergie,2));
for i=1:size(freqWithEnergie,2)
    subplot(2,11,i)
    results.SM1_oneSided = 2*results.SM1(:,freqWithEnergie(i));
    results.SM1_oneSided(1) = results.SM1(1,freqWithEnergie(i));
    res = reshape(abs(results.SM1_oneSided(:)),[system.sizeXSM,system.sizeYSM]);
    imagesc(system.xSM,system.ySM,res);
    maxSM(i) = max(max(abs(results.SM1_oneSided(:))));
    axis square
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    title(sprintf('%i',freqWithEnergie(i)))
    colormap('gray')
end

for i=1:size(freqWithEnergie,2)
    subplot(2,11,i+11)
    results.SM2_oneSided = 2*results.SM2(:,freqWithEnergie(i));
    results.SM2_oneSided(1) = results.SM2(1,freqWithEnergie(i));
    res = reshape(abs(results.SM2_oneSided(:)),[system.sizeXSM,system.sizeYSM]);
    imagesc(system.xSM,system.ySM,res);
    maxSM(i) = max(max(abs(results.SM2_oneSided(:))));
    axis square
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    title(sprintf('%i',freqWithEnergie(i)))
    colormap('gray')
end


%% 10. Make the phantom
% we use the same particle as for the SM
disp('10. Make the phantom.')
phantom.shape = [system.sizeXPH system.sizeYPH system.sizeZPH];
phantom.particleDiameter = system.particleDiameter;
phantom.concentrationPartiMax =  system.concentrationPartiMax;
phantom.shapeScaled = createResolutionPhantomGael4(phantom.shape, 8);

%Make sure of the scaling of the phantom
phantom.shapeScaled = phantom.shapeScaled/max(phantom.shapeScaled(:));
phantom.volumeSample = calculation.dxPH*calculation.dyPH*calculation.dzPH; % [m^3] volume of a voxel
phantom.MsatSinglePart = system.MsatSinglePart;

figure;imagesc(system.xPH,system.yPH,phantom.shapeScaled);colormap('gray'); axis image

% apply the field calculation limitation of the spherical harmonics to the
% phantom shape
for i=1:system.sizeXPH
    for j=1:system.sizeYPH
        if sqrt(system.xPH(i)^2+system.yPH(j)^2)>calculation.radiusFoV
            phantom.shapeScaled(i,j) = 0;
        end
    end
end

%% 11. Calculate the fields for the phantom measurement.

disp('11. Calculate the fields for the phantom measurement.')
Selection_Z.B = RebuildField7(Selection_Z.bc,Selection_Z.bs,Selection_Z.rhoReference,system.xPH,system.yPH,system.zPH,'sch');
Quadru_0.B = RebuildField7(Quadru_0.bc,Quadru_0.bs,Quadru_0.rhoReference,system.xPH,system.yPH,system.zPH,'sch');
Quadru_45.B = RebuildField7(Quadru_45.bc,Quadru_45.bs,Quadru_45.rhoReference,system.xPH,system.yPH,system.zPH,'sch');
Drive_X.B = RebuildField7(Drive_X.bc,Drive_X.bs,Drive_X.rhoReference,system.xPH,system.yPH,system.zPH,'sch');
Drive_Y.B  = RebuildField7(Drive_Y.bc,Drive_Y.bs,Drive_Y.rhoReference,system.xPH,system.yPH,system.zPH,'sch');

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

%% 12. Create the time-varying magnetic flux-density for the phantom measurement
tic
disp('12. Create the time-varying magnetic flux-density for the phantom measurement.')
Bx =system.coefSelection_Z.'*   Selection_Z.current*Selection_Z.B(1,:)+...
    system.coefQuadru_0.'   *	Quadru_0.current*Quadru_0.B(1,:)+...
    system.coefQuadru_45.'  *	Quadru_45.current*Quadru_45.B(1,:)+...
    system.coefDrive_X.'    *	Drive_X.current*Drive_X.B(1,:)+...
    system.coefDrive_Y.'    *	Drive_Y.current*Drive_Y.B(1,:);

By =system.coefSelection_Z.'*   Selection_Z.current*Selection_Z.B(2,:)+...
    system.coefQuadru_0.'   *	Quadru_0.current*Quadru_0.B(2,:)+...
    system.coefQuadru_45.'  *	Quadru_45.current*Quadru_45.B(2,:)+...
    system.coefDrive_X.'    *	Drive_X.current*Drive_X.B(2,:)+...
    system.coefDrive_Y.'    *	Drive_Y.current*Drive_Y.B(2,:);

Bz =system.coefSelection_Z.'*   Selection_Z.current*Selection_Z.B(3,:)+...
    system.coefQuadru_0.'   *	Quadru_0.current*Quadru_0.B(3,:)+...
    system.coefQuadru_45.'  *	Quadru_45.current*Quadru_45.B(3,:)+...
    system.coefDrive_X.'    *	Drive_X.current*Drive_X.B(3,:)+...
    system.coefDrive_Y.'    *	Drive_Y.current*Drive_Y.B(3,:);

Babs = sqrt(Bx.^2+By.^2+Bz.^2);

Bx_dt = system.coefSelection_Z_dt.'*   Selection_Z.current*Selection_Z.B(1,:)+...
        system.coefQuadru_0_dt.'   *	Quadru_0.current*Quadru_0.B(1,:)+...
        system.coefQuadru_45_dt.'  *	Quadru_45.current*Quadru_45.B(1,:)+...
        system.coefDrive_X_dt.'    *	Drive_X.current*Drive_X.B(1,:)+...
        system.coefDrive_Y_dt.'    *	Drive_Y.current*Drive_Y.B(1,:);

By_dt = system.coefSelection_Z_dt.'*   Selection_Z.current*Selection_Z.B(2,:)+...
        system.coefQuadru_0_dt.'   *	Quadru_0.current*Quadru_0.B(2,:)+...
        system.coefQuadru_45_dt.'  *	Quadru_45.current*Quadru_45.B(2,:)+...
        system.coefDrive_X_dt.'    *	Drive_X.current*Drive_X.B(2,:)+...
        system.coefDrive_Y_dt.'    *	Drive_Y.current*Drive_Y.B(2,:);

Bz_dt = system.coefSelection_Z_dt.'*   Selection_Z.current*Selection_Z.B(3,:)+...
        system.coefQuadru_0_dt.'   *	Quadru_0.current*Quadru_0.B(3,:)+...
        system.coefQuadru_45_dt.'  *	Quadru_45.current*Quadru_45.B(3,:)+...
        system.coefDrive_X_dt.'    *	Drive_X.current*Drive_X.B(3,:)+...
        system.coefDrive_Y_dt.'    *	Drive_Y.current*Drive_Y.B(3,:);

Babs_dt = sqrt(Bx_dt.^2+By_dt.^2+Bz_dt.^2);
    
clear('c1','c2','c3','c4','c5','c1_dt','c2_dt','c3_dt','c4_dt','c5_dt')
fprintf('Time taken %2.0f s.\n', toc)


%% 13. Calculate the time-varying magnetic moment for the phantom measurment.
% Often using the Langevin model. 
% Note that the magnetic flux density matrix are used to stored the 
% obtained values, in order to save memory. If you have enough memory on 
% your system, you can of-course save the magnetic moment in another matrix.

disp('13. Calculate the time-varying magnetic moment for each point of the phantom.')
tic
for i=1:system.numberOfTimePoints,
    [Bx(i,:),By(i,:),Bz(i,:)] = langevinParticle4(Bx(i,:),By(i,:),Bz(i,:),Babs(i,:),phantom.particleDiameter,phantom.MsatSinglePart,system.Tempe,phantom.concentrationPartiMax,phantom.volumeSample);
    [Bx_dt(i,:),By_dt(i,:),Bz_dt(i,:)] = langevinParticle4(Bx_dt(i,:),By_dt(i,:),Bz_dt(i,:),Babs_dt(i,:),phantom.particleDiameter,phantom.MsatSinglePart,system.Tempe,phantom.concentrationPartiMax,phantom.volumeSample);
end
mx = Bx;
my = By;
mz = Bz;
mx_dt = Bx_dt;
my_dt = By_dt;
mz_dt = Bz_dt;
clear('Bx','By','Bz','Babs','Bx_dt','By_dt','Bz_dt','Babs_dt')
fprintf('Time taken %2.0f s.\n', toc)

%% 14. Multiplying the time-varying magnetic moment for each point of the 
% phantom with the phantom's tracer concentration.
disp('14. Multiplying the time-varying magnetic moment for each point of the phantom with the phantom''s tracer concentration.');
tic
for i=1:system.numberOfTimePoints,
    mx(i,:) = mx(i,:).*phantom.shapeScaled(:).';
    my(i,:) = my(i,:).*phantom.shapeScaled(:).';
    mz(i,:) = mz(i,:).*phantom.shapeScaled(:).';
    mx_dt(i,:) = mx_dt(i,:).*phantom.shapeScaled(:).';
    my_dt(i,:) = my_dt(i,:).*phantom.shapeScaled(:).';
    mz_dt(i,:) = mz_dt(i,:).*phantom.shapeScaled(:).';
end

clear('i');
fprintf('Time taken %2.0f s.\n', toc)

%% 15. Calculate the induced voltage for the phantom measurement.
% Here the FFT is calculated, reduced to his one-sided form and the energy is corrected.
disp('15. Calculate the induced voltage for the phantom measurement.')

tic
results.u1_2=zeros(system.numberOfTimePoints,1);
results.u2_2=zeros(system.numberOfTimePoints,1);
for i=1:system.numberOfTimePoints
    results.u1_2(i) = sum(((calculation.s1x.*(mx_dt(i,:)-mx(i,:))/calculation.dt)+(calculation.s1y.*(my_dt(i,:)-my(i,:))/calculation.dt)+(calculation.s1z.*(mz_dt(i,:)-mz(i,:))/calculation.dt)))+ normrnd(0,noise.maxAmplitude);
    results.u2_2(i) = sum(((calculation.s2x.*(mx_dt(i,:)-mx(i,:))/calculation.dt)+(calculation.s2y.*(my_dt(i,:)-my(i,:))/calculation.dt)+(calculation.s2z.*(mz_dt(i,:)-mz(i,:))/calculation.dt)))+ normrnd(0,noise.maxAmplitude);
end
clear('mx','my','mz','mx_dt','my_dt','mz_dt');

results.signal1FFT = fft(results.u1_2).'/system.numberOfTimePoints; %We store in the line each point, in the column the time point
results.signal1FFT_oneSided = 2*results.signal1FFT(1:calculation.numberOfFrequencies);% We now remove the folded part of the spectrum
results.signal1FFT_oneSided(:,1) = results.signal1FFT_oneSided(:,1)/2;% Correct the energy
results.signal1FFT_oneSided(:,end) = results.signal1FFT_oneSided(:,end)/2;% Correct the energy
results.signal1AbsFFT_oneSided = abs(results.signal1FFT_oneSided);
calculation.S1 = 1; %as we are using a rectangular windows
results.signal1PowerSpectrum = results.signal1AbsFFT_oneSided.^2/calculation.S1^2; % According to report GH_FFT from G. Heinzel
results.signal1AmplitudeSpectrum = sqrt(results.signal1PowerSpectrum);

results.signal2FFT = fft(results.u2_2).'/system.numberOfTimePoints;
results.signal2FFT_oneSided = 2*results.signal2FFT(1:calculation.numberOfFrequencies);
results.signal2FFT_oneSided(:,1) = results.signal2FFT_oneSided(:,1)/2;% Correct the energy
results.signal2FFT_oneSided(:,end) = results.signal2FFT_oneSided(:,end)/2;% Correct the energy
results.signal2AbsFFT_oneSided = abs(results.signal2FFT_oneSided);
results.signal2PowerSpectrum = results.signal2AbsFFT_oneSided.^2/calculation.S1^2; % According to report GH_FFT from G. Heinzel
results.signal2AmplitudeSpectrum = sqrt(results.signal2PowerSpectrum);

figure('Name','Signal')

subplot(5,1,1)
hold all
plot(calculation.time*1000,system.coefDrive_X*Drive_X.current)
plot(calculation.time*1000,system.coefDrive_Y*Drive_Y.current)
xlabel('Time / ms')
ylabel('Current amplitude / A')
title('Drive fields current')

subplot(5,1,2)
plot(calculation.time*1000,results.u1_2)
xlabel('Time / ms')
ylabel('Voltage amplitude / V')
title('Induced voltage in the x canal by the particles only.')

subplot(5,1,3)
plot(calculation.time*1000,results.u2_2)
xlabel('Time / ms')
ylabel('Voltage amplitude / V')
title('Induced voltage in the y canal by the particles only.')

%find the fDx
[v c] = min(abs(system.freq-system.frequencyDrive));
firstharmonicIndex = c;
distanceBetweenharmonics = c-1;
nbrHarmonicsRemoved = 0.5;
nbrcolum = size(results.signal1FFT,2);

results.signal1FilteredFFT  = fft(results.u1_2).'/system.numberOfTimePoints;
results.signal1FilteredFFT(1:(firstharmonicIndex+nbrHarmonicsRemoved*distanceBetweenharmonics))=0;
results.signal1FilteredFFT(nbrcolum-distanceBetweenharmonics-nbrHarmonicsRemoved*distanceBetweenharmonics+1:nbrcolum)=0;

results.signal2FilteredFFT  = fft(results.u2_2).'/system.numberOfTimePoints;
results.signal2FilteredFFT(1:(firstharmonicIndex+nbrHarmonicsRemoved*distanceBetweenharmonics))=0;
results.signal2FilteredFFT(nbrcolum-distanceBetweenharmonics-nbrHarmonicsRemoved*distanceBetweenharmonics+1:nbrcolum)=0;

subplot(5,1,4)
plot(calculation.time*1000,ifft(results.signal1FilteredFFT))
xlabel('Time / ms')
ylabel('Voltage /?')
title('Induced voltage in the x canal (signal filtered up to 1.5 f_D)')


subplot(5,1,5)
plot(calculation.time*1000,ifft(results.signal2FilteredFFT))
xlabel('Time / ms')
ylabel('Voltage /?')
title('Induced voltage in the x canal (signal filtered up to 1.5 f_D)')

%% 16. Calculate the SNR.
disp('16. Calculate the SNR.')
noise.uNoise = zeros(system.numberOfTimePoints,1);
for i=1:system.numberOfTimePoints
	noise.uNoise(i) = normrnd(0,noise.maxAmplitude);
end
noise.noiseFFT = fft(noise.uNoise).'/system.numberOfTimePoints;   
noise.noiseFFT_oneSided= 2*noise.noiseFFT(1:calculation.numberOfFrequencies); % We now remove the folded part of the spectrum and multiply by 2 to keep the same energy
noise.noiseFFT_oneSided(:,1) = noise.noiseFFT_oneSided(:,1)/2;% Correct the energy
noise.noiseFFT_oneSided(:,end) = noise.noiseFFT_oneSided(:,end)/2;% Correct the energy
noise.noiseAbsFFT_oneSided = abs(noise.noiseFFT_oneSided);
calculation.S1 = 1; %as we are using a rectangular windows
noise.noisePowerSpectrum = noise.noiseAbsFFT_oneSided.^2/calculation.S1^2;
noise.noiseAmplitudeSpectrum = sqrt(noise.noisePowerSpectrum);
noise.stdNoiseAmplitudeSpectrum = std(noise.noiseAmplitudeSpectrum);

% Calculating the SNR of the phantom signal
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

%% 17. Truncating signals and SM to remove as much unnecessary information as possible
disp('17. Truncating signals and SM.')
index = 1;
results.tSNR1 = zeros(1,results.nbrGoodFrequency1-1);
results.tFreq1 = zeros(1,results.nbrGoodFrequency1-1);
results.tSM1 = zeros(size(results.SM1,1),results.nbrGoodFrequency1-1);
results.tSignal1FFT_oneSided = zeros(1,results.nbrGoodFrequency1-1);
for i=2:size(results.signal1FFT_oneSided,2)
    if results.signal1SNR(i) > system.SNRLimits && system.freq(i) > system.fSMmin && system.freq(i) < system.fSMmax
        results.tSNR1(index) = results.signal1SNR(i); % used in a figure
        results.tSignal1FFT_oneSided(index) = results.signal1FFT_oneSided(i); % used in reco
        results.tSM1(:,index) = results.SM1(:,i); % used in reco
        results.tFreq1(index) = system.freq(i); % used in a figure
        index = index+1;
    end
end

% second signal
index = 1;
results.tSNR2 = zeros(1,results.nbrGoodFrequency2-1);
results.tFreq2 = zeros(1,results.nbrGoodFrequency2-1);
results.tSM2 = zeros(size(results.SM2,1),results.nbrGoodFrequency2-1);
results.tSignal2FFT_oneSided = zeros(1,results.nbrGoodFrequency2-1);
for i=2:calculation.numberOfFrequencies
    if results.signal2SNR(i) > system.SNRLimits && system.freq(i) > system.fSMmin && system.freq(i) < system.fSMmax
        results.tSNR2(index) = results.signal2SNR(i); % used in a figure
        results.tSignal2FFT_oneSided(index) = results.signal2FFT_oneSided(i); % used in reco
        results.tSM2(:,index) = results.SM2(:,i); % used in reco
        results.tFreq2(index) = system.freq(i); % used in a figure
        index = index+1;
    end
end

clear('index','i','j');
fprintf('Time taken %2.0f s.\n', toc)


figure('Name','Signal')
subplot(2,1,1)
stem(system.freq/system.frequencyDrive,results.signal1SNR,'Marker','None')
set(gca,'yscal','log')
hold all
stem(results.tFreq1/system.frequencyDrive,results.tSNR1,'Marker','None')
set(gca,'yscal','log')
xlabel('# Harmonic')
ylabel('SNR')
title('SNR x channel (based on the amplitude spectrum)')
xlim([0 10])
ylim([1 10^3]);

subplot(2,1,2)
stem(system.freq/system.frequencyDrive,results.signal2SNR,'Marker','None')
set(gca,'yscal','log')
hold all
stem(results.tFreq2/system.frequencyDrive,results.tSNR2,'Marker','None')
set(gca,'yscal','log')
xlabel('# Harmonic')
ylabel('SNR')
title('SNR y channel (based on the amplitude spectrum)')
xlim([0 10])
ylim([1 10^3]);

%% 18. Assemble the channel and reconstruct the truncated signals. (SNR thresholding)

disp('18. Assemble the channel and reconstruct the truncated signals.')
tic

S = [results.tSM1,results.tSM2].';
u = [results.tSignal1FFT_oneSided,results.tSignal2FFT_oneSided].';
[results.C,~,~] = artGael(S,u,system.maxIterationReco);

% least square solution
% lambda0 = trace(S'*S)/size(S,2);
% lambdaRela = 0.1;
% lambda = lambdaRela*lambda0;
% 
% [results.C,~,~] = artGael(S'*S+lambda*eye(size(S,2)),S'*u,system.maxIterationReco);

figure
for i=1:system.maxIterationReco
    subplot(1,2,1)
    imagesc(phantom.shapeScaled)
    axis square
    subplot(1,2,2)
    res = reshape(results.C(:,i),[system.sizeXSM,system.sizeYSM]);
    imagesc(system.xSM,system.ySM,real(res));
    colormap('gray')
    axis square
    title(sprintf('i=%i',i));
    pause(1/25)
end

clear('res')
fprintf('Time taken %2.0f s.\n', toc)

%% Figure
disp('display the results')
figure('Name','Signal')

subplot(3,1,1)
hold all
plot(calculation.time*1000,system.coefDrive_X*Drive_X.current)
plot(calculation.time*1000,system.coefDrive_Y*Drive_Y.current)
xlabel('Time / ms')
ylabel('Current amplitude / A')
title('Drive fields current')

subplot(3,1,2)
plot(calculation.time*1000,results.u1_2)
xlabel('Time / ms')
ylabel('Voltage amplitude / V')
title('Induced voltage in the x canal by the particles only.')

subplot(3,1,3)
stem(system.freq/system.frequencyDrive,results.signal1SNR,'Marker','None')
set(gca,'yscal','log')
hold all
stem(results.tFreq1/system.frequencyDrive,results.tSNR1,'Marker','None')
set(gca,'yscal','log')
xlabel('# Harmonic')
ylabel('SNR')
title('SNR (based on the amplitude spectrum)')
xlim([0 10])
ylim([1 10^3]);


figure('Name','SM 1 - absolute value of the first frequency components')
freqWithEnergie = [248,250,252,254,499,501,503,748,750,752,754];
maxSM = zeros(1,size(freqWithEnergie,2));
for i=1:size(freqWithEnergie,2)
    subplot(2,11,i)
    results.SM1_oneSided = 2*results.SM1(:,freqWithEnergie(i));
    results.SM1_oneSided(1) = results.SM1(1,freqWithEnergie(i));
    res = reshape(abs(results.SM1_oneSided(:)),[system.sizeXSM,system.sizeYSM]);
    imagesc(system.xSM,system.ySM,res);
    maxSM(i) = max(max(abs(results.SM1_oneSided(:))));
    axis square
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    title(sprintf('%i',freqWithEnergie(i)))
    colormap('gray')
end
subplot(2,11,12:22)
stem(results.signal1SNR,'Marker','None')
set(gca,'yscal','log')
xlabel('Frequency Components')
ylabel('SNR')
title('SNR (based on the amplitude spectrum)')
xlim([0 10])
ylim([1 10^3]);
xlim([1 1000]);

figure('Name','SM 2 - absolute value of the first frequency components')
freqWithEnergie = [248,250,252,254,499,501,503,748,750,752,754];
for i=1:size(freqWithEnergie,2)
    subplot(2,11,i)
    results.SM2_oneSided = 2*results.SM2(:,freqWithEnergie(i));
    results.SM2_oneSided(1) = results.SM2(1,freqWithEnergie(i));
    res = reshape(abs(results.SM2_oneSided(:)),[system.sizeXSM,system.sizeYSM]);
    imagesc(system.xSM,system.ySM,res);
    maxSM(i) = max(max(abs(results.SM2_oneSided(:))));
    axis square
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    title(sprintf('%i',freqWithEnergie(i)))
    colormap('gray')
end
subplot(2,11,12:22)
stem(results.signal2SNR,'Marker','None')
set(gca,'yscal','log')
xlabel('Frequency Components')
ylabel('SNR')
title('SNR (based on the amplitude spectrum)')
xlim([0 10])
ylim([1 10^3]);
xlim([1 1000]);

clear('freqIndex','maxSM','i','res')
