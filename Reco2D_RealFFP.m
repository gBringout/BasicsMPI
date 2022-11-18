%% 1. Loading the required external functions
% Based on https://github.com/MagneticParticleImaging/MDF/blob/master/matlab/reco.m
% Latest commit 632cf42
% The script is using the Reco algo from the MDF repo

addpath(genpath(fullfile('.')))
addpath(genpath(fullfile('..','MDF')))
%addpath(genpath(fullfile('..','ScannerDesign')))

clear all
close all


%% 2. Download measurement and systemMatrix from http://media.tuhh.de/ibi/mdf/

filenameSM = 'systemMatrix.mdf';
filenameMeas = 'measurement.mdf';

if exist('systemMatrix.mdf')==0
    websave(filenameSM,'http://media.tuhh.de/ibi/mdfv2/systemMatrix_V2.mdf')
end
if exist('measurement.mdf')==0
    websave(filenameMeas,'http://media.tuhh.de/ibi/mdfv2/measurement_V2.mdf')
end

%% 3. Define all the other system parameters.
analogFilterMinFreq = 80e3;


%% 4. Calculate the field for the system matrix.
%% 5. Define the time-varying amplitude applied on the different coils. Also named sequence.
%% 6. Create the time-varying magnetic flux-density for the SM
%% 7. The magnetic flux-density can be displayed in a 2D plane for all time t,
% to check your MPI-signal generating volume's shape
%% 8. Calculate the time-varying magnetization for each point in space for the SM.
%% 9. Loading the SM data
% For the System matrix (later named SM)
% to obtain infos on the file, use the command: infoSM = h5info(filename_SM);
% or read the format documentation

% read the data, saved as real numbers
S = h5read(filenameSM, '/measurement/data');

% reinterpret as complex numbers
S = complex(S.r,S.i);
% get rid of background frames
isBG = h5read(filenameSM, '/measurement/isBackgroundFrame');
S = S(isBG == 0,:,:,:);
%% 10. Make the phantom
%% 11. Calculate the fields for the phantom measurement.
%% 12. Create the time-varying magnetic flux-density for the phantom measurement
%% 13. Calculate the time-varying magnetization for the phantom measurment.
%% 14. Multiplying the time-varying magnetization for each point of the 
% phantom with the phantom's tracer concentration.
%% 15. Loading the Measurements data
% For the measurements
% read and convert the data as complex numbers
% note that these data contain 500 measurements
u = h5read(filenameMeas, '/measurement/data');
%u = squeeze(u(1,:,:,:) + 1i*u(2,:,:,:));
u = fft(cast(u,'double'));
u = u(1:(size(u,1)/2+1),:,:,:);
%% 16. Calculate the SNR.
%% 17. Pre-process - Remove the frequencies which are lower than 80 kHz, as they are unreliable due to the anologue filter in the scanner

% generate frequency vector
numFreq = h5read(filenameMeas, '/acquisition/receiver/numSamplingPoints')/2+1;
rxBandwidth = h5read(filenameMeas, '/acquisition/receiver/bandwidth');
freq = linspace(0,1,numFreq) .* rxBandwidth;

% we supose that the same frequencies are measured on all channel for 
% the SM and the measurements. use only x/y receive channels
idxFreq = freq > analogFilterMinFreq;
numFreq_truncated = sum(idxFreq);
freq_truncated = freq(idxFreq);
S_truncated = S(:,idxFreq,1:2);
u_truncated = u(idxFreq,1:2,:);

figure
subplot(2,2,1)
semilogy(freq,abs(real(u(:,1,:,42))),'*')
hold all
semilogy(freq,abs(real(u(:,2,:,42))),'*')
semilogy(freq,abs(real(u(:,3,:,42))),'*')
legend("1","2","3")
title("FFT of the acquired signal u")
ylabel("|Re(u)|")
xlabel("Frequency /Hz")

subplot(2,2,2)
semilogy(freq,abs(imag(u(:,1,:,42))),'*')
hold all
semilogy(freq,abs(imag(u(:,2,:,42))),'*')
semilogy(freq,abs(imag(u(:,3,:,42))),'*')
legend("1","2","3")
title("FFT of the acquired signal u")
ylabel("|Im(u)|")
xlabel("Frequency /Hz")

subplot(2,2,3)
semilogy(freq_truncated,abs(real(u_truncated(:,1,42))),'*')
hold all
semilogy(freq_truncated,abs(real(u_truncated(:,2,42))),'*')
legend("1","2")
title("FFT of the truncated acquired signal u")
ylabel("|Re(u)|")
xlabel("Frequency /Hz")

subplot(2,2,4)
semilogy(freq_truncated,abs(imag(u_truncated(:,1,42))),'*')
hold all
semilogy(freq_truncated,abs(imag(u_truncated(:,2,42))),'*')
legend("1","2")
title("FFT of the truncated acquired signal u")
ylabel("|Im(u)|")
xlabel("Frequency /Hz")


%% 18.1. Merge frequency and receive channel dimensions
S_truncated = reshape(S_truncated, size(S_truncated,1), size(S_truncated,2)*size(S_truncated,3));
u_truncated = reshape(u_truncated, size(u_truncated,1)*size(u_truncated,2), size(u_truncated,3));

%% New Step. Averaged the measurement used for the reconstruction over all temporal frames
u_mean_truncated = mean(u_truncated,2);

%% 18.2 Make two simple reconstructions
% a normalized regularized kaczmarz approach
%c_normReguArt = artGael(S_truncated',...
%                        u_mean_truncated,...
%                        10);

c_normReguArt = kaczmarzReg(S_truncated(:,:),...
                        u_mean_truncated(:),...
                        1,1*10^-6,0,1,1);
% Display an image
% read the original size of an image
number_Position = h5read(filenameSM, '/calibration/size');

figure
imagesc(real(reshape(c_normReguArt(:,end),number_Position(1),number_Position(2))));
colormap(gray); axis square
%title({'Regularized and modified ART - 3 channels';'1 iterations / lambda = 10^{-6} / real part'})


%% Figure
% disp('display the results')

figure('Name','SM 1 - real part of the first frequency components')
maxSM = zeros(1,100);
for i=2:100
    subplot(10,10,i)
    %units have to be checked
    SM_oneSided = S_truncated(:,i);
    %SM_oneSided(1) = results.SM1(1,i);
    res = reshape(real(SM_oneSided(:)),[number_Position(1),number_Position(2)]);
    imagesc(res);
    maxSM(i) = max(max(abs(SM_oneSided(:))));
    axis square
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    title(sprintf('%i',i-1))
    colormap('gray')
end
subplot(10,10,1)
plot(maxSM/max(maxSM))