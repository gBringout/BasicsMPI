function [mx,my,mz,a] = langevinParticle4(Bx,By,Bz,Babs,D,Msat,T,ironConcentrationSample,vSample)
% Calculate the magnetic moment as used in simple simulation of MPI
% The particle magnetization is calculated using the langevin function 


% We have to provide
% Bx : [T] magnetic flux density in the x direction
% By : [T] magnetic flux density in the y direction
% Bz : [T] magnetic flux density in the z direction
% Babs : [T] absolute value of the magnetic flux density
% D : [m] particle core diameter 
% Msat : [A/m] [A.m^2/m^3] saturation magnetization (~0.6/mu0 = 4.7e+5)
% T : [K] absolute temperature (~310)
% ironConcentrationSample  : [mol(Fe)/m^3] concentration of usefull iron in Resovist
% vSample  : [m^3] volume of material in the sample

% output
% mx : [Am^2] resulting magnetic moment of the sample in x direction
% my : [Am^2] resulting magnetic moment of the sample in y direction 
% mz : [Am^2] resulting magnetic moment of the sample in z direction
% a : energy given to the particle

%% Calcuation of the needed values
kB  = 1.380650424e-23;      % [J/K] Boltzmann constant 
mu0 = 4*pi*1e-7;            % [N/A^2] permeability of free space  

Vpart  = 4/3*pi*(D/2)^3;    % [m^3] single particle volume     
mPartSat = Vpart * Msat;             % [A.m^2] magnetic moment at saturation for a single particle


% The idea is here to scale from the IRON concentration to the real amount
% of particle with the given radius.
% Molar mass of Fe3O4 = 3*55.845 + 4*15.999 = 231.531  g/mol
molarMassMagnetite = 231.531/1000; % [Kg/mol(Fe3O4)] it's equal to 231.531*10^-3 kg/mol (checked, it's ok)
densityMagnetite = 5.17*1000; % [Kg/m^3(Fe3O4)] see google "density Iron(II,III) oxide"
concentrationMagnetite = ironConcentrationSample/3;% [mol(Fe3O4)/m^3] Concentration of Magnetite / particle
massConcentrationMagnetite = concentrationMagnetite*molarMassMagnetite; % [Kg(Fe3O4)/m^3] mass of magnetite per volume of the sample
weightPart = Vpart*densityMagnetite; % [Kg] weight of a single particle
weigthMagnetiteSample = vSample*massConcentrationMagnetite; % [Kg] weigth of magnetit in our sample
nPart = weigthMagnetiteSample/weightPart; %[] number of particle in the sample
mSat = mPartSat*nPart; % [A.m^2] magnetic moment at saturation of our sample

%% Change the matrixes as vector
shapeFields = size(Bx);
vBx = Bx(:);
vBy = By(:);
vBz = Bz(:);
normB = Babs(:);   % norm of the field [T]

%% Computation of the Langevin function
a = mPart * normB/ (kB * T);

% we have to approximate the field if x is small (due to the coth function
if abs(a)<10^-9
% see http://ocw.nctu.edu.tw/upload/classbfs1209042706186759.pdf page 16
% or make a series expansion of the function yourself :)
	L = a/3; % Approximation
else
    L = coth(a) - 1./a; % Langevin function
end

%% Project back the magnetization

% Here we may have a division per zero
% So when normB is zero, the magnetization should be zero
% find where we have zero
indexZero = find(normB==0);

% Then we calculate everythings normally
% resulting magnetic moment  [Am^2]
vmx = mSat*L.*vBx./normB;            
vmy = mSat*L.*vBy./normB;
vmz = mSat*L.*vBz./normB;

% and we replace the value per zero where it should be zero
vmx(indexZero) = 0;
vmy(indexZero) = 0;
vmz(indexZero) = 0;

% reshape the vector as matrix
mx = reshape(vmx,shapeFields);
my = reshape(vmy,shapeFields);
mz = reshape(vmz,shapeFields);
end
