function [Mx,My,Mz,a] = langevinParticle4(Bx,By,Bz,Babs,D,Ms,T,ironConcentrationSample,vSample)
% Implement the langevin function used to simulate the MPI signal

% We have to provide
% Bx : [T] magnetic flux density in the x direction
% By : [T] magnetic flux density in the y direction
% Bz : [T] magnetic flux density in the z direction
% Babs : [T] absolute value of the magnetic flux density
% D : [m] particle core diameter 
% Ms : [A/m] [A.m^2/m^3] saturation magnetization (~0.6/mu0 = 4.7e+5)
% T : [K] absolute temperature (~310)
% ironConcentrationSample  : [mol(Fe)/m^3] concentration of usefull iron in Resovist
% vSample  : [m^3] volume of material in the sample

% output
% Mx : [Am^2/m^3] resulting particle magnetic moment in x direction per m^3 
% My : [Am^2/m^3] resulting particle magnetic moment in y direction per m^3 
% Mz : [Am^2/m^3] resulting particle magnetic moment in z direction per m^3
% a : energy given to the particle

%% Calcuation of the needed values
kB  = 1.380650424e-23;      % [J/K] Boltzmann constant 
mu0 = 4*pi*1e-7;            % [N/A^2] permeability of free space  

Vpart  = 4/3*pi*(D/2)^3;    % [m^3] single particle volume     
m = Vpart * Ms;             % [Am^2] magnetic moment at saturation for a single particle


% The idea is here to scale from the IRON concentration to the real amount
% of particle with the given radius.
% Molar mass of Fe3O4 = 3*55.845 + 4*15.999 = 231.531  g/mol
molarMassMagnetite = 231.531/1000; % [Kg/mol] it's equal to 231.531*10^-3 kg/mol (checked, it's ok)
densityMagnetite = 5.17*1000; % [Kg/m^3] see google "density Iron(II,III) oxide"
concentrationMagnetite = ironConcentrationSample/3;% [mol(Fe3O4)/m^3] Concentration of Magnetite / particle
massConcentrationMagnetite = concentrationMagnetite*molarMassMagnetite; % [Kg/m^3] mass of magnetite
M0 = m*vSample*massConcentrationMagnetite/(Vpart*densityMagnetite); % [Kg/m^3] mass of magnetite in our sample

%% Change the matrixes as vector
shapeFields = size(Bx);
vBx = Bx(:);
vBy = By(:);
vBz = Bz(:);
normB = Babs(:);   % norm of the field [T]


%% Computation of the Lnagevin function

a = m * normB/ (kB * T);

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
vMx = M0*L.*vBx./normB;            
vMy = M0*L.*vBy./normB;
vMz = M0*L.*vBz./normB;

% and we replace the value per zero where it should be zero
vMx(indexZero) = 0;
vMy(indexZero) = 0;
vMz(indexZero) = 0;

% reshape the vector as matrix
Mx = reshape(vMx,shapeFields);
My = reshape(vMy,shapeFields);
Mz = reshape(vMz,shapeFields);
end