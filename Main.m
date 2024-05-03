addpath(pwd)
addpath('ObjectFiles') 
clear classes
clear all
InitialParameters.Size=[1000,2^11]; %space,time
InitialParameters.CentralWavelength=820*1e-9;
InitialParameters.MinWavelength=480*1e-9; %min wavelength present in lambda vector
InitialParameters.FWHM=15e-15; 
InitialParameters.Width=5*1e-3;
InitialParameters.MaxSpace=4*InitialParameters.Width; %span of the space vector
InitialParameters.Energy=10*1e-3;
InitialParameters.DispersionCoefficients=[1200*1e-30,-2400*1e-45,0]; %GDD,TOD,FOD
field{1}=Field(InitialParameters);

%medium1:
length = 0.01;
splitsteps = 5;
runge_kutta_steps = 4;
n2 = 3e-20;
fusedSilica1=Material('FusedSilica',length,splitsteps,runge_kutta_steps,n2,field{1}.WavelengthVector,field{1}.CentralWavelength);

%lens1:
length = 0.001;
splitsteps = 2;
runge_kutta_steps = 2;
n2 = 0;
f = 1;
lens1=Lens('FusedSilica',length,splitsteps,runge_kutta_steps,n2,field{1}.WavelengthVector,field{1}.CentralWavelength,f);

%medium2:
length = 1;
splitsteps = 15;
runge_kutta_steps = 2;
n2 = 0;
vacuum1=Material('Vacuum',length,splitsteps,runge_kutta_steps,n2,field{1}.WavelengthVector,field{1}.CentralWavelength);

%compression1
GDD = -1200*1e-30;
TOD = 2400*1e-45;
FOD = 0;
compression1=Dispersion('compression1',GDD,TOD,FOD);


OpticalPath1={fusedSilica1,lens1,compression1,vacuum1};
figure1=figure;
figure1.Position = [100 100 800 800];
propagate(field{1},OpticalPath1,figure1);

