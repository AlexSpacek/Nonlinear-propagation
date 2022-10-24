function [W] = keldyshMPI(w,me,Eg,n,I)
%code taken from
% Function to calculate the Keldysh tunneling rate (eq 41 from ...
% Keldysh (1965))
%
% Inputs:
% ? w: radial frequency of light (omega) [rad/s]
% ? me: effective electron mass [kg]
% ? ?: bandgap of material [J]
% ? n: refractive index [unitless]
% ? I: Laser Irradiance [W/mˆ2]
% Outputs:
% ? W: Keldysh MPI Rate [electrons/s/mˆ3]
%
% Note: There are two ways to calculate the dawson integral within...
% this function. One is through mfun, which is supplied with ...
% the symbolic math toolbox. The other is through dawson.m, ...
% which is a file from Matlab File Exchange that is faster.
%

% Authors: Chris Ferris, Troy Anderson
% Last modified on 4/14/2014

%% Constants:
% Electron Charge [C]
e = 1.6e-19;
% Planck Constant [J s]
hbar = 1.054e-34;
% Speed of light (m/s)
c = 3e8;
% Permittivity of Free Space (F/m)
ep0 = 8.85e-12;

 %% Calculations
 % Electric Field Strength [V/m]
 F = sqrt((2*I)./(c*n*ep0));

Eg_bar = Eg+(e^2.*F.^2)./(4*me*w^2);

X = fix(Eg_bar./(hbar.*w) +1);

Wmpi1 = (2*w)/(9*pi()) * ((me*w)/hbar)^(3/2);

% Check to see if user?supplied 'dawson.m' exists in the path. ...
% This function is a faster implementation of the dawson integral...
% than the 'mfun' implementation.
if exist('dawson.m','file') == 2
Wmpi2 = dawson(((2.*X-(2.*Eg_bar)./(hbar*w)).^(1/2)));
else
%     tic 
Wmpi2 = mfun('dawson',((2.*X-(2.*Eg_bar)./(hbar*w)).^(1/2)));
% toc
end


Wmpi3 = exp(2.*X.*(1-(e^2.*F.^2)./(4*me*w^2*Eg)));
Wmpi4 = ((e^2.*F.^2)./(16*me*w^2*Eg)).^X;

W = Wmpi1 .* Wmpi2 .* Wmpi3 .* Wmpi4; % [electrons/s/mˆ3]

end