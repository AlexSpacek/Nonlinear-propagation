function W = KeldyshFullRate(w,me,Eg,n,I)
% Function to calculate the full Keldysh rate (eq 37 from Keldysh ...
% (1965))
% Inputs:
% ? w: radial frequency of light (omega) [rad/s]
% ? me: effective electron mass [kg]
% ? : bandgap of material [J]
 % ? n: refractive index [unitless]
 % ? I: Laser Irradiance [W/mˆ2]
 % Outputs:
 % ? W: Keldysh photoionization rate [electrons/s/mˆ3]
 %

 % Note: There are two ways to calculate the dawson integral within...
% this function. One is through mfun, which is supplied with ...
% the symbolic math toolbox. The other is through dawson.m, ...
% which is a file from Matlab File Exchange that is faster.
 %
 % Authors: Troy Anderson, Chris Ferris
 % Last modified on 4/14/2014

 %% Constants:
 % Speed of light [m/s]
 c = 3e8;
 % Electron Charge [C]
 e = 1.6e-19;
 % Permittivity of Free Space [F/m]
 ep0 = 8.85e-12;
 % Planck Constant [J s]
 hbar = 1.054*10^(-34);

 %% Calculations
 % Electric Field Strength
 F = sqrt((2*I)./(c*n*ep0)); %[V/m]
 % Gamma, Keldysh Parameter
 gamma = (w./(e.*F)).*sqrt(me*Eg); % [unitless]
 % Create variables for common terms
 gg = gamma.^2./(1+gamma.^2);
 g1 = 1./(1+gamma.^2);

 % Elliptic Integrals
 % Elliptic Integrals: Keldysh expressions assume modulus k (m = k...
% ˆ2). Since ellipke uses modulus m, the expressions for gg and ...
% g1 are squared relative to those found in Keldysh
tic 
 [Kg,Egg] = ellipke(gg,eps*1e6);
 [K1,E1] = ellipke(g1,eps*1e6);

 Egtau = 2*Eg*sqrt(1+gamma.^2).*E1./(pi()*gamma);
 X = floor(Egtau./(hbar*w)+1);

 Wf1 = 2*w/(9*pi()) .* (sqrt(1+gamma.^2) * me * w ./ (gamma * ...
hbar)).^(3/2);
 Wf2 = Qfun(gamma,Egtau./(hbar*w));
 Wf3 = exp(-pi().*X.*(Kg-Egg)./E1);

 W = Wf1 .* Wf2 .* Wf3; % [electrons/s/mˆ3]

 % Set all NaN values to 0. NaNs can occur if the value of the ...
% intensity is too small. In this case, the photoionization rate...
% is negligible.
 W(isnan(W)) = 0;

 %% Nested Subfunctions
 % Q function (from Keldysh)

 function Q = Qfun(gamma,x)
 Q1 = sqrt(pi()./(2.*K1));
 Q2 = zeros(1,length(gamma));

 for i = 1:length(gamma)
 j = 0;
 tol = 1e-3;
 err = 1;
 OldQ2 = 0;
 while err > tol


 % Check to see if user?supplied 'dawson.m' exists ...
% in the path. This function is a faster ...
% implementation of the dawson integral than the ...
% 'mfun' implementation.
%  if exist('Dawson5.m','file') == 2
 Q2(i) = Q2(i) + exp(-pi() .* (Kg(i)-Egg(i)) .* ...
j ./ E1(i)) .* Dawson5(sqrt(pi()^2.*(2*floor...
(x(i)+1)-2.*x(i) + j) ./(2*K1(i) .* E1(i)))...
);
%  else
%  Q2(i) = Q2(i) + exp(?pi() .* (Kg(i)?Eg(i)) .* ...
% j ./ E1(i)) .* mfun('dawson',sqrt(pi()ˆ2....
% *(2*floor(x(i)+1)?2.*x(i) + j) ./(2*K1(i) ....
% * E1(i))));
%  end
 err = abs(Q2(i)-OldQ2);
 j = j + 1;
 OldQ2 = Q2(i);
 end
 end
 Q = Q1.*Q2;
 end
 end
