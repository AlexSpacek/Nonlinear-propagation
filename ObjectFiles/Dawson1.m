function y = Dawson1(x)

% 1st order rational approximation of the Dawson integral function
% See www.dx.doi.org/10.3247/SL4Soft12.001

y = x*x;
p = 1.0 + y*(0.4582332073);
q = 1.0 + y*(0.8041350741 + 2*0.4582332073*y);
y = x*(p/q);
