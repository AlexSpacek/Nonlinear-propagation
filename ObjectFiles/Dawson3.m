function y = Dawson3(x)

% 3rd order rational approximation of the Dawson integral function
% See www.dx.doi.org/10.3247/SL4Soft12.001

y = x*x;
p = 1.0 + y*(0.1349423927 + y*(0.0352304655 ...
        + y*(0.0138159073)));
q = 1.0 + y*(0.8001569104 + y*(0.3190611301 ...
        + y*(0.0540828748 + 2*0.0138159073*y)));
y = x*(p/q);
