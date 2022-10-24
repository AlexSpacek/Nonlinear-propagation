function y = Dawson2(x)

% 2-nd order rational approximation of the Dawson integral function
% See www.dx.doi.org/10.3247/SL4Soft12.001

y = x*x;
p = 1.0 + y*(0.1329766230 + y*(0.0996005943));
q = 1.0 + y*(0.8544964660 + y*(0.2259838671 + 2*0.0996005943*y));
y = x*(p/q);
