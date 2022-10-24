function y = Dawson4(x)

% 4-th order rational approximation of the Dawson integral function
% See www.dx.doi.org/10.3247/SL4Soft12.001

y = x*x;
p = 1.0 + y*(0.1107817784 + y*(0.0437734184 ...
        + y*(0.0049750952 + y*(0.0015481656))));
q = 1.0 + y*(0.7783701713 + y*(0.2924513912 ...
        + y*(0.0756152146 + y*(0.0084730365 + 2*0.0015481656*y))));
y = x*(p/q);
