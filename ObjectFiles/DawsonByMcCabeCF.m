function y = DawsonByMcCabeCF(z)

% Dawson integral function via the McCabe's continued fraction
% See www.dx.doi.org/10.3247/SL4Soft12.001

a2  = 2*z*z;

am1  = 1;
bm1  = 1+a2;
a    = -2*a2;
b    = 3+a2;
newa = b;
newb = b*bm1 + a;
fn   = newa/newb;
an   = newa;
bn   = newb;
n   = 2;
    
while n < 100
    n = n + 1;
    a = a - 2*a2;
    b = b + 2;
    newa = b*an + a*am1;
    newb = b*bn + a*bm1;
    f = newa/newb;
    if f==fn;
        break
    end
    fn  = f;
    am1 = an;
    bm1 = bn;
    an  = newa;
    bn  = newb;
end

y = z*f;
