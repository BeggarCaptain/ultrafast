function c = metrics_contrast(env, x, z, P)
%METRICS_CONTRAST  Cyst contrast and CNR over the nine anechoic cysts.
%   C   = 20 log10(mu_in/mu_out)  (mu from envelope SQUARED = intensity)
%   CNR = |mu_out-mu_in| / sqrt(var_out+var_in)
%   The nine cysts are the Cartesian product of P.cystX and P.cystZ.

I = double(env).^2;
[XX, ZZ] = meshgrid(x, z);
Cx = cyst_centres(P);
n = size(Cx, 1);
C = zeros(1,n);  CNR = zeros(1,n);
for k = 1:n
    r  = sqrt((XX - Cx(k,1)).^2 + (ZZ - Cx(k,2)).^2);
    in = r <= 2e-3;
    ou = r >= 4e-3 & r <= 6e-3;
    mi = mean(I(in));  mo = mean(I(ou));
    vi = var(I(in));   vo = var(I(ou));
    C(k)   = 20*log10(mi/mo);
    CNR(k) = abs(mo - mi) / sqrt(vo + vi);
end
c.C = C;  c.CNR = CNR;
c.Cmean = mean(abs(C));  c.Cstd = std(abs(C));
c.CNRmean = mean(CNR);   c.CNRstd = std(CNR);
end
