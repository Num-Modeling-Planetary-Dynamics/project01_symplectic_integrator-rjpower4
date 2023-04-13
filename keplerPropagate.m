function xOut = keplerPropagate(gm, dt, x0)
arguments
    gm
    dt (:, 1)
    x0
end
rv0 = x0(:, 1:3);
vv0 = x0(:, 4:6);

rmag0 = sqrt(sum(rv0.^2, 2));
vmag2 = sum(vv0.^2, 2);
h = cross(rv0, vv0, 2);
hmag2 = dot(h, h, 2);
a = 1 ./ (2 ./ rmag0 - vmag2 ./ gm);
ecc = sqrt(1 - hmag2 ./ (gm * a));

n = sqrt(gm ./ a.^3);

if ecc > 1e-12
    E0 = acos(-(rmag0 - a) ./ (a .* ecc));
else
    E0 = 0.0;
end

msk = dot(rv0, vv0, 2) < -1e-12;
E0(msk) = 2 .* pi - E0(msk);


M0 = E0 - ecc .* sin(E0);
M = M0 + n .* dt;

E = zeros(size(M));
for k = 1:size(E, 1)
    ecck = min(length(ecc), k);
    E(k) = oak.tbp.solveKepler(ecc(ecck), M(k));
end
dE = E - E0;

f = a ./ rmag0 .* (cos(dE) - 1.0) + 1.0;
g = dt + (sin(dE) - dE) ./ n;

rv = f .* rv0 + g .* vv0;
rmag = vecnorm(rv, 2);

fdot = -a.^2 ./ (rmag .* rmag0) .* n .* sin(dE);
gdot = a ./ rmag .* (cos(dE) - 1.0) + 1.0;

vv = fdot .* rv0 + gdot .* vv0;

xOut = [rv, vv];
end