function rotmats = rotationMatrices(x)
rs = x(:, 1:3);
vs = x(:, 4:6);

rmags = sqrt(sum(rs.^2, 2));
vmags = sqrt(sum(vs.^2, 2));

rhats = rs ./ rmags;
vhats = vs ./ vmags;
hhats = cross(rhats, vhats, 2);

rotmats = zeros(3, 3, size(rs, 1));
rotmats(1, :, :) = permute(rhats, [3, 2, 1]);
rotmats(3, :, :) = permute(hhats, [3, 2, 1]);
rotmats(2, :, :) = permute(cross(hhats, rhats, 2), [3, 2, 1]);
end