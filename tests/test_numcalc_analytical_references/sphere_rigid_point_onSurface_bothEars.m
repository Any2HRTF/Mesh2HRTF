% [1], Section 12.4
% AKtools revision 191
%
% [1] L. Beranek and T. Mellow, Acoustics. Sound fields, transducers and
%     vibration, 2nd ed. London et al.: Academic Press, 2019.
%
% Fabian Brinkmann 11/2020
clear; clc

% speed of sound
c = 343;
% density of air
rho_0 = 1.1839;
% frequency
f_Hz = 4e3;
% wave number
k = 2*pi*f_Hz / c;
% Point source volume velocity
U0 = 1;
% point source position [x y z] in m
xyz_src = [.0875 0 0];


% point source distnace
d = sqrt(sum(xyz_src.^2));
% sphere radius
a = .0875;
% sphere surface area
S = 4 * pi * a^2;
% max SH order
N_sh = 100;

% points at which he sound field is evaluated
r = (.09 : .01 : .5)';
az = (0:2:360-1)';

R = repelem(r, numel(az), 1);
AZ = repmat(az, [numel(r) 1]);

% Left ear grid for which the sound field is calculated
XYZ = [cosd(AZ).*R sind(AZ).*R zeros(size(R))];

% calculate the sound field
kr = k * R;
ka = k * a;

n = 0:N_sh;

% incident field (Eq. 12.24) ----------------------------------------------
r1 = sqrt((XYZ(:,1) - xyz_src(1)).^2 + ...
          (XYZ(:,2) - xyz_src(2)).^2 + ...
          (XYZ(:,3) - xyz_src(3)).^2);

p_i = 1j * k * rho_0 * c * U0 * exp(-1j * k * r1) ./ (4 * pi * r1);

% scattered field (Eq. 12.34) ---------------------------------------------
pre_factor = -k^2*a^2 * rho_0*c*U0 / S * (2*n + 1);

radial_term = AKshRadial(k*d, 'hankel', 2, n) ...
    .* ( ...
    AKshRadial(ka, 'bessel', 1, n, 'derived') ...
    ./ ...
    AKshRadial(ka, 'hankel', 2, n, 'derived') ...
    ) ...
    .* ...
    AKshRadial(kr, 'hankel', 2, n);

legendre_term = nan(numel(az), numel(n));
for nn = 0:N_sh
    legendre_tmp = legendre(nn, cosd(az));
    legendre_term(:, nn+1) = legendre_tmp(1,:).';
end
legendre_term = repmat(legendre_term, [numel(r) 1]);

p_s = sum(pre_factor .* radial_term .* legendre_term, 2);

% total field (Eq. 12.27) -------------------------------------------------
p = p_i + p_s;

% mirror sound field for right ear source ---------------------------------
p_left = zeros(size(p));

n_az = numel(az);

for rr = 1:numel(r)
    id_a = (rr-1)*n_az+1;
    id_b = rr*n_az;
    p_left(id_a:id_b) = circshift(p(id_a:id_b), n_az/2);
end

% save sound pressure for each source type --------------------------------
p_total_rightear = p;
p_total_leftear = p_left;
p_total_bothears = cat(3, p_left, p);

save('ref_rigid_leftear',  'p_total_leftear',  'f_Hz', 'c', 'a', 'R', 'AZ', 'XYZ')
save('ref_rigid_rightear', 'p_total_rightear', 'f_Hz', 'c', 'a', 'R', 'AZ', 'XYZ')
save('ref_rigid_bothears', 'p_total_bothears', 'f_Hz', 'c', 'a', 'R', 'AZ', 'XYZ')

%%

% AKf(20,8)
% 
% c_lim = prctile(db(abs(p_total_bothears(:,:,1))), 95);
% 
% for nn = 1:2
%     subplot(1,2,nn)
%     scatter(XYZ(:,1), XYZ(:,2), 75, db(p_total(:,:,nn)), 'Marker', '.')
%     cb = colorbar;
%     set(gca, 'CLim', [c_lim-40 c_lim + 20])
%     ylabel(cb, 'pressure in dB (re 1, normalized)')
%     axis equal
%     grid off
%     xlabel 'x in m'; ylabel 'y in m'
% end