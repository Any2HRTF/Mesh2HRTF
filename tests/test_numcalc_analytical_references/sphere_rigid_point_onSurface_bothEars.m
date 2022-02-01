% [1], Section 12.4
% AKtools revision 191
%
% [1] L. Beranek and T. Mellow, Acoustics. Sound fields, transducers and
%     vibration, 2nd ed. London et al.: Academic Press, 2019.
%
% Fabian Brinkmann 11/2020
close all; clear; clc

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
xyz_src_leftear  = [-0.0875 0 0];
xyz_src_rightear = [ 0.0875 0 0];
% point source distance (equal for left and right)
d = sqrt(sum(xyz_src_leftear.^2));
% sphere radius
a = 0.0875;
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

% incident fields (Eq. 12.24) ----------------------------------------------
r1_leftear = sqrt((XYZ(:,1) - xyz_src_leftear(1)).^2 + ...
               (XYZ(:,2) - xyz_src_leftear(2)).^2 + ...
               (XYZ(:,3) - xyz_src_leftear(3)).^2);
      
p_i_leftear = 1j * k * rho_0 * c * U0 * exp(-1j * k * r1_leftear) ./ ...
           (4 * pi * r1_leftear);


r1_rightear = sqrt((XYZ(:,1) - xyz_src_rightear(1)).^2 + ...
                (XYZ(:,2) - xyz_src_rightear(2)).^2 + ...
                (XYZ(:,3) - xyz_src_rightear(3)).^2);
 
p_i_rightear = 1j * k * rho_0 * c * U0 * exp(-1j * k * r1_rightear) ./ ...
            (4 * pi * r1_rightear);
        
        
r1_bothears = r1_leftear + r1_rightear;

p_i_bothears = 1j * k * rho_0 * c * U0 * exp(-1j * k * r1_bothears) ./ ...
            (4 * pi * r1_bothears);

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
p_total_leftear  = p_i_leftear  + p_s;
p_total_rightear = p_i_rightear + p_s;
p_total_bothears = p_i_bothears + p_s;

save('ref_rigid_leftear',  'p_total_leftear',  'f_Hz', 'c', 'a', 'R', 'AZ', 'XYZ')
save('ref_rigid_rightear', 'p_total_rightear', 'f_Hz', 'c', 'a', 'R', 'AZ', 'XYZ')
save('ref_rigid_bothears', 'p_total_bothears', 'f_Hz', 'c', 'a', 'R', 'AZ', 'XYZ')

%%
AKf(20,24)
subplot(3,2,1)
scatter(XYZ(:,1), XYZ(:,2), 75, db(p_total_leftear), 'Marker', '.')
cb = colorbar;
ylabel(cb, 'pressure in dB (re 1, normalized)')
axis equal
grid off
xlabel 'x in m'; ylabel 'y in m'
title 'Point source on left ear'

subplot(3,2,2)
scatter(XYZ(:,1), XYZ(:,2), 75, db(p_total_rightear), 'Marker', '.')
cb = colorbar;
ylabel(cb, 'pressure in dB (re 1, normalized)')
axis equal
grid off
xlabel 'x in m'; ylabel 'y in m'
title 'Point source on right ear'

%%
subplot(3,2,3)
scatter(XYZ(:,1), XYZ(:,2), 75, db(sum(legendre_term,2)), 'Marker', '.')
cb = colorbar;
ylabel(cb, 'pressure in dB (re 1, normalized)')
axis equal
grid off
xlabel 'x in m'; ylabel 'y in m'
title 'Legendre term of scattered field'

subplot(3,2,4)
scatter(XYZ(:,1), XYZ(:,2), 75, db(p_s), 'Marker', '.')
cb = colorbar;
ylabel(cb, 'pressure in dB (re 1, normalized)')
axis equal
grid off
xlabel 'x in m'; ylabel 'y in m'
title 'Scattered field'

%%
subplot(3,2,5)
scatter(XYZ(:,1), XYZ(:,2), 75, db(p_i_leftear), 'Marker', '.')
cb = colorbar;
ylabel(cb, 'pressure in dB (re 1, normalized)')
axis equal
grid off
xlabel 'x in m'; ylabel 'y in m'
title 'Left ear incident sound field'

subplot(3,2,6)
scatter(XYZ(:,1), XYZ(:,2), 75, db(p_i_rightear), 'Marker', '.')
cb = colorbar;
ylabel(cb, 'pressure in dB (re 1, normalized)')
axis equal
grid off
xlabel 'x in m'; ylabel 'y in m'
title 'Right ear incident sound field'