% [1], Section 6.10.3
% AKtools revision 191
%
% [1] E. G. Williams, Fourier Acoustics. Sound radiation and nearfield
%     acoustical holography, 1st ed. Academic Press, 1999.
%
% Fabian Brinkmann 11/2020
close all; clear; clc

% speed of sound
c = 343;
% frequency
f_Hz = 4e3;
% wave number
k = 2*pi*f_Hz / c;
% plane wave amplitude
phase = 0;
P0 = cosd(phase) + 1j*sind(phase);
% sphere radius
a = .0875;
% max SH order
N_sh = 100;

% points at which he sound field is evaluated
r = (.09 : .01 : .5)';
az = (0:2:360-1)';
az = mod(az+90, 360); % to match plane wave direction in simulation

R = repelem(r, numel(az), 1);
AZ = repmat(az, [numel(r) 1]);

XYZ = [cosd(AZ).*R sind(AZ).*R zeros(size(R))];

% calculate the sound field
kr = k * R;
ka = k * a;

n = 0:N_sh;

% single terms
pre_factor = (2*n + 1) .* 1j.^n;
radial_term = AKshRadial(kr, 'bessel', 1, n) ...
              - ( ...
                  AKshRadial(ka, 'bessel', 1, n, 'derived') ...
                  ./ ...
                  AKshRadial(ka, 'hankel', 1, n, 'derived') ...
                 ) ...
              .* ...
              AKshRadial(kr, 'hankel', 1, n);


legendre_term = nan(numel(az), numel(n));
for nn = 0:N_sh
    legendre_tmp = legendre(nn, cosd(az));
    legendre_term(:, nn+1) = legendre_tmp(1,:).';
end
legendre_term = repmat(legendre_term, [numel(r) 1]);

% combine
p_total = P0 * sum(pre_factor .* radial_term .* legendre_term, 2);

% align grid coordinates with simulation results
XYZ(:,1:2) = [XYZ(:,2), -XYZ(:,1)];

save('ref_rigid_plane', 'p_total', 'f_Hz', 'c', 'a', 'R', 'AZ', 'XYZ')
