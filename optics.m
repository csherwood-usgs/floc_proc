% optics.m
% calculate area of suspended particles

% d = diameter
% m = mass concentration
% r = density
% A = area (summed over sizes)

% N = total volume per m3 = mass/m3 / density
%     -------------------
%     volume per particle
%
% A = N*area per particle 

% area of flocs
d= squeeze(repmat(fdiam(1:NST),1,nz,nt));
r= squeeze(repmat(rhos(1:NST),1,nz,nt));
Af =  squeeze( sum ((3./4.)*( m ./r)./(0.5*d) ));
clear d r

figure; clf
pcolorjw( s2d*tz, h+z_w, Af)
colorbar
caxis([0 100])
title('Floc Optical Response (\mum)')

% area of sand grains
d= squeeze(repmat(fdiam(NST+NNN+1:NST+NNN+NND),1,nz,nt));
r= squeeze(repmat(rhos(NST+NNN+1:NST+NNN+NND),1,nz,nt));
As =  squeeze( sum ((3./4.)*( snd ./r)./(0.5*d) ));
clear d r

figure; clf
pcolorjw( s2d*tz, h+z_w, As)
colorbar
caxis([0 100])
title('Sand Optical Response (\mum)')
