% test_acoustics_v - Script to test v returned from acoustics function
nf = 2.0
Dp0 = 4.e-6
Dmax = 1500.e-6
Df = logspace( log10(40e-6), log10(Dmax), 5)'

rhow = 1000;
rhos = 2650;
f = 3.e6;
r = .2;
rhof = max(rhos+20,min(rhos,rhow + (rhos-rhow)*(Df./Dp0).^(nf-3)));
rhof = 2650*ones(size(rhof));
M = .1*ones(size(Df))./length(Df);
v = acoustics( Df, M, rhof, r, f )