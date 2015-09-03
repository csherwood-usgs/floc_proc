% floc density - Compare estimates for fractal dimensions with power law
nf = [1.9 2.0 2.2]
Dp0 = 4.e-6
Dmax = 1500.e-6
Df = logspace( log10(1.e-6), log10(Dmax), 61)

rhow = 1000
rhos = 2650
co=1480;

rhof = NaN*ones(length(nf),length(Df));
for i=1:3
rhof(i,:) = max(1020,min(rhos,rhow + (rhos-rhow)*(Df./Dp0).^(nf(i)-3)));
end
% Estimated floc densities of MacDonald et al., JGR 118:2581
delrho = 5.829e-4*(Df./2).^(-1.1238)
delrho(delrho>(rhos-rhow))=rhos-rhow;
delrho(delrho<20)=20;
% Second, simpler estimate from Thorne et al. 2014 CSR paper, below Fig 2.
delrho2 = 0.001*(Df./2).^(-1.)
delrho2(delrho2>(rhos-rhow))=rhos-rhow;
delrho2(delrho2<20)=20;

figure(1); clf
hh=semilogx( Df./2, rhof-rhow );
set(hh,'linewidth',2)
hold on
hk=semilogx( Df./2, delrho,'-k')
set(hk,'color',[.4 .4 .4])
set(hk,'linewidth',2)
hk2=semilogx( Df./2, delrho2,'-k')
set(hk2,'color',[.2 .2 .2])
set(hk2,'linewidth',2)
ylabel('Effective density \rho_e (kg/m^3)')
xlabel('Particle radius {\ita} (m)')
legend([hh;hk;hk2],'nf=1.9','nf=2','nf=2.2','0.0006/(a^{1.1})','0.001/a')
%%
f = [1. 2.5 4].*1e6
for i=1:length(f)
   k = 2*pi*f(i)/co
end