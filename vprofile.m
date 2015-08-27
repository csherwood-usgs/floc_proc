% vprofile - Make a profile at the end of the simulation
i=3;
j=3;
kappa=0.41;
netcdf_load('ocean_his.nc')
N = length(s_rho)
nt = length(ocean_time);
[zr]=set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 1, h, squeeze(zeta(:,:,nt)));
[zw]=set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 5, h, squeeze(zeta(:,:,nt)));
[zu]=set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 3, h, squeeze(zeta(:,:,nt)));
%%
h = h(i,j);
zr = h+squeeze(zr(i,j,:))
zw = h+squeeze(zw(i,j,:))
zu = h+squeeze(zu(i,j,:))

AKv = squeeze(AKv(i,j,:,nt))
zwcalc = zw(1)+[0;cumsum(diff(zu))];
AKvi = interp1(zw,AKv,zwcalc)
u = squeeze(u(i,j,:,nt))

dudz = -diff(u)./diff(zu)
taub = bustr(i,j,end)
taus = sustr(i,j,end)
rhoav = mean(rho(:))+1000
tau = [taub; rhoav*AKvi(2:end).*dudz; -taus]

dudzt = (taus/rhoav)./AKv;

%%
ustr = sqrt(abs(taub)/rho0);
K = kappa*ustr*zw
K = kappa*ustr*zw.*(1-zw./h)
up = (ustr./kappa)*log(zu./Zob)
%
figure(1);clf
h3=plot([0;0],[log10(Zob);log10(Zob)],'ok')
set(h3,'markerfacecolor',[0 0 0],'markersize',6)
hold on
h1=plot(abs(u),log10(zu),'linewidth',2)
hold on
h2=plot([0;up],[log10(Zob);log10(zu)],'--k')
%h2=plot([up],[log10(zu)],'--k')
axis([-.001  max(up) -2.5 log10(max(zw))])

xlabel('Speed (m/s)')
ylabel('log_{10} [ z (m) ]')
s1 = '\itu'
s2 = '\itu_* / \kappa ln(\it{z/z_0})'
s3 = '\itz_{ob}'
legend([h1;h2;h3],s1,s2,s3,'location','southeast')
title('Log Profile')
%%
figure(2);clf
subplot(121)
plot(AKv,zw,'.')
hold on
hh1=plot(AKv,zw,'-b')
hh2=plot(K,zw,'--k')
ylabel(' z (m)')
xlabel('AKv (m^2/s)')
title('Eddy Viscosity and Reynolds Stress Profiles')
sa='Model {\itA_{Kv}}'
sb='{\itA_{Kv} = \kappa u_* z} ({\it1-z/h} )'
legend([hh1;hh2],sa,sb)
subplot(122)
plot(tau,zw,'.')
hold on
hh1=plot(tau,zw,'-b')
%plot([taub;0],[0; h],'--k')
plot([taub;-taus],[0; h],'--k')
xlabel('Stress (Pa)')
sa = '{\it\tau = \rho A_{Kv}} d{\itu}/d{\itz}'
legend(hh1,sa)
%%
dudzt = (taub*(1-zw./h)/rhoav)./AKv;
figure(3); clf
hr=plot(dudzt(2:end-1),zw(2:end-1),'-r','linewidth',2)
hold on
hm=plot(dudzt(2:end-1),zwcalc(2:end),'+k','linewidth',2)
hb=plot(dudz,zwcalc(2:end),'--b','linewidth',2)
sb = '\Deltau/\Deltaz_u'
sr = '\tau_b(1-\it{z_w / h})/(\rho \itA_{Kv})'
sm = '\tau_b(1-\it{z_{wcalc} / h})/(\rho \itA_{Kv})'
xlabel( '(s^{-1})')
ylabel('z (m)')
legend([hr;hm;hb],sr,sm,sb)
title('Shear Profile')
%%
figure(4); clf
plot(ocean_time./3600,squeeze(bustr(i,j,:)))
hold on
plot(ocean_time./3600,squeeze(sustr(i,j,:)))
xlabel('Time (h)')
ylabel('Stress {\it{\tau}}-velocity (m^2/^2s)')
%%
figure(5); clf
plot(squeeze(mud_01(i,j,:,end)),zr)

