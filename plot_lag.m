%function plot_lag(x,y,namex,namey,dt,dtunits)
% 
x = rho0*ustar_av.ustrr.^2;
y = ustrcwm;
namex = 'u*cw data';
namey = 'u*cw model';
dt = 1;
dtunits = 'hrs';

figure(5); clf;

ok = find(~isnan(x+y));
x = x(ok);
y = y(ok);
[C,lags]=xcorr(x,y,'coeff');
plot(lags,C)
izero = find(lags==0);
zeroC = C(izero);
[maxC,imaxC]=max(C);
fprintf(1,'C(0) = %5.3f, Cmax = %5.3f at lag = %5.3f %s\n',...
   zeroC,maxC,(imaxC-izero)*dt,dtunits)
if(imaxC-izero)>0
   fprintf(1,'%s leads %s by %5.3f %s\n',namey,namex,(imaxC-izero)*dt,dtunits);
elseif(imaxC-izero)<0
   fprintf(1,'%s leads %s by %5.3f %s\n',namex,namey,-(imaxC-izero)*dt,dtunits);
end
%%
x = rho0*ustar_av.ustrc.^2;
y = ustrc;
dt = 1;
dtunits = 'hrs';
namex = 'u*c data';
namey = 'u*c model';
ok = find(~isnan(x+y));
x = x(ok);
y = y(ok);
[C,lags]=xcorr(x,y,'coeff');
hold on
plot(lags,C)
ylabel('Cross-correlation \itr')
xlabel(['Lag (',dtunits,')'])
axis([-24 24 0 1])
izero = find(lags==0);
zeroC = C(izero);
[maxC,imaxC]=max(C);
fprintf(1,'C(0) = %5.3f, Cmax = %5.3f at lag = %5.3f %s\n',...
   zeroC,maxC,(imaxC-izero)*dt,dtunits)
if(imaxC-izero)>0
   fprintf(1,'%s leads %s by %5.3f %s\n',namey,namex,(imaxC-izero)*dt,dtunits);
elseif(imaxC-izero)<0
   fprintf(1,'%s leads %s by %5.3f %s\n',namex,namey,-(imaxC-izero)*dt,dtunits);
end