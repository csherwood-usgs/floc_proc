dgrab = 1e-3*dmm(-2:1:11);
% lab results from grab samples (pd3.xls)
% last two rows are switched so sample with fines can be pulled out more
% easily
frgrab =[0	0	0.09	0.63	6.4	69.45	19.19	2.3	0.89	0.42	0.31	0.18	0.09	0.04;...
0.29	0.12	0.37	1.38	7.51	67.62	18.57	1.95	1.06	0.47	0.31	0.2	0.1	0.04;...
0	0	0.07	0.75	7.62	68.26	16.76	1.81	1.36	1.22	0.96	0.68	0.35	0.16;...
0	0	0.12	1.78	13.55	66.27	13.4	0.97	0.92	1.02	0.91	0.61	0.3	0.15;...
0	0	0.59	6.37	9	63.77	15.09	0.65	1.18	1.19	1.12	0.64	0.28	0.12;...
0	0	0.11	1.04	7.33	67.04	18.65	3.12	1.54	0.58	0.28	0.18	0.08	0.04];

% sum(frgrab,2) % close enough
dgrab = fliplr(dgrab);
frgrab = fliplr(frgrab);
dfr = 1e3*dgrab';
fr = mean(frgrab(1:5,:))';
%
di = [.03 .0625 .125 .25]'
fi = interp1(dfr,fr,di);
fi = 100*fi./sum(fi)
fprintf(1,'D (mm) f\n')
disp( [dfr, fr] );


figure(1); clf
subplot(211)
h1=semilogx(1e3*dgrab,frgrab,'linewidth',1.5,'color',[.5 .5 .5]);
hold on
h2=semilogx( 1e3*dgrab,mean(frgrab(1:5,:)),'linewidth',2.5,'color',[.2 .2 .2]);
h3=plot(di,fi,'or');
set(h3,'markersize',9,'markerfacecolor',[.7 .7 .7]);
legend([h1(1);h2;h3],'All Grabs','Mean','Model')
grid on
xlim([0.01 1])
ylabel('Percent')

subplot(212)
h1=semilogx( 1e3*dgrab,cumsum(frgrab,2),'linewidth',1.5,'color',[.5 .5 .5]);
hold on
h2=semilogx( 1e3*dgrab,cumsum( mean(frgrab(1:5,:)) ),'linewidth',2.5,'color',[.2 .2 .2]);
grid on
xlim([0.01 1])
subplot(212)
h3=plot(di,cumsum(fi),'-r','linewidth',2);
xlabel('Grain Diameter (mm)')
ylabel('Percent Finer')
legend([h1(1);h2;h3],'All Grabs','Mean','Model')

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
     'PaperPosition',newpos)
print('-dpng', 'model_size_dist.png', '-r300');
drawnow
set(gcf,'Units',oldscreenunits,...
     'PaperUnits',oldpaperunits,...
     'PaperPosition',oldpaperpos)

% calculate critical shear stress and settling velocity
p = soulsby_particle(1e-3*di,2650,1025,1.5e-6)
%% ROMS-like input
fprintf(1,'Initial fractional distribution.\n\n')
for i = 1:length(di), fprintf(1,' %6.2fd0',fi(i)),end
fprintf(1,'\n\n')

fprintf(1,'! Median sediment grain diameter (mm).\n\n') 
fprintf(1,' SAND_SD50 == ')
for i = 1:length(di), fprintf(1,' %6.4fd0',di(i)),end
fprintf(1,'\n\n')
fprintf(1,'! Sediment concentration (kg/m3).\n\n')
fprintf(1,' SAND_CSED == ')
for i = 1:length(di), fprintf(1,' %6.4fd0',0.),end
fprintf(1,'\n\n')
fprintf(1,'! Sediment grain density (kg/m3).\n\n')
fprintf(1,' SAND_SRHO == ')
for i = 1:length(di), fprintf(1,' %6.1fd0',2650.),end
fprintf(1,'\n\n')
fprintf(1,'! Particle settling velocity (mm/s).\n\n')
fprintf(1,' SAND_WSED == ')
for i = 1:length(di), fprintf(1,' %6.4fd0',1e3*p.ws(i)),end
fprintf(1,'\n\n')
fprintf(1,'! Surface erosion rate (kg/m2/s).\n\n')
fprintf(1,' SAND_ERATE == ')
for i = 1:length(di), fprintf(1,' %s ','5.0d-4'),end
fprintf(1,'\n\n')
fprintf(1,'! Critical shear for erosion and deposition (N/m2).\n\n')
fprintf(1,' SAND_TAU_CE == ')
for i = 1:length(di), fprintf(1,' %6.4fd0',p.tau_crit(i)),end
fprintf(1,'\n')
fprintf(1,' SAND_TAU_CD == ')
for i = 1:length(di), fprintf(1,'   %s ','0.0d0'),end
fprintf(1,'\n\n')








