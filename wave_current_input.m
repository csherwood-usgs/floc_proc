clear
hA=12;

load('ustar_av')
% ustar_av.dn=ustar_av.dn(2:end);
% ustar_av.u=ustar_av.u(2:end);
% ustar_av.v=ustar_av.v(2:end);
% ustar_av.ubr=ustar_av.ubr(2:end);
% ustar_av.Tr=ustar_av.Tr(2:end);
% ustar_av.ustrc=ustar_av.ustrc(2:end);
% ustar_av.ustrr=ustar_av.ustrr(2:end);
% ustar_av.ustwm=ustar_av.ustwm(2:end);
% ustar_av.zoa=ustar_av.zoa(2:end);
%ustar_av.dn(1)=ustar_av.dn(2);
ustar_av.u(1)=ustar_av.u(2);
ustar_av.v(1)=ustar_av.v(2);
ustar_av.ubr(1)=ustar_av.ubr(2);
ustar_av.Tr(1)=ustar_av.Tr(2);
ustar_av.ustrc(1)=ustar_av.ustrc(2);
ustar_av.ustrr(1)=ustar_av.ustrr(2);
ustar_av.ustwm(1)=ustar_av.ustwm(2);
ustar_av.zoa(1)=ustar_av.zoa(2);

AwavH=.7*ones(size(ustar_av.Tr))+.4*((ustar_av.ubr-mean(ustar_av.ubr))/std(ustar_av.ubr));

[ub,Tbav] = ubspecfun2(AwavH,ustar_av.Tr,hA);
for i=1:length(ustar_av.dn)
    rl(i) = rwavelength( ustar_av.Tr(i), hA, abs(complex(ustar_av.u(i),ustar_av.v(i))), 20*pi/180 );
end

for i=1:length(ustar_av.dn)
    m =  soulsby_mvco( ub(i), 2*pi./ustar_av.Tr(i), abs(complex(ustar_av.u(i),ustar_av.v(i))), 1, 20*pi/180,ustar_av.zoa(i)  );
    ustrwm(i)=m.ustrwm;
    ustrw(i)=m.ustrw;
    ustrc(i)=m.ustrc;
end
%%
figure(2);clf
subplot 211
h1=line(ustar_av.dn,ustrc);
h2=line(ustar_av.dn,ustrw,'Color','r');axtt
set(gca,'XTick',datenum(2013,12+[-2:6],1),'TickDir','out');grid
datetick('x','mmmyy','keeplimits','keepticks');title('Bottom stress');legend('current','wave',0)
subplot 212
line(ustar_av.dn,ustrwm);
line(ustar_av.dn,ustar_av.ustwm,'Color','r');axtt
set(gca,'XTick',datenum(2013,12+[-2:6],1),'TickDir','out');grid
datetick('x','mmmyy','keeplimits','keepticks');title('wave-current Bottom stress')
%%
figure(3);clf
h1=line(ustar_av.dn,ub);
line(ustar_av.dn,ustar_av.ubr,'Color','r');axtt
set(gca,'XTick',datenum(2013,12+[-2:6],1),'TickDir','out');grid
datetick('x','mmmyy','keeplimits','keepticks');title('orbital velocity')
legend('recalc','obs')
%print -dpng -painters WC_bottom_stress_MVCO.png
%%
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
[ix iy]=size(x_rho);
Lm=6;       %<--- else put size of grid here, from mod_param.F
Mm=5;        %<--- else put size of grid here, from mod_param.F
LP = Lm+2;    %don't change this.
MP = Mm+2;    %don't change this.


Fname=['wave_force_MVCO.nc'];

wtime=ustar_av.dn-ustar_av.dn(1);
Ub=repmat(ustar_av.ubr,[1,MP,LP]);
Dwave=repmat(20,[1,MP,LP]);
Pwave=repmat(ustar_av.Tr,[1,MP,LP]);
Hwave=repmat(AwavH,[1,MP,LP]);
Lwave=repmat(rl',[1,MP,LP]);

write_roms_wave_ub_forcing_mvco(LP,MP,Fname,wtime);

ncwrite(Fname,'wave_time',wtime);
ncwrite(Fname,'Hwave',permute(Hwave,[3,2,1]));
ncwrite(Fname,'Dwave',permute(Dwave,[3,2,1]));
ncwrite(Fname,'Pwave_bot',permute(Pwave,[3,2,1]));
ncwrite(Fname,'Lwave',permute(Lwave,[3,2,1]));
ncwrite(Fname,'Uwave_rms',permute(Ub,[3,2,1]));



return
%%
