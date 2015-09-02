% model_floc_paper_crs.m
% CRS comments on original code from:

%This software is in support of the paper
% 'Modelling acoustic scattering by suspended flocculating sediments'
%by Peter D. Thorne1, Iain T. MacDonald and Christopher E. Vincent in Continental Shelf Research

% The program is for educational use and not for comercial purposes.
% The code is unsupported.

%1 The form function is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%2 The average form function is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%3 The normalise total scattering and viscous cross-section is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%4 The average normalise total scattering and viscous cross-section is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%5 Plots of the results are given


close all
clear all

f=3e6; % frequency
aolist = logspace(-9, -3, 10 );
rhoflist = 1000+1e-3./(aolist.^1); %choosen for illustration
dsigma=0.5;
for i=1:length(aolist)
   ao = aolist(i);
   rhof = rhoflist(i);
% Input
%    f - frequency in Hz
%    ao - average particle radius
%    rhof - particle density
%    dsigma - width of lognormal distribution

% define primary parameters
co=1480;   % speed of sound in water
rho0=1000; % density of water
ck=5550;   % speed of sound in particle
rhok=2600; % density of particle
Co=1/(rho0*co^2);
Ck=1/(rhok*ck^2); %compressibility
gk=rhok/rho0;
hk=ck/co;

nz1=1e5;
a=logspace(-10,-2,nz1);
x=2*pi*a*f/co;
dx=diff(x);
% trim to match lenght of dx:
x=x(1:nz1-1); a=a(1:nz1-1);

xo=2*pi*ao*f ./co;

% effective density
rhoex = rhof-rho0

% gamma_o  zeta_o values
beta1=1.02; % gamma_o
beta2=1.02; % zeta_o

% maximum physical value for effective density
zz=find(rhoex>(rhok-rho0)); rhoex(zz)=(rhok-rho0);
% not letting effective density go below zero
zz=find(rhoex<0); rhoex(zz)=0;
gam=1+rhoex/(rho0);
zz=find(gam<beta1); gam(zz)=beta1;
phi=(gk-gam)/(gk-1);

% Wood's expression - Eqn 9
vw=((phi*rho0+(1-phi)*rhok).*(phi*Co+(1-phi)*Ck)).^(-0.5);
h=vw/co;
zz=find(h<beta2); h(zz)=beta2;

% **************************   form function   ************************
% form function f_fi
% unnumbered Clay & Medwin eqn on p. 65
e=(gam.*h.^2);
kf=2*abs((e-1)./(3*e)+(gam-1)./(2*gam+1)); 
alpha1=1.2; % called epsilon_1 on p. 85
ffi=(kf.*x.^2)./(1+alpha1*x.^2); % Eqn. 6a

% Average form function for lognormal distribution f_ho

% of all x around each xo
   sigma=dsigma*xo;
   mu=log((xo.^2)./sqrt(xo.^2+sigma.^2));
   sigman=sqrt(log((sigma/xo).^2+1));
   plognorm=(1./(x*sigman*sqrt(2*pi))).*exp(-((log(x)-mu).^2)/(2*sigman.^2));
   axlng(i)=sum(x.*plognorm.*dx);
   ax3=sum((x.^3).*plognorm.*dx);
   axf=sum(((x.^2).*(ffi.^2).*dx).*plognorm);
   aflng(i)=sqrt((axlng(i).*axf)./ax3);

end

figure(1),orient tall
subplot(2,1,1),
%loglog(xlist,ffilist,'k'), hold on
xlabel('x_o','fontsize',15), ylabel('f_{ho}','fontsize',15)
axis([1e-4 1e1 1e-7 1e-0])
set(gca,'fontsize',15)
text(5,0.4,'a','fontsize',15)
loglog(axlng,aflng,'--k')
hh=legend(' \delta=0.0','\delta=0.5')
set(hh,'fontsize',12)

%% **************************   atten   ************************
% chi_fi scattering
% second Clay & Medwin eqn. on p. 85
kfalfa=2*(((e-1)./(3*e)).^2+(1/3)*((gam-1)./(2*gam+1)).^2);
chifi=(kfalfa.*x.^4)./(1-1.0*x+1.5*x.^2+kfalfa.*x.^4); % Eqn. 6b

%viscous atten
% chi_hs scattering

rho=gam; %rhok/rho0;
v=1e-6; c=1480;

da=c*dx/(2*pi*f);
a=a(1:nz1-1);
w=2*pi*f;
k=2*pi*f/c;
beta=sqrt(w/(2*v));
delta=0.5*(1+(9./(2*beta*a)));
s=(9./(4*beta*a)).*(1+(1./(beta*a)));
e1=((k*(rho-1).^2)/2);
e2=s./(s.^2+(rho+delta).^2);

%converted to chi_hv values
e12=e1.*e2;
chiv=(4*a/3).*e12;



%combine scat and visc atten chi_h
chisvo=chifi+chiv;
subplot(2,1,2)
loglog(x,chisvo,'k'), hold on
axis([1e-4 1e1 1e-7 1e-0])
xlabel('x_o','fontsize',mn), ylabel('\chi_{ho}','fontsize',mn)
set(gca,'fontsize',mn)
%title('Normalised total scattering cross-section  - attenuation','fontsize',mn)
text(5,0.4,'b','fontsize',mn)


% ave scattering for lognormal distribution chi_ho
dsigma=0.5;
%Lognormal
for jj=1:nz;
   xo=xx(jj);
   sigma=dsigma*xo;
   mu=log((xo.^2)./sqrt(xo.^2+sigma.^2));
   sigman=sqrt(log((sigma/xo).^2+1));
   plognorm=(1./(x*sigman*sqrt(2*pi))).*exp(-((log(x)-mu).^2)/(2*sigman.^2));
   axlng(jj)=sum(x.*plognorm.*dx);
   ax3=sum((x.^3).*plognorm.*dx);
   axchi=sum((x.^2).*chisvo.*plognorm.*dx);
   achilngo(jj)=(axlng(jj).*axchi)./ax3;
end

loglog(axlng,achilngo,'--k')



