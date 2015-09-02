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

% defining primary parameters
co=1480;
rho0=1000;
ck=5550;
rhok=2600;
Co=1/(rho0*co^2);
Ck=1/(rhok*ck^2); %compressibility
gk=rhok/rho0;
hk=ck/co;

f=3e6; % frequency

nz1=1e5;
a=logspace(-10,-2,nz1);
x=2*pi*a*f/co;
dx=diff(x); x=x(1:nz1-1); a=a(1:nz1-1);
nz=100; aa=logspace(-9,-3,nz);
xx=2*pi*aa*f/co;

% effective density
rhoex=1e-3./(a.^1); %choosen for illustration

% gamma_o  zeta_o values
beta1=1.02; % gamma_o
beta2=1.02; % zeta_o

% maximum physical value for effective density
zz=find(rhoex>(rhok-rho0)); rhoex(zz)=(rhok-rho0);
% not letting effective density go below zero
zz=find(rhoex<0); rhoex(zz)=0;
g=1+rhoex/(rho0);
zz=find(g<beta1); g(zz)=beta1;
phi=(gk-g)/(gk-1);


%wood expression
vw=((phi*rho0+(1-phi)*rhok).*(phi*Co+(1-phi)*Ck)).^(-0.5);
h=vw/co;
zz=find(h<beta2); h(zz)=beta2;
e=(g.*h.^2);


figure(1),orient tall


% **************************   form function   ************************
% form function f_fi
% g is gamma - this is unnumbered Clay & Medwin eqn on p. 65
kf=2*abs((e-1)./(3*e)+(g-1)./(2*g+1)); 
alpha1=1.2; % called epsilon_1 on p. 85
ffi=(kf.*x.^2)./(1+alpha1*x.^2); % Eqn. 6a

subplot(2,1,1),
loglog(x,ffi,'k'), hold on
xlabel('x_o','fontsize',15), ylabel('f_{ho}','fontsize',15)
axis([1e-4 1e1 1e-7 1e-0])
set(gca,'fontsize',15)
text(5,0.4,'a','fontsize',15)



% ave form function for lognormal distribution f_ho
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
   axf=sum(((x.^2).*(ffi.^2).*dx).*plognorm);
   aflng(jj)=sqrt((axlng(jj).*axf)./ax3);
end

loglog(axlng,aflng,'--k')


hh=legend(' \delta=0.0','\delta=0.5')
set(hh,'fontsize',12)

%% **************************   atten   ************************
% chi_fi scattering
% second Clay & Medwin eqn. on p. 85
kfalfa=2*(((e-1)./(3*e)).^2+(1/3)*((g-1)./(2*g+1)).^2);
chifi=(kfalfa.*x.^4)./(1-1.0*x+1.5*x.^2+kfalfa.*x.^4); % Eqn. 6b

%viscous atten
% chi_hs scattering

rho=g; %rhok/rho0;
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



