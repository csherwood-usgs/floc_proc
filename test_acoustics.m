% test_acoustics - Script to test f_chi_func routine
global thorne
thorne = 1; % sets physical params to match Thorne's paper
f=3e6; % frequency
a = logspace(-9, -3, 200 )'; % radius (m)
Df = 2*a;                    % diameter (m)
if(0) % hybrid - this produces same intrinsic results as Thornes code and Fig. 3
   rhow = 1000;
   % effective density according to Macdonald
   Cf = 0.001; % kg /m2
   m = 1.;
   rhoef = Cf./(a.^m); %choosen for illustration
   rhoef(rhoef<20)=20;
   rhof = rhoef+rhow;
   rhof(rhof>2600)=2600;
elseif(0) % solid
   rhow = 1000;
   rhof = 2600*ones(size(a));
elseif(0) % water
   rhow = 1000;
   rhof = rhow*ones(size(a));
elseif(1) % default to MVCO parameters
   thorne = 0;
   rhow = 1025;
   rhos = 2650;
   nf = 1.9;
   Dp0 = 4.e-6
   rhof = max(rhow+20,min(rhos,rhow + (rhos-rhow)*(Df./Dp0).^(nf-3))); 
end

rhoef = rhof-rhow;
gam = 1.+rhoef/rhow;

ffi = zeros(size(a));
chisvo = zeros(size(a));
x = zeros(size(a));

for i = 1:length(a)
   [ffi(i), chisvo(i), x(i)]=f_chi_func( a(i), f, rhof(i) );
end

figure(1); hold on
% subplot(311)
% semilogx(x,gam)
% xlim([1e-3 1])

subplot(211)
loglog( x, ffi )
% ylim([1e-6 10^1.2])
% xlim([1e-3 1])

subplot(212)
loglog( x, chisvo )
% xlim([1e-3 1])
% ylim([1e-8 10^1.2])
