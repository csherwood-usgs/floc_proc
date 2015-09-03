% acoustics - Calculate o
global mvco
mvco = 1;
% Need:
%  f = acoustic frequency
%  a = Df/2 - particle radius
%  rhof = particle density
%  M = mass concentration
%  r = range (m)

Df = 0.5e-3
M = 0.1; % 100 mg/l
r = 0.2;
f = 1e6
rhof = 1500;
nsed = length(Df);
a = Df./2;
N = M./ (rhof .* (4/3)*pi.*a.^3);

ffi = zeros(size(a));
chisvo = zeros(size(a));
x = zeros(size(a));
for i = 1:length(a)
   [ffi(i), chisvo(i), x(i)]=f_chi_func( a(i), f, rhof(i) );
end
sumN = sum(N)
ao = sum(a.*N)./sumN
chisvo
chio = (sum( a.*N )/sumN) * (sum(a.^2 .*chisvo .*N)/sumN) / (sum( a.^3 .*N)/sumN)
ffi
fo = sqrt( ((sum(a.*N)/sumN) * (sum(a.^2 .*ffi.^2 .*N)/sumN) ) / (sum( a.^3 .*N)/sumN) )
xo = sum(x.*N)./sumN

rhofo = sum(rhof.*N)/sumN % weighted by number
% or:
rhofo = sum(rhof.*M)./sum(M) % weighted by mass

%% Acoustic response at fixed range r
c = 1500.               % (m/s) speed of sound in water
k = 2*pi*f./c           % (1/m) wave number
ka = k*a  % dimensionless particle size

if(0)
   % old versions of floc scattering functions
   [ff,chi]=gfp( 1e6*a ) % 1 MHz only
   % Need K and Eta for our flocs
   K = ff./(sqrt( ao*rhofo ))
   Eta = 3*chi./(4*ao*rhofo)
end
K = ffi/sqrt(a*rhof)
Eta = 3*chisvo/(r*a*rhof)

% transducer radius (should these be the same for all f?
at=5.0*1e-3;
%system constant from measurement given by R in the paper
meankt=0.0175;
meankt=10

% water attenuation
alphaw=f.*f*31e-15; % alpha_w in the paper at 15C for distilled water
% near field correction
k=2*pi*f/c;
Rstar=r./((k/2).*at.^2);
Psi = 1./(1 - 1./(1+1.35.*Rstar+(2.5.*Rstar).^3.2))'

% Sum water attenuation over all sizes
alphas = r*Eta*sum(M)

att = alphas+r*alphaw;
% backscatter signal, eqn A1.4 in the paper, squared
v=(sum( (((meankt*K.*sqrt(sum(M)))./(r.*Psi)).*exp(-2*att)).^2 ) );
