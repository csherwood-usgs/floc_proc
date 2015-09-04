% optical response
Dfv = fdiam(1:NST);
r_ac9f=zeros(nz,nt);

if(exist('isfloc','var')~=1),isfloc=1;,end
% import D, nf, c_ac9, c_mass, c_lisst, and same for sand
if(exist('c_ac9','var')~=1),
   load c_coeffs
   load c_coeffs_sand
end
% skip next step because nf=2 is column 1
% c_ac9_nf = interp1(nf',c_ac9',fnf)';

c_ac9_i = interp1(D*1e-6',c_ac9(:,1),Dfv')';
r_ac9=sum(mv.*c_ac9_i);
for jj=1:nt
   for ii=1:nz
      mv = squeeze(m(:,ii,jj));
      r_ac9f(ii,jj)=sum(mv.*c_ac9_i);
   end
end