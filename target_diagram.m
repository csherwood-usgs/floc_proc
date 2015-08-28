function [RMSD_star,BIAS,R] = target_diagram( model, data, iplot, ou )
% [RMSD_star, BIAS, R] = target_diagram( model, data, iplot, ou)
%
% Input:
%    model - array or vector of modeled results
%    data  - array or vector of observations, data, or reference
%            (model and data must be same size, shape, and units)
%            (NaNs in either will be ignored)
%    iplot - 0 = no plot (default), 1 = plot
%    ou    - observational uncertainty in normalized units (optional)
% Returns:
%    RMDS_star - signed unbiased RMS difference (eqn 8)
%    BIAS      - normalized bias (eqn 9)
%    R         - correlation coeffcient (eqn 1)

% Jolliff et al., 2009, Summary diagrams for coupled hydrodynamic-ecosystem model
% skill assessment. Journal of Marine Systems 76 (2009) 64-82.
% DOI: 10.1016/j.jmarsys.2008.05.014

% csherwood@usgs.gov
% 27 August 2015
if(exist('iplot','var')~=1),iplot = 0; end
if(exist('ou','var')~=1),ou = NaN; end
ok = find( ~isnan( model(:)+data(:) ) );
N = length(model(ok));
sigm = std(model(ok))
sigd = std(data(ok))
sig_star = sigm/sigd
meanm = mean(model(ok))
meand = mean(data(ok))
BIAS = (meanm - meand)/sigd
R = mean( (model(ok)-meanm).*(data(ok)-meand) )/(sigm*sigd)
RMSD = rms( model(ok)-data(ok) )
RMSD_prime = sqrt( mean( ((model(ok)-meanm)-(data(ok)-meand)).^2 ) )
RMSD_star = sign(sigm-sigd)*sqrt( 1. + sig_star^2 - 2.*sig_star*R )
% RR = corrcoef(model(ok),data(ok))
% C = RR(1,2)
if(iplot)
   plot([-2 2],[0 0],'-k')
   hold on
   plot([0 0],[-2 2],'-k')
   h=circle(1,0,0,'-k');
   if(~isnan(ou))
      h=circle(ou,0,0,'--k');
   end
   h=scatter(RMSD_star,BIAS,54,R,'filled');
   set(h,'markeredgecolor',[.6 .6 .6])
   caxis([-1,1])
   axis([-2 2 -2 2])
   axis square
   xlabel('$$sign(\sigma_{model}-\sigma_{data}) \cdot RMSD^{\prime}/\sigma_{data}$$',...
            'interpreter','latex','fontsize',14)
   ylabel('$$(\overline{model} - \overline{data})/\sigma_{data}$$',...
      'interpreter','latex','fontsize',14)
   ht=text(-1.9,-1.75,'\downarrow Model bias low');
   ht=text(-1.9,-1.9,'\leftarrow Model variance low');
   ht=text(+1.9,+1.9,'Model variance high \rightarrow');
   set(ht,'HorizontalAlignment','right')
   ht=text(+1.9,+1.75,'Model bias high \uparrow');
   set(ht,'HorizontalAlignment','right')
   % xlabel('RMSD*''')
   colorbar
end
   
