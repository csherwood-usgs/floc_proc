res = load('sens.txt')
cas = res(:,1);
s1 = res(:,2:4);
BIAS = res(:,5)
RMSD_star = res(:,6);
R = res(:,7);
s2 = res(:,8:10);
BIAS2 = res(:,11);
RMSD_star2 = res(:,12);
R2 = res(:,13);
ou = NaN;

fconcen = res(:,18);
fsize = res(:,19);

c = [.7 .4 .4];
plot([-2 2],[0 0],'-k');
   hold on
   plot([0 0],[-2 2],'-k');
   h=circle(1,0,0,'-k');
   if(~isnan(ou))
      h=circle(ou,0,0,'--k');
   end
   h=scatter(RMSD_star,BIAS,54,R,'filled');
   %get(h)
   c = [.7 .4 .4];
   set(h,'markeredgecolor',c,'linewidth',1.5);
      h=scatter(RMSD_star2,BIAS2,54,R2,'filled');
   %get(h)
   c = [.4 .4 .8]
   set(h,'markeredgecolor',c,'linewidth',1.5);
   caxis([0,.5]);

   axis([-2 2 -2 2]);

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
