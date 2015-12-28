% proc_all - run a series of floc_proc and floc_plot calcs

clear
fdyn % load measurements
iplot = 0;
topdir = pwd;
% caslist = [ 99; 101; 105; 106; 107 ];
caslist = [99; 101; 103; 104; 105; 107; 108; 109; 110 ];
elevlist = [0.5; 1.0; 2.]
fid = fopen('sens.txt','a+');

for ics = 1:length(caslist)
   cas = caslist(ics)
   fprintf(fid,'%d ,',cas);
   casdir = ['case_',num2str(cas)]
   mkdir(casdir)
   cd(casdir)
   
   floc_proc
   floc_plot

   cd(topdir)
end
fclose(fid)

