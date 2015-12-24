% proc_all - run a series of floc_proc and floc_plot calcs

% clear
% fdyn % load measurements
iplot = 1;
caslist = [ 99; 101; 107 ];
caslist = [ 101; 105; 106; 107 ];
elevlist = [0.5; 1.0; 2.]
for ics = 1:length(caslist)
   cas = caslist(ics)
   casdir = ['case_',num2str(cas)]
   mkdir(casdir)
   cd(casdir)
   logfn = [casdir,'.txt']
   fid = fopen(logfn,'w');
   fprintf(fid,'%s\n',casdir);
   
   floc_proc
   floc_plot
   
   fclose(fid)
   cd('..')
end
   
