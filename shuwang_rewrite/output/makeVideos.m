addpath ..
%set(gcf,'Position',get(0,'Screensize'));
set(0,'DefaultAxesFontSize',22)
options.savefig = false;
clf;

irate = 1000; % controls the speed of the visualization

%file = 'Chi0_ra065_beta1_conc0p3.bin';
%file = 'Chi4_ra07564_beta0_conc0p48_MoreEL.bin';
%file = 'Chi5_ra065_beta0p1_conc0.bin';
file = 'Chi200_RA0p85_Conc0p3_Beta0_y0p1_eps0p04_n20.bin';
 [posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1] = loadFile(file);
  
istart = 1;
iend = numel(time1); 
ntime = numel(time1);
time = time1(istart:iend);
count = 1;
for k = istart:irate:iend

    xx1 = posx1(:,:,k);
    yy1 = posy1(:,:,k);
    tt =  conc1(:,:,k);
    vec1 = [xx1(:,:);xx1(1,:)];
    vec2 = [yy1(:,:);yy1(1,:)];
    N = length(tt);
    rbn = 1 * (ones(N,1) - tt) + 0.1*tt;
    vec3 = [rbn(:,:);rbn(1,:)];
    h = cline(vec1,vec2,vec3);
    set(h,'linewidth',3)
    %h = figure(1);clf;
    %plot(vec1,vec2,'r-','linewidth',3)
    title(['t = ' num2str(time(k),'%4.2e')])
    axis equal
    pause(0.1)
    %set(h,'Units','Inches');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'xcolor','white')
    set(gca,'ycolor','white')
    %set(gca,'visible','off')
    
   
    %pos = get(h,'Position');
    %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    filename = ['relaxVideo_' num2str(count,'%04.f') '.pdf'];
  
    print(gcf, filename,'-dpdf')
    pause(0.01);clf;
    count = count + 1;
end
