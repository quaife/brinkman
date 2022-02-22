addpath ..

ax = [-50 50 -4 4];

file1 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_right.bin';
[posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);
file2 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_left.bin';
[posx2,posy2,conc2,ea2,el2,time2,xvel2,yvel2,ten2] = loadFile(file2);
file3 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_right.bin';
[posx3,posy3,conc3,ea3,el3,time3,xvel3,yvel3,ten3] = loadFile(file3);
file4 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_right.bin';
[posx4,posy4,conc4,ea4,el4,time4,xvel4,yvel4,ten4] = loadFile(file4);

bmin = 0.1;
if 1
    Nbd = 768;
    geomCenter = [0;0];
    wallGeometry = 'contracting';
    oc = curve(Nbd);
    [~,Xwalls] = oc.initConfig(Nbd,false,...
                 'scale', 1, ...
                 'center', geomCenter, 'geometry', wallGeometry);
    [xwalls,ywalls] = oc.getXY(Xwalls);
else 
    xwalls = 0;
    ywalls = 0;
end
 irate = 100; 
 istart = 1;
 iend = min([numel(time1),numel(time2),numel(time3),numel(time4)]);
 ntime = iend;
 ntime = numel(time4);
 N = length(posx1(:,:,1));
 %time1 = time1(istart:iend);
 oc = curve(N);
 
  for k = istart:irate:iend
       clf; 
       % ----- case 1
       xx1 = posx1(:,:,k);
       yy1 = posy1(:,:,k);
       tt1 = conc1(:,:,k);
       
       xxx1 = [xx1(:,:);xx1(1,:)];
       yyy1 = [yy1(:,:);yy1(1,:)];
       ttt1 = [tt1(:,:);tt1(1,:)];

       % ----- case 2       
       xx2 = posx2(:,:,k);
       yy2 = posy2(:,:,k);
       tt2 = conc2(:,:,k);
       
       xxx2 = [xx2(:,:);xx2(1,:)];
       yyy2 = [yy2(:,:);yy2(1,:)];
       ttt2 = [tt2(:,:);tt2(1,:)];

       % ----- case 3       
       xx3 = posx3(:,:,k);
       yy3 = posy3(:,:,k);
       tt3 = conc3(:,:,k);
       
       xxx3 = [xx3(:,:);xx3(1,:)];
       yyy3 = [yy3(:,:);yy3(1,:)];
       ttt3 = [tt3(:,:);tt3(1,:)];
       
       % ----- case 4       
       xx4 = posx4(:,:,k);
       yy4 = posy4(:,:,k);
       tt4 = conc4(:,:,k);
       
       xxx4 = [xx4(:,:);xx4(1,:)];
       yyy4 = [yy4(:,:);yy4(1,:)];
       ttt4 = [tt4(:,:);tt4(1,:)];
       
       
       % ----- case 1       
       N = numel(conc1(:,:,k));
       rbn = (bmin - 1)/2*tanh(3*(tt1 - 0.5)) + (1 + bmin)/2;
       rbn1 = [rbn(:,:);rbn(1,:)];
       subplot(3,2,1)
       h = cline(xxx1,yyy1,rbn1);
       hold on
       plot(xwalls,ywalls,'k','linewidth',2);
       axis equal
       axis([min(xxx1)-1, max(xxx1)+1,min(yyy1)-1, max(yyy1)+1 ])
       title("init conc rand")
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       hold off
       
%        % ----- case 2       
       N = numel(conc2(:,:,k));
       rbn = (bmin - 1)/2*tanh(3*(tt2 - 0.5)) + (1 + bmin)/2;
       rbn2 = [rbn(:,:);rbn(1,:)];
       subplot(3,2,2)
       h = cline(xxx2,yyy2,rbn2);
       title("init conc right")
       hold on
       plot(xwalls,ywalls,'k','linewidth',2);
       axis equal
       axis([min(xxx2)-1, max(xxx2)+1,min(yyy2)-1, max(yyy2)+1 ])
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       hold off
%        
%        % ----- case 3       
       N = numel(conc3(:,:,k));
       rbn = (bmin - 1)/2*tanh(3*(tt3 - 0.5)) + (1 + bmin)/2;
       rbn3 = [rbn(:,:);rbn(1,:)];
       subplot(3,2,3)
       h = cline(xxx3,yyy3,rbn3);
       hold on
       plot(xwalls,ywalls,'k','linewidth',2);
       axis equal
       axis([min(xxx3)-1, max(xxx3)+1,min(yyy3)-1, max(yyy3)+1 ])
       title("init conc left")
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       hold off
% %        set(h,'linewidth',1)
% %        plot(xxx3,yyy3,'g')
% %        plot(xxx3(1,:),yyy3(1,:),'g.','markersize',20)
% 
%        % ----- case 4       
       N = numel(conc4(:,:,k));
%        rbn = 1 * (ones(N,1) - conc4(:,:,k)) + 0.01*conc4(:,:,k);
%        rbn4 = [rbn(:,:);rbn(1,:)];
       rbn = (bmin - 1)/2*tanh(3*(tt4 - 0.5)) + (1 + bmin)/2;
       rbn4 = [rbn(:,:);rbn(1,:)];
       subplot(3,2,4)
       h = cline(xxx4,yyy4,rbn4);
       hold on
       plot(xwalls,ywalls);
       axis equal
       axis([min(xxx4)-1, max(xxx4)+1,min(yyy4)-1, max(yyy4)+1 ])
%        set(h,'linewidth',1)
%        plot(xxx4,yyy4,'m')
%        plot(xxx4(1,:),yyy4(1,:),'m.','markersize',20)
% %        
% %        % ----- walls
% %        plot(xwalls,ywalls,'k','linewidth',2)
%     
%        % ----- axis properties
%      
%        %axis equal
%        %axis([min(xxx1)-1, max(xxx1)+1,min(yyy2)-1, max(yyy2)+1 ])
%        title("No conc ")
%        set(gca,'xtick',[])
%        set(gca,'ytick',[])
%        set(gca,'xcolor','white')
%        set(gca,'ycolor','white')
%        hold off
%        
%        % ----- legend
%        %legend('', 'No Concentration','','0.3 floppy','','0.5 floppy', '')
%        % legend('', 'a = 50','','a = 100','','a = 200', '', 'a = 1000', '')
%        
       subplot(3,2,[5 6])
       plot(xxx1,yyy1);
       hold on
       plot(xxx2,yyy2);
       plot(xxx3,yyy3);
       plot(xxx4,yyy4);
       plot(xwalls,ywalls,'k','linewidth',2);
       axis equal
       axis([min(xxx1)-1, max(xxx1)+1,min(yyy1)-1, max(yyy1)+1 ])
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       hold off
       legend('random','right','left','no conc')
       pause(0.1);
  end
 
