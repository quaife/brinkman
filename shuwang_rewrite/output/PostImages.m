
addpath ..
ax = [-50 50 -4 4];

 file1 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt5en5_bmax1_bmin0p1_eps0p04_a100_contracting_left.bin';
 file2 = 'Chip25_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt5en5_bmax1_bmin0p1_eps0p04_a100_contracting_right.bin';
% %file1 = 'Chi400_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en6_bmax1_bmin0p1_eps0p04_a100_longchoke.bin';
% 
 bmin = 0.1;
[posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);
[posx2,posy2,conc2,ea2,el2,time2,xvel2,yvel2,ten2] = loadFile(file2);

% posx2 = posx2(:,:,1:5:end);
% posy2 = posy2(:,:,1:5:end);
% conc2 = conc2(:,:,1:5:end);
% time2 = time2(1:5:end);
if 1
    Nbd = 1024;
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
 iend = numel(time1)-1;
 ntime = numel(time1);
 N = length(posx1(:,:,1));
 %time1 = time1(istart:iend);
 oc = curve(N);
 [~, ~, Len] = oc.geomProp([posx1(:,:,1);posy1(:,:,1)]);
 a = 100;
 eps = 0.04;
 i = 1;
  for k = istart:irate:iend
       figure(2)
       clf; 
       
       xx1 = posx1(:,:,k);
       yy1 = posy1(:,:,k);
       tt1 = conc1(:,:,k);
       
       xxx1 = [xx1(:,:);xx1(1,:)];
       yyy1 = [yy1(:,:);yy1(1,:)];
       ttt1 = [tt1(:,:);tt1(1,:)];
       
       [~,area(i),~] = oc.geomProp([xx1;yy1]);
%        
       xx2 = posx2(:,:,k);
       yy2 = posy2(:,:,k);
       tt2 = conc2(:,:,k);
       
       xxx2 = [xx2(:,:);xx2(1,:)];
       yyy2 = [yy2(:,:);yy2(1,:)];
       ttt2 = [tt2(:,:);tt2(1,:)];

       % ----- plot title 
       titleStr = ['t = ' num2str(time1(k),'%4.2e') ...
           ' eA = ' num2str(ea1(k),'%4.2e') ...
           ' eL = ' num2str(el1(k),'%4.2e')];
       title(titleStr);


       % ----- vesicle       
       N = numel(conc1(:,:,k));
       %rbn = 1 * (ones(N,1) - conc1(:,:,k)) + 0.001*conc1(:,:,k);
       rbn = (bmin - 1)/2*tanh(3*(tt1 - 0.5)) + (1 + bmin)/2;
       rbn1 = [rbn(:,:);rbn(1,:)];
       
       fu = 0.25*conc1(:,:,k).^2.*(1-conc1(:,:,k)).^2;
       du = oc.diffFT(conc1(:,:,k))/Len;
       tens = -ten1(:,:,k) - a/eps*(fu-(eps^2/2)*du.^2);
       tens1 = [tens;tens(1)];
       
       h = cline(xxx1,yyy1,rbn1);
       hold on
       plot(xxx1(1,:),yyy1(1,:),'k.','markersize',30);
       % second curve
       hold on
       rbn2 = (bmin - 1)/2*tanh(3*(tt2 - 0.5)) + (1 + bmin)/2;
       rbn2 = [rbn2(:,:);rbn2(1,:)];
       
       fu = 0.25*conc2(:,:,k).^2.*(1-conc2(:,:,k)).^2;
       du = oc.diffFT(conc2(:,:,k))/Len;
       tens2 = -ten2(:,:,k) - a/eps*(fu-(eps^2/2)*du.^2);
       tens2 = [tens2;tens2(1)];
       
       h = cline(xxx2,yyy2,rbn2);
      plot(xxx1,yyy1,'k')
       set(h,'linewidth',3)
       colorbar
       hold on
       plot(xxx1(1,:),yyy1(1,:),'k.','markersize',20)
       plot(xxx2(1,:),yyy2(1,:),'k.','markersize',20)
%        % ----- walls
        plot(xwalls,ywalls,'k','linewidth',2)
%     
%        % ----- axis properties
       axis equal
       axis([min(xxx1)-1, max(xxx1)+1,min(yyy1)-1, max(yyy1)+1 ])
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       set(gca,'xcolor','white')
       set(gca,'ycolor','white')
       hold off
       
       % ----- legend
       %legend('', 'No Concentration','','0.3 floppy','','0.5 floppy', '')

       pause(0.1);
       i = i +1;
  end
 
