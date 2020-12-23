addpath ..
set(0,'DefaultAxesFontSize',22)
options.savefig = false;
clf;
irate = 10; % controls the speed of the visualization
% 
% conc = [0];%,0.1,0.3,0.5,0.7];
% beta = [0,0.1,1,10];
% chi = [0];%,10,20];
% colors =['r','b','k','g'];
% 
% for i=1:numel(chi)
%     for j = 1:numel(conc)
%         for k = 1:numel(beta)
%           str = ['Chi' num2str(chi(i)) '_ra065_beta' ...
%               num2str(beta(k)) '_conc' num2str(conc(j))];
%           
%           str = strrep(str,'.','p');
%           str = strrep(str,'-0','m');
%           str = strrep(str,'+0','');
%           file = [str '.bin'];
%      
%           [posx,posy,concc,ea,el,time,xvel,yvel] = loadFile(file);
%           
%           istart = 1;
%           iend = numel(time);
%           ntime = iend;
%             
           oc = curve;
%           [~,~,L] = oc.geomProp([posx(:,1,1);posy(:,1,1)]);
%           area = zeros(ntime,1);
%           length = zeros(ntime,1);
%           ra = zeros(ntime,1);
%           incAng = zeros(ntime,1);
%             
%           for l = istart:1:iend
%             [ra(l),area(l),length(l)] = oc.geomProp([posx(:,1,l);posy(:,1,l)]);
%           end
%           figure(j);
%           plot(time,ra,colors(k),'linewidth',3)
%           hold on
%         end
%         hold off
%         legend("\beta = 0", "\beta = 0.1", "\beta = 1", "\beta = 10")
%         %legend("Conc = 0", "Conc = 0.1", "Conc = 0.3", "Conc = 0.5", "Conc = 0.7") 
%         title(['Relaxation'])
%         ylabel('Reduced Area')
%         xlabel('time')
%     end
% end


% 
% 
 file1 = 'Chi80_ra08146_beta01_conc03_dt1e3.bin';
%  file2 = 'Chi10_ra085_beta0p1_conc0.bin';
%  file3 = 'Chi10_ra099_beta0p1_conc0.bin';
% % 
 [posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1] = loadFile(file1);
  
%   [posx2,posy2,conc2,ea2,el2,time2,xvel2,yvel2] = loadFile(file2);
%   [posx3,posy3,conc3,ea3,el3,time3,xvel3,yvel3] = loadFile(file3);
% % % % % load positions, curvature, errors, time, and velocities

%  istart = 1;
%  iend = numel(time1);
%  ntime = iend;
%  ntime = numel(time1);
% % 
% % % %time = time(istart:iend);
% % % % % 
% % % % oc = curve;
% % % % [~,~,L] = oc.geomProp([posx1(:,1,1);posy1(:,1,1)]);
% % % % area = zeros(ntime,1);
% % % % length = zeros(ntime,1);
% % % % ra = zeros(ntime,1);
% % % % incAng = zeros(ntime,1);
% % % % 
%  %for k = istart:1:iend
% % %  [ra1(k),area1(k),length1(k)] = oc.geomProp([posx1(:,1,k);posy1(:,1,k)]);
% %   %[ra2(k),area2(k),length2(k)] = oc.geomProp([posx2(:,1,k);posy2(:,1,k)]);
% %  %end
% % %
% % %  
% % %   istart = 1;
% % %   iend = numel(time3);
% % %   ntime = iend;
% % %  for k = istart:1:iend
% % %  [ra3(k),area3(k),length3(k)] = oc.geomProp([posx3(:,1,k);posy3(:,1,k)]);
% % % end
% % %[length1(1), length2(1), length3(1)]
% % %figure(1);clf;
% % % plot(time1,ra1, 'r', 'linewidth', 3)
% % % hold on
% % % plot(time2,ra2, 'b', 'linewidth', 3)
% % % plot(time3,ra3, 'g', 'linewidth', 3)
% % % hold off
% % %  title(['Shear Flow: \gamma = 10'])
% % %  ylabel('Reduced Area')
% % %  xlabel('time')
% % %  legend('RA_0 = 0.65','RA_0 = 0.85','RA_0 = 0.99')
% % %   
% % % count = 1;
% for k = istart:irate:iend
% % % % % %  xx = interpft(posx(:,:,k),256); yy = interpft(posy(:,:,k),256);  
%       xx1 = posx1(:,:,k);
%       yy1 = posy1(:,:,k);
% % % % % %   xx2 = posx2(:,:,k);
% % % % % %   yy2 = posy2(:,:,k);
%      tt = conc1(:,:,k);
%       vec1 = [xx1(:,:);xx1(1,:)];
%       vec2 = [yy1(:,:);yy1(1,:)];
%       vec3 = [tt(:,:);tt(1,:)];
% % % % % %   vec3 = [xx2(:,:);xx2(1,:)];
% % % % % %   vec4 = [yy2(:,:);yy2(1,:)];
%     clf;
%      h = cline(vec1,vec2,vec3);
%      set(h,'linewidth',3)
%      colorbar
% % % %    h = figure(2);clf;
% % % %    plot(vec1,vec2,'r-','linewidth',3)
% % % % % % %   hold on   
% % % % % %   plot(vec3,vec4,'r-','linewidth',3)
% % % % % %  hold off
%     axis equal
%     pause(0.1)
% % %   %set(h,'Units','Inches');
% % %   set(gca,'xtick',[])
% % %   set(gca,'ytick',[])
% % %   set(gca,'xcolor','white')
% % %   set(gca,'ycolor','white')
% % %   titleStr = ['t = ' num2str(time1(k),'%4.2e')];% ...
% % %     % ' eA = ' num2str(ea1(k),'%4.2e') ...
% % %     % ' eL = ' num2str(el1(k),'%4.2e')];
% % %   title(titleStr);
% % %   
% % %   %pos = get(h,'Position');
% % %   %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % %   filename = ['relaxVideo_' num2str(count,'%04.f') '.pdf'];
% % % %   
% % % %   print(gcf, filename,'-dpdf')
% % %   pause(0.01);clf;
% % % %   count = count + 1;
%     end
% % % % % 
