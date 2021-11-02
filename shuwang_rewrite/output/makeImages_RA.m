addpath ..
set(0,'DefaultAxesFontSize',22)
options.savefig = false;
clf;
irate = 10; % controls the speed of the visualization

% conc = [0.3,0.5];%[0,0.1,0.3,0.5,0.7];
% beta = [0,0.1,1,10];
% chi = 0;%[0,10,20];
conc = 0.3;%[0.3,0.5];%[0,0.1,0.3,0.5,0.7];
beta = [0,0.01,0.1];%[0,0.1,1,10];
chi = 5;%[1,10];%[0,10,20];

colors =['r','b','k','g'];
 

for i=1:numel(chi)
    for j = 1:numel(conc)
        for k = 1:numel(beta)
          str = ['Chi' num2str(chi(i)) '_ra065_beta' ...
              num2str(beta(k)) '_conc' num2str(conc(j))];
          
          str = strrep(str,'.','p');
          file = [str '.bin'];
     
          [posx,posy,concc,ea,el,time,xvel,yvel,ten] = loadFile(file);
          
          istart = 1;
          iend = numel(time);
          ntime = iend;
            
           oc = curve;
          [~,~,L] = oc.geomProp([posx(:,1,1);posy(:,1,1)]);
          area = zeros(ntime,1);
          %length = zeros(ntime,1);
          ra = zeros(ntime,1);
          incAng = zeros(ntime,1);
            
          for l = istart:1:iend
            [ra(l),area(l),len(l)] = oc.geomProp([posx(:,1,l);posy(:,1,l)]);
          end
          figure(j);
          plot(time,ra,colors(k),'linewidth',3)
          hold on
        end
        hold off
        set(gcf,'units','inch','position',[0,0,10,2*10/3]);
        legend("\beta = 0", "\beta = 0.01", "\beta = 0.1")
        %legend("Conc = 0", "Conc = 0.1", "Conc = 0.3", "Conc = 0.5", "Conc = 0.7") 
        %title(['Varying Beta: Chi' num2str(chi(i)) ', ra058,' ' conc' num2str(conc(j))])
        title(['Shear rate ', num2str(chi(i)), ' Initial Concentration ',num2str(conc(j))])
        %title(['Relaxation, Initial Concentration ',num2str(conc(j))])
        ylabel('Reduced Area')
        xlabel('time')
        %figureLabel =['Varying Beta_Chi' num2str(chi(i)) '_ra058' '_conc' num2str(conc(j))];
        figureLabel = ['Shear rate ', num2str(chi(i)), ' Initial Concentration ',num2str(conc(j))];
        %figureLabel = ['Relaxation Initial Concentration ',num2str(conc(j))];
        figureLabel = strrep(figureLabel,'.','p'); 
        %saveas(gcf,[figureLabel '.jpg'])
    end
 
end
