addpath ..
set(0,'DefaultAxesFontSize',22)
options.savefig = false;
clf;
irate = 10; % controls the speed of the visualization

conc = [0,0.1,0.3,0.5,0.7];
beta = [0,0.1,1,10];
chi = [0,10,20];
colors =['r','b','k','g','c'];
 
for i=1:numel(chi)
    for j = 1:numel(beta)
        for k = 1:numel(conc)
          str = ['Chi' num2str(chi(i)) '_ra065_beta' ...
              num2str(beta(j)) '_conc' num2str(conc(k))];
          
          str = strrep(str,'.','p');
          file = [str '.bin'];
     
          [posx,posy,concc,ea,el,time,xvel,yvel,ten] = loadFile(file);
          
          istart = 1;
          iend = numel(time);
          ntime = iend;
            
           oc = curve;
          [~,~,L] = oc.geomProp([posx(:,1,1);posy(:,1,1)]);
          area = zeros(ntime,1);
          length = zeros(ntime,1);
          ra = zeros(ntime,1);
          incAng = zeros(ntime,1);
            
          for l = istart:1:iend
            [ra(l),area(l),length(l)] = oc.geomProp([posx(:,1,l);posy(:,1,l)]);
          end
          figure(j);
          plot(time,ra,colors(k),'linewidth',3)
          hold on
        end
        hold off
        set(gcf,'units','inch','position',[0,0,10,2*10/3]);
        legend("Conc = 0", "Conc = 0.1", "Conc = 0.3", "Conc = 0.5", "Conc = 0.7") 
        title(['Varying Concentration: Chi' num2str(chi(i)) ', ra058,' ' beta' num2str(beta(j))])
        ylabel('Reduced Area')
        xlabel('time')
        figureLabel =['Varying conc_Chi' num2str(chi(i)) '_ra058' '_beta' num2str(beta(j))];
        figureLabel = strrep(figureLabel,'.','p'); 
        saveas(gcf,[figureLabel '.jpg'])
    end
end