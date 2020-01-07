function mov=movie2o(x,y,t,k,rconn)
 
m = length(x(1,:)); % number of time steps
n = length(x(:,1)); % number of spatial points
if nargin<5,
  for i=1:k:m
    hold off;
    figure(1);clf;hold on;
    plot(x(:,i),y(:,i),'.','Markersize',10);
    plot(x(1,i),y(1,i),'ro');

    axis image;
    axis([-1.1 +1.1 -1.1 +1.1]);
    set(gca,'FontSize',30,'FontName','Times') 
    set(gca,'XTick',[0,2*pi]);
    set(gca,'yTick',[-2,0,2]);
    title(strcat('t=',num2str(t(i))),'FontSize',30,'FontName','Times');
    axis on;
    mov((i-1)/k+1)=getframe(gcf);
  end
else 
  for i=1:k:m
    figure(1); clf; hold on
    for iplot=1:n
      if rconn(iplot,i)>1
        rconn(iplot,i)=1;
      elseif rconn(iplot,i)<0
        rconn(iplot,i)=0;
      end
      plot(x(iplot,i),y(iplot,i),'.','markersize',6,...
          'Color',[rconn(iplot,i),0,1-rconn(iplot,i)]); 
    end  
    plot(x(1,i),y(1,i),'gs');
    axis image;
    axis([-1.1 1.1 -1.1 1.1]);
    set(gca,'FontSize',30,'FontName','Times') 
    title(strcat('t=',num2str(t(i))),'FontSize',30,'FontName','Times');
    mov((i-1)/k+1)=getframe(gcf);
    axis on;
  end   

end
