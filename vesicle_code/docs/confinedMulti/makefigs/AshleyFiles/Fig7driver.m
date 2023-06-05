figure(1); clf;
ha = tight_subplot(3,1,0.02,0,[0 0.02]);
MakeSnaps('Stenosis_RAp4_SC.mat',1,false,ha);
MakeSnaps('Stenosis_RAp4_SCp55.mat',2,false,ha);
MakeSnaps('Stenosis_RAp4_MCp5.mat',3,true,ha);

f = figure(1);
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)])
set(f,'Position',[0 0 19 6.4]);
print(f,'-dpdf','STENOSIS_RAp4MCp5.pdf','-r0')
