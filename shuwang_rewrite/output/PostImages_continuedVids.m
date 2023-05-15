clc;clear
addpath ..
%addpath 'D:\reserach\Runsforpaper\stenosis'
ax = [-1 1 -1 1];

file1 = 'Chi2p5_shax4p37_scL0p397_Conc0p5_Beta0_n1024_nbd1024_dt5en6_bmax1_bmin0p1_eps0p04_a100_contracting_rand.bin';
file2 = 'Chi2p5_shax4p37_scL0p397_Conc0p5_Beta0_n1024_nbd1024_dt5en6_bmax1_bmin0p1_eps0p04_a100_contracting_rand_contFromMid.bin';
% Merge files
[posx,posy,conc,ea,el,time,xvel,yvel,ten] = loadFile(file1);
[posx2,posy2,conc2,ea2,el2,time2,xvel2,yvel2,ten2] = loadFile(file2);

% Runs that restarted from t - 5000
posx = posx(:,:,1:end-5001); 
posy = posy(:,:,1:end-5001);
conc = conc(:,:,1:end-5001);
ea = ea(1:end-5001);
el = el(1:end-5001);
time = time(1:end-5001);
xvel = xvel(:,:,1:end-5001);
yvel = yvel(:,:,1:end-5001);
ten = ten(:,:,1:end-5001);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Runs that restarted bring back one
%k = 618;%for whatever reason, this reatrt is missing timesteps
posx2 = posx2(:,:,2:end); 
posy2 = posy2(:,:,2:end);
conc2 = conc2(:,:,2:end);
ea2 = ea2(2:end);
el2 = el2(2:end);
time2 = time2(2:end);
xvel2 = xvel2(:,:,2:end);
yvel2 = yvel2(:,:,2:end);
ten2(:,:,2) = ten(:,:,end);
ten2 = ten2(:,:,2:end);
% posx2 = posx2(:,:,2:end); 
% posy2 = posy2(:,:,2:end);
% conc2 = conc2(:,:,2:end);
% ea2 = ea2(2:end);
% el2 = el2(2:end);
% time2 = time2(2:end);
% xvel2 = xvel2(:,:,2:end);
% yvel2 = yvel2(:,:,2:end);
% ten2 = ten2(:,:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bmin = 0.1;
bmax = 1;

[sx11,sx12,sx13] = size(posx);
[sx21,sx22,sx23] = size(posx2);
posx(1:sx11, sx12, 1:sx13) = posx;
posx(1:sx11, sx12, sx13+1:sx23+sx13) = posx2;

[sy11,sy12,sy13] = size(posy);
[sy21,sy22,sy23] = size(posy2);
posy(1:sy11, sy12, 1:sy13) = posy;
posy(1:sy11, sy12, sy13+1:sy23+sy13) = posy2;

[sc11,sc12,sc13] = size(conc);
[sc21,sc22,sc23] = size(conc2);
conc(1:sc11, sc12, 1:sc13) = conc;
conc(1:sc11, sc12, sc13+1:sc23+sc13) = conc2;

[st11,st12,st13] = size(ten);
[st21,st22,st23] = size(ten2);
ten(1:st11, st12, 1:st13) = ten;
ten(1:st11, st12, st13+1:st23+st13) = ten2;

[sy11,sy12,sy13] = size(yvel);
[sy21,sy22,sy23] = size(yvel2);
yvel(1:sy11, sy12, 1:sy13) = yvel;
yvel(1:sy11, sy12, sy13+1:sy23+sy13) = yvel2;
[sy11,sy12,sy13] = size(xvel);
[sy21,sy22,sy23] = size(xvel2);
xvel(1:sy11, sy12, 1:sy13) = xvel;
xvel(1:sy11, sy12, sy13+1:sy23+sy13) = xvel2;

time2 = time2 + time(end);
time = [time;time2];
ea = [ea;ea2+ea(end)];
el = [el;el2+el(end)];

posx1 = posx; posy1 = posy; conc1 = conc; ten1 = ten; yvel1 = yvel; xvel1 = xvel; time1 = time; ea1 = ea; el1 = el;
pause

file3 = 'Chi2p5_shax4p37_scL0p397_Conc0_Beta0_n1024_nbd1024_dt5en6_bmax1_bmin0p1_eps0p04_a100_contracting_contFromEnd_contFromEnd.bin';
[posx3,posy3,conc3,ea3,el3,time3,xvel3,yvel3,ten3] = loadFile(file3);

%Runs that restarted bring back one
%in case the vesicle collides with the wall
k = 24000; %length(posx3);
posx3 = posx3(:,:,2:k); 
posy3 = posy3(:,:,2:k);
conc3 = conc3(:,:,2:k);
ea3 = ea3(2:k);
el3 = el3(2:k);
time3 = time3(2:k);
xvel3 = xvel3(:,:,2:k);
yvel3 = yvel3(:,:,2:k);
ten3(:,:,2) = ten2(:,:,end);
ten3 = ten3(:,:,2:k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[sx11,sx12,sx13] = size(posx);
[sx21,sx22,sx23] = size(posx3);
posx4(1:sx11, sx12, 1:sx13) = posx;
posx4(1:sx11, sx12, sx13+1:sx23+sx13) = posx3;

[sy11,sy12,sy13] = size(posy);
[sy21,sy22,sy23] = size(posy3);
posy4(1:sy11, sy12, 1:sy13) = posy;
posy4(1:sy11, sy12, sy13+1:sy23+sy13) = posy3;

[sc11,sc12,sc13] = size(conc);
[sc21,sc22,sc23] = size(conc3);
conc4(1:sc11, sc12, 1:sc13) = conc;
conc4(1:sc11, sc12, sc13+1:sc23+sc13) = conc3;

[st11,st12,st13] = size(ten);
[st21,st22,st23] = size(ten3);
ten4(1:st11, st12, 1:st13) = ten;
ten4(1:st11, st12, st13+1:st23+st13) = ten3;

[sy11,sy12,sy13] = size(yvel);
[sy21,sy22,sy23] = size(yvel3);
yvel4(1:sy11, sy12, 1:sy13) = yvel;
yvel4(1:sy11, sy12, sy13+1:sy23+sy13) = yvel3;
[sy11,sy12,sy13] = size(xvel);
[sy21,sy22,sy23] = size(xvel3);
xvel4(1:sy11, sy12, 1:sy13) = xvel;
xvel4(1:sy11, sy12, sy13+1:sy23+sy13) = xvel3;

time3 = time3 + time(end);
time4 = [time;time3];
ea4 = [ea;ea3+ea(end)];
el4 = [el;el3+el(end)];

posx1 = posx4; posy1 = posy4; conc1 = conc4; ea1 = ea4; el1 = el4; time1 = time4;
xvel1 = xvel4; yvel1 = yvel4; ten1 = ten4;

% file4 = 'BendMod0p55_Chip25_shax5p7_scL0p311_Conc0_Beta0_n1024_nbd1024_dt1en4_bmax0p55_bmin0p3025_eps0p04_a100_longchoke_rand_contFromEnd_contFromEnd.bin';
% [posx5,posy5,conc5,ea5,el5,time5,xvel5,yvel5,ten5] = loadFile(file4);
% 
% %Runs that restarted bring back one
% posx5 = posx5(:,:,2:end); 
% posy5 = posy5(:,:,2:end);
% conc5 = conc5(:,:,2:end);
% ea5 = ea5(2:end);
% el5 = el5(2:end);
% time5 = time5(2:end);
% xvel5 = xvel5(:,:,2:end);
% yvel5 = yvel5(:,:,2:end);
% ten5 = ten5(:,:,2:end);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [sx11,sx12,sx13] = size(posx4);
% [sx21,sx22,sx23] = size(posx5);
% posx(1:sx11, sx12, 1:sx13) = posx4;
% posx(1:sx11, sx12, sx13+1:sx23+sx13) = posx5;
% 
% [sy11,sy12,sy13] = size(posy4);
% [sy21,sy22,sy23] = size(posy5);
% posy(1:sy11, sy12, 1:sy13) = posy4;
% posy(1:sy11, sy12, sy13+1:sy23+sy13) = posy5;
% 
% [sc11,sc12,sc13] = size(conc4);
% [sc21,sc22,sc23] = size(conc5);
% conc(1:sc11, sc12, 1:sc13) = conc4;
% conc(1:sc11, sc12, sc13+1:sc23+sc13) = conc5;
% 
% [st11,st12,st13] = size(ten4);
% [st21,st22,st23] = size(ten5);
% ten(1:st11, st12, 1:st13) = ten4;
% ten(1:st11, st12, st13+1:st23+st13) = ten5;
% 
% [sy11,sy12,sy13] = size(yvel4);
% [sy21,sy22,sy23] = size(yvel5);
% yvel(1:sy11, sy12, 1:sy13) = yvel4;
% yvel(1:sy11, sy12, sy13+1:sy23+sy13) = yvel5;
% [sy11,sy12,sy13] = size(xvel4);
% [sy21,sy22,sy23] = size(xvel5);
% xvel(1:sy11, sy12, 1:sy13) = xvel4;
% xvel(1:sy11, sy12, sy13+1:sy23+sy13) = xvel5;
% 
% time5 = time5 + time4(end);
% time = [time4;time5];
% ea = [ea4;ea5+ea4(end)];
% el = [el4;el5+el4(end)];
% 
% file5 = 'BendMod0p55_Chip25_shax5p7_scL0p311_Conc0_Beta0_n1024_nbd1024_dt1en4_bmax0p55_bmin0p3025_eps0p04_a100_longchoke_rand_contFromEnd_contFromEnd_contFromEnd.bin';
% [posx6,posy6,conc6,ea6,el6,time6,xvel6,yvel6,ten6] = loadFile(file5);
% 
% %Runs that restarted bring back one
% posx6 = posx6(:,:,2:end); 
% posy6 = posy6(:,:,2:end);
% conc6 = conc6(:,:,2:end);
% ea6 = ea6(2:end);
% el6 = el6(2:end);
% time6 = time6(2:end);
% xvel6 = xvel6(:,:,2:end);
% yvel6 = yvel6(:,:,2:end);
% ten6 = ten6(:,:,2:end);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [sx11,sx12,sx13] = size(posx);
% [sx21,sx22,sx23] = size(posx6);
% posx1(1:sx11, sx12, 1:sx13) = posx;
% posx1(1:sx11, sx12, sx13+1:sx23+sx13) = posx6;
% 
% [sy11,sy12,sy13] = size(posy);
% [sy21,sy22,sy23] = size(posy6);
% posy1(1:sy11, sy12, 1:sy13) = posy;
% posy1(1:sy11, sy12, sy13+1:sy23+sy13) = posy6;
% 
% [sc11,sc12,sc13] = size(conc);
% [sc21,sc22,sc23] = size(conc6);
% conc1(1:sc11, sc12, 1:sc13) = conc;
% conc1(1:sc11, sc12, sc13+1:sc23+sc13) = conc6;
% 
% [st11,st12,st13] = size(ten);
% [st21,st22,st23] = size(ten6);
% ten1(1:st11, st12, 1:st13) = ten;
% ten1(1:st11, st12, st13+1:st23+st13) = ten6;
% 
% [sy11,sy12,sy13] = size(yvel);
% [sy21,sy22,sy23] = size(yvel6);
% yvel1(1:sy11, sy12, 1:sy13) = yvel;
% yvel1(1:sy11, sy12, sy13+1:sy23+sy13) = yvel6;
% [sy11,sy12,sy13] = size(xvel);
% [sy21,sy22,sy23] = size(xvel6);
% xvel1(1:sy11, sy12, 1:sy13) = xvel;
% xvel1(1:sy11, sy12, sy13+1:sy23+sy13) = xvel6;
% 
% time6 = time6 + time(end);
% time1 = [time;time6];
% ea1 = [ea;ea6+ea(end)];
% el1 = [el;el6+el(end)];
% 

% posx1 = posx4; posy1 = posy4; conc1 = conc4; ea1 = ea4; el1 = el4; time1 = time4;
% xvel1 = xvel4; yvel1 = yvel4; ten1 = ten4;

% pause
% %generate walls
% if 1
%     Nbd = 1024;
%     geomCenter = [0;0];
%     wallGeometry = 'longchoke';
%     oc = curve(Nbd);
%     [~,Xwalls] = oc.initConfig(Nbd,false,...
%                  'scale', 1, ...
%                  'center', geomCenter, 'geometry', wallGeometry);
%     [xwalls,ywalls] = oc.getXY(Xwalls);
% else
%     xwalls = 0;
%     ywalls = 0;
% end
% 
% irate = 500; 
% istart = 1;%floor(numel(time)/2);
% iend = numel(time1);
% 
% N = length(posx1(:,:,1));
% oc = curve(N);
% [~, ~, Len] = oc.geomProp([posx1(:,:,1);posy1(:,:,1)]);
% a = 100;
% eps = 0.04;
% 
% i = 1;
% 
% for k = istart:irate:iend
%   figure(2)
%   clf; 
%        
%   xx1 = posx1(:,:,k);
%   yy1 = posy1(:,:,k);
%   tt1 = conc1(:,:,k);
% 
%   xxx1 = [xx1(:,:);xx1(1,:)];
%   yyy1 = [yy1(:,:);yy1(1,:)];
%   ttt1 = [tt1(:,:);tt1(1,:)];
% 
%   [~,area(i),~] = oc.geomProp([xx1;yy1]);
% 
%   % ----- plot title 
%   titleStr = ['t = ' num2str(time1(k),'%4.2e') ...
%       ' eA = ' num2str(ea1(k),'%4.2e') ...
%       ' eL = ' num2str(el1(k),'%4.2e')];
%   title(titleStr);
% 
% 
%   % ----- vesicle       
%   N = numel(conc1(:,:,k));
%   %rbn = 1 * (ones(N,1) - conc1(:,:,k)) + 0.001*conc1(:,:,k);
%   rbn = (bmin - bmax)/2*tanh(3*(tt1 - 0.5)) + (bmax + bmin)/2;
%   rbn1 = [rbn(:,:);rbn(1,:)];
% 
%   fu = 0.25*conc1(:,:,k).^2.*(1-conc1(:,:,k)).^2;
%   du = oc.diffFT(conc1(:,:,k))/Len;
%   tens = -ten1(:,:,k) - a/eps*(fu-(eps^2/2)*du.^2);
%   tension = [tens;tens(1)];
%   h = cline(xxx1,yyy1,rbn1);
%   hold on
%   plot(xxx1(1,:),yyy1(1,:),'k.','markersize',30);
%   plot(xxx1,yyy1,'k')
%   set(h,'linewidth',3)
%   q = colorbar;
%   set(q,'Limits',[0 1]) 
%   % ----- walls
%    plot(xwalls,ywalls,'k','linewidth',2)
%     
%   % ----- axis properties
%   axis equal
%   axis([min(xxx1)-1, max(xxx1)+1,min(yyy1)-1, max(yyy1)+1 ])
%   set(gca,'xtick',[])
%   set(gca,'ytick',[])
%   set(gca,'xcolor','white')
%   set(gca,'ycolor','white')
%   hold off
% 
%   pause(0.1);
%   i = i +1;
% end
%  
