addpath ~/projects/brinkman/vesicle_code/examples/output
addpath ~/projects/brinkman/vesicle_code/src
scale = 2;

file{1} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em5/shear1VesData.bin';
file{2} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em5/shear1VesData_Part2.bin';
file{3} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em5/shear1VesData_Part3.bin';
file{4} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em5/shear1VesData_Part4.bin';
file{5} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em5/shear1VesData_Part5.bin';
file{6} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em5/shear1VesData_Part5.bin';
file{7} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em5/shear1VesData.bin';

file{8} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em4p5/shear1VesData.bin';
file{9} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em4p5/shear1VesData.bin';
file{10} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em4p5/shear1VesData.bin';
file{11} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em4p5/shear1VesData_Part2.bin';
file{12} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em4p5/shear1VesData_Part2.bin';
file{13} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em4p5/shear1VesData_Part3.bin';
file{14} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em4p5/shear1VesData_Part3.bin';

file{15} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em4/shear1VesData.bin';
file{16} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em4/shear1VesData.bin';
file{17} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em4/shear1VesData.bin';
file{18} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em4/shear1VesData.bin';
file{19} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em4/shear1VesData_Part2.bin';
file{20} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em4/shear1VesData_Part2.bin';
file{21} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em4/shear1VesData.bin';

file{22} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em3p5/shear1VesData.bin';
file{23} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em3p5/shear1VesData.bin';
file{24} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em3p5/shear1VesData.bin';
file{25} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em3p5/shear1VesData.bin';
file{26} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em3p5/shear1VesData.bin';
file{27} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em3p5/shear1VesData.bin';
file{28} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em3p5/shear1VesData.bin';

file{29} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em3/shear1VesData.bin';
file{30} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em3/shear1VesData.bin';
file{31} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em3/shear1VesData.bin';
file{32} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em3/shear1VesData.bin';
file{33} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em3/shear1VesData.bin';
file{34} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em3/shear1VesData.bin';
file{35} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em3/shear1VesData.bin';

file{36} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em2p5/shear1VesData.bin';
file{37} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em2p5/shear1VesData.bin';
file{38} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em2p5/shear1VesData.bin';
file{39} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em2p5/shear1VesData.bin';
file{40} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em2p5/shear1VesData.bin';
file{41} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em2p5/shear1VesData.bin';
file{42} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em2p5/shear1VesData.bin';

file{43} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em2/shear1VesData.bin';
file{44} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em2/shear1VesData.bin';
file{45} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em2/shear1VesData.bin';
file{46} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em2/shear1VesData.bin';
file{47} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em2/shear1VesData.bin';
file{48} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em2/shear1VesData.bin';
file{49} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em2/shear1VesData.bin';

file{50} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em1p5/shear1VesData.bin';
file{51} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em1p5/shear1VesData.bin';
file{52} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em1p5/shear1VesData.bin';
file{53} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em1p5/shear1VesData.bin';
file{54} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em1p5/shear1VesData.bin';
file{55} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em1p5/shear1VesData.bin';
file{56} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em1p5/shear1VesData.bin';

file{57} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em1/shear1VesData.bin';
file{58} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em1/shear1VesData.bin';
file{59} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em1/shear1VesData.bin';
file{60} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em1/shear1VesData.bin';
file{61} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em1/shear1VesData.bin';
file{62} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em1/shear1VesData.bin';
file{63} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em1/shear1VesData.bin';

file{64} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0em0p5/shear1VesData.bin';
file{65} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0em0p5/shear1VesData.bin';
file{66} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0em0p5/shear1VesData.bin';
file{67} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0em0p5/shear1VesData.bin';
file{68} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0em0p5/shear1VesData.bin';
file{69} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0em0p5/shear1VesData.bin';
file{70} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0em0p5/shear1VesData.bin';

file{71} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em1_ra065_beta1p0e0/shear1VesData.bin';
file{72} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0em0p5_ra065_beta1p0e0/shear1VesData.bin';
file{73} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0_ra065_beta1p0e0/shear1VesData.bin';
file{74} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e0p5_ra065_beta1p0e0/shear1VesData.bin';
file{75} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1_ra065_beta1p0e0/shear1VesData.bin';
file{76} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e1p5_ra065_beta1p0e0/shear1VesData.bin';
file{77} = '~/projects/brinkman/vesicle_code/results/shear1Ves/Chi1p0e2_ra065_beta1p0e0/shear1VesData.bin';

xshift = [2 2 2 2 2 2 2];
xshift = [xshift xshift+6 xshift+12 xshift+18 xshift+24 ...
                xshift+30 xshift+36 xshift+42 xshift+48 ...
                xshift+54 xshift+60];

yshift = (2:4:26);
yshift = repmat(yshift,1,11);

betas = {'1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5',...
         '1e-1','1e-0.5','1e+0'};
chis = {'1e-1','1e-0.5','1e0','1e0.5','1e1','1e1.5','1e2'};

fid = fopen('data.dat','w');
%for j = 1:numel(file)
for j = 15:35
  disp(j)
  [posx,posy,~,~,~,~,~,~,n,~] = loadFile(file{j});
% some files did not print correctly to file becuase they were either
% too big or the walltime ran out.
%  if j == 11 
%    posx = posx(:,:,1:3400000);
%    posy = posy(:,:,1:3400000);
%  elseif j == 12
%    posx = posx(:,:,1:1517319);
%    posy = posy(:,:,1:1517319);
%  elseif j == 13
%    posx = posx(:,:,1:1487510);
%    posy = posy(:,:,1:1487510);
%  end
  x = [posx(:,:,end);posx(1,:,end)];
  y = [posy(:,:,end);posy(1,:,end)];
%  clear 'posx' 'posy'
  x0 = xshift(j); y0 = yshift(j);
  % x0 is a beta shift and y0 is a shear shift

  str = ['%% beta = ',betas{ceil(j/7)},...
      ',shear rate = ',chis{mod(j-1,7)+1},'\n'];
  fprintf(fid,str,'%s');
  str = '\\addplot[red,line width=0.5pt] coordinates{\n';
  fprintf(fid,str,'%s');

  for k = 1:n+1
    str = ['(' num2str(x(k)/scale+x0,'%6.4e') ',' num2str(y(k)/scale+y0,'%6.4e') ')\n']; 
    fprintf(fid,str,'%s'); 
  end
  str = '};\n\n';
  fprintf(fid,str,'%s');

  clf;
  plot(x,y,'r')
  pause(.01)
end
fclose(fid); 






