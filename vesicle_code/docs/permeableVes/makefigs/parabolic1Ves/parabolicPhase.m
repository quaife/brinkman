addpath ~/projects/brinkman/vesicle_code/examples/output
addpath ~/projects/brinkman/vesicle_code/src
scale = 1.5;

file{1} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u0p1B1em5Data.bin';
file{2} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1em0p5B1em5Data.bin';
file{3} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1B1em5Data.bin';
file{4} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e0p5B1em5Data.bin';
file{5} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u10B1em5Data.bin';
file{6} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e1p5B1em5bData.bin';
%file{7} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em5Data.bin';
file{7} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em5aData.bin';

file{8} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u0p1B1em4Data.bin';
file{9} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1em0p5B1em4Data.bin';
file{10} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1B1em4Data.bin';
file{11} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e0p5B1em4Data.bin';
file{12} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u10B1em4Data.bin';
%file{13} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e1p5B1em4Data.bin';
file{13} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e1p5B1em4bData.bin';
file{14} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em4Data.bin';

file{15} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u0p1B1em3Data.bin';
file{16} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1em0p5B1em3Data.bin';
file{17} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1B1em3Data.bin';
file{18} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e0p5B1em3Data.bin';
file{19} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u10B1em3Data.bin';
file{20} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e1p5B1em3Data.bin';
file{21} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em3Data.bin';
file{21} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em3aData.bin';

file{22} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u0p1B1em2Data.bin';
file{23} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1em0p5B1em2Data.bin';
file{24} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1B1em2Data.bin';
file{25} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e0p5B1em2Data.bin';
file{26} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u10B1em2Data.bin';
file{27} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e1p5B1em2Data.bin';
file{28} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em2Data.bin';

file{29} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u0p1B1em1Data.bin';
file{30} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1em0p5B1em1Data.bin';
file{31} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1B1em1Data.bin';
file{32} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e0p5B1em1Data.bin';
file{33} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u10B1em1Data.bin';
file{34} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e1p5B1em1Data.bin';
file{35} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em1Data.bin';

file{36} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u0p1B1em0Data.bin';
file{37} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1em0p5B1em0Data.bin';
file{38} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1B1em0Data.bin';
file{39} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e0p5B1em0Data.bin';
file{40} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u10B1em0Data.bin';
file{41} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u1e1p5B1em0Data.bin';
file{42} = '~/projects/brinkman/vesicle_code/results/parabolic/pflowR10u100B1em0Data.bin';

xshift = [2 2 2 2 2 2 2];
xshift = [xshift xshift+6 xshift+12 xshift+18 xshift+24 xshift+30];

yshift = (2:4:26);
yshift = repmat(yshift,1,6);

betas = {'1e-5','1e-4','1e-3','1e-2','1e-1','1e+0'};
chis = {'1e-1','1e-0.5','1e0','1e0.5','1e1','1e1.5','1e2'};

fid = fopen('data.dat','w');
%for j = 1:numel(file)
for j = 7:7
%  if (j ~= 14 && j~=21 && j ~= 35 && j ~= 42)
    disp(j)
    [posx,posy,~,~,~,~,~,~,n,~] = loadFile(file{j});
    x = [posx(:,:,end);posx(1,:,end)];
    y = [posy(:,:,end);posy(1,:,end)];
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
%  end
end
fclose(fid); 





