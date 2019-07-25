clc;clear;
m = [5];%,64,128,256,2056];
ms = string({'5'});%,'64','128','256','2056'});
beta = [1];%, 0.1, 1, 10];
betas = string({'1'});%,'pt1','1','10'});

for i = 1:length(m)
    for j = 1:length(beta)
        [Xfinal] = relaxation1Ves(beta(j), m(i), betas{j}, ms{i});
    end
end