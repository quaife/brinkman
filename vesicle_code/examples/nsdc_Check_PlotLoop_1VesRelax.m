m = [32,64,128,256,2056];
ms = string({'32','64','128','256','2056'});
beta = [0, 0.1, 1, 10];
betas = string({'0','pt1','1','10'});

for i = 1:length(m)
    for j = 1:length(beta)
        file = sprintf('%s','output/relaxation1Ves',ms,'_',betas,'.log');
        
    end
end