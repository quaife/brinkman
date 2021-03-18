conc = [0,0.3];       %[0,0.1,0.3,0.5,0.7];
beta = 0;       %[0,0.01, 0.1];%[0,0.1,1,10];
chi  = [600];
shortaxis = [0.324, 0.405, 0.51, 0.691];
RA =        [0.6, 0.75, 0.85, 0.95];
% L        [4.4350 4.6141 4.8695 5.3569]
scaleL = [1.2079, 1.1610, 1.1001, 1];

for l =1:numel(shortaxis)
    for i=1:numel(chi)
        for j = 1:numel(conc)
            for k = 1:numel(beta)
              str = ['Parabolic_RA' num2str(RA(l)) '_Conc' num2str(conc(j)) '_Chi' num2str(chi(i)) '_beta' ...
                  num2str(beta(k))];
              str = strrep(str,'.','p');
              shear1Ves(beta(k),chi(i),conc(j),shortaxis(l), scaleL(l),str);
              end
          end
    end
end

