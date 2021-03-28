conc = [0,0.3];
beta = 0; 
chi  = [800];
scaleL = [0.538, 0.635, 0.714, 0.827];
RA =     [0.6, 0.75, 0.85, 0.95];
shortaxis = [3.45, 2.48, 1.96, 1.45];

for l =1:numel(shortaxis)
    for i=1:numel(chi)
        for j = 1:numel(conc)
            for k = 1:numel(beta)
              str = ['Parabolic_RA' num2str(RA(l)) '_Conc' num2str(conc(j)) '_Chi' num2str(chi(i)) '_beta' ...
                  num2str(beta(k))];
              str = strrep(str,'.','p');
              parabolic1Ves(beta(k),chi(i),conc(j),shortaxis(l), scaleL(l),str);
              end
          end
    end
end

