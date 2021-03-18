conc = [0.3];%[0,0.1,0.3,0.5,0.7];
beta = 0.1;%[0,0.01, 0.1];%[0,0.1,1,10];
chi = 5;%[1,10];%[0,10,20];

for i=1:numel(chi)
    for j = 1:numel(conc)
        for k = 1:numel(beta)
          str = ['Chi' num2str(chi(i)) '_ra065_beta' ...
              num2str(beta(k)) '_conc' num2str(conc(j))];
          str = strrep(str,'.','p');
          shear1Ves(beta(k),chi(i),conc(j),str);
        end
    end
end


