% [beta,chi] = meshgrid(0.1:0.1:2,3.2:.2:5);
% beta = beta'; chi = chi';
% beta = beta(:);
% chi = chi(:);
beta = [0.1,1,10];
chi = [10,0.1,1];

for k = 1:numel(beta)
  str = ['Chi' num2str(chi(k),'%2.1e') '_ra065_beta' ...
      num2str(beta(k),'%2.1e')];

  str = strrep(str,'.','p');
  str = strrep(str,'-0','m');
  str = strrep(str,'+0','');
  shear1Ves(beta(k),chi(k),str);
end

