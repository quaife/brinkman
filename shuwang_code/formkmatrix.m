function kmatrix=formkmatrix(m)

% speedup
kmatrix = ones(2,2*m);
kmatrix(1,1:2:end) = 0;
kmatrix(2,2:2:end) = 0;
kmatrix = repmat(kmatrix,m,1);

%tic
%kmatrix=ones(2*m,2*m);
%for i=1:2*m
%  for j=1:2*m
%    if mod(abs(i-j),2)<0.5
%      kmatrix(i,j)=0;
%    end
%  end
%end
%toc


end
