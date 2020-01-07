function kmatrix=formkmatrix(m)

kmatrix=ones(2*m,2*m);
for i=1:2*m
  for j=1:2*m
  if mod(abs(i-j),2)<0.5
    kmatrix(i,j)=0;
  end
end

end

end
