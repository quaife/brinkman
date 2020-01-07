%  subroutine to compute SPECTRAL derivatives in space for periodic data
function sx = fd1(x,n)

sx = fdiff(x,n);
sx(n+1) = sx(1);

end
