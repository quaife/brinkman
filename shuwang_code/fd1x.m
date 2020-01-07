% remove linear part
function sx = fd1x(x,n)

N = (0:1:n-1); 
x1(1,1:n) = x(1,1:n)-2*pi*N/n;
x1(1,n+1) = x(1,1);
sx = fdiff(x1,n);
sx = sx + pi*2;
sx(1,n+1) = sx(1,1);
       
end
