function XI = synth1(E,N,Z)

N2 = N/2;
N22 = N2-2;

CXP1 = exp(1i*E);
CXP2 = conj(CXP1);
U1 = Z(N2);
U2 = Z(N2+2);

for j = 1:N22
  ZJ1 = Z(N2-j);
  % ZJ1
  ZJ2 = Z(N2+j+2);
  % ZJ2
  U1 = ZJ1 + CXP1*U1;
  U2 = ZJ2 + CXP2*U2;
end
XI = real(Z(1))+real(Z(N2+1))*cos(N2*E) + ...
             real(CXP1*U1+CXP2*U2);  

end
