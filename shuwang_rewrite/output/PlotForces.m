clc;clear;
addpath ..
file1 = 'Chi2p5_shax4p37_scL0p397_Conc0p8_Beta0_n1024_nbd1024_dt5en6_bmax1_bmin0p1_eps0p04_a100_contracting_left_contFromMid.bin';
[posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);
%Define a and epsilon from run
a = 100; % a is always 100
eps = 0.04;
%define bmin from file
 bmin = 0.1;
 bmax = 1;

irate = 1; 
istart = 1;
iend = numel(time1);

N = length(posx1(:,:,1));
oc = curve(N);

count = 1;
Eb = zeros(size(time1(istart:irate:iend)));
Et1 = zeros(size(time1(istart:irate:iend)));
Et2 = zeros(size(time1(istart:irate:iend)));
Etension = zeros(size(time1(istart:irate:iend)));
for k = istart:irate:iend    
    
%Compute b(u) in eq(10)
%bu = 1 * (ones(N,1) - conc1(:,:,k)) + bmin*conc1(:,:,k);
bu = (bmin - bmax)/2*tanh(3*(conc1(:,:,k) - 0.5)) + (bmax + bmin)/2;
%bu = bu*1e-19;
%Compute the curvature
X = [posx1(:,:,k);posy1(:,:,k)];
[~,~,cur] = oc.diffProp(X);
%Compute the length
[~, ~, Len] = oc.geomProp(X);
%Compute the double-well potential
fu = 0.25*conc1(:,:,k).^2.*(1-conc1(:,:,k)).^2;


%Compute the bending energy using the Trapezoid rule
Eb(count) = 0.5 * Len/N * sum(bu.*cur.^2);
EbStiff(count) = 0.5 * Len/N * sum(bu.*(bu>0.55).*cur.^2);
EbFloppy(count) = 0.5 * Len/N * sum(bu.*(bu<0.55).*cur.^2);
%Compute the line energy using the Trapezoid rule
du = oc.diffFT(conc1(:,:,k))/Len;
Et1(count) = a/eps * Len/N * sum(fu + 0*(eps^2 /2)*du.^2);
Et2(count) = a/eps * Len/N * sum(0*fu + (eps^2 /2)*du.^2);
ET = Et1 + Et2;
%Etension(count) = sum(-ten1(:,:,count) - a/eps*(fu-(eps^2/2)*du.^2))* Len/N;
Etension(count) = sum(-ten1(:,:,count))*Len/N;
count = count + 1;
end
% figure(7)
% plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Eb)
% hold on
% %semilogy(squeeze(mean(posx1(:,:,istart:irate:iend))),Et1)
% %semilogy(squeeze(mean(posx1(:,:,istart:irate:iend))),Et2)
% plot(squeeze(mean(posx1(:,:,istart:irate:iend))),ET)
% plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Etension)
% %xlim([0 ])
% legend('Eb', 'Ep', 'Et')
% %  title('Multi-component')
% title('Single-component')
% %title(titlestr)
% 
% posxRAp4NoConc = posx1;
% posyRAp4NoConc = posy1;
% EbRAp4NoConc = Eb;
% ETRAp4NoConc = ET;
% EtensionRAp4NoConc = Etension;
%save('EnergiesRAp4NoConc','posxRAp4NoConc','posyRAp4NoConc', 'EbRAp4NoConc','ETRAp4NoConc','EtensionRAp4NoConc')


       
