%clc;clear;
addpath ..

file1 = 'Chi400_shax3p45_scL0p49_Conc0p3_Beta0_n1024_nbd1024_dt1en6_bmax1_bmin0p001_eps0p04_longchoke_newBending.bin';
titlestr = 'Chi = 0, RA_0 = 0.6, ubar = 0.3, beta = 0, bmin = 0.1';
[posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);

%Define a and epsilon from run
a = 100; % a is always 100
eps = 0.04;
%define bmin from file
bmin = 0.1;

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
% bu = 1 * (ones(N,1) - conc1(:,:,k)) + bmin*conc1(:,:,k);
bu = (bmin - 1)/2*tanh(3*(conc1(:,:,k) - 0.5)) + (1 + bmin)/2;
%Compute the curvature
X = [posx1(:,:,k);posy1(:,:,k)];
[~,~,cur] = oc.diffProp(X);
%Compute the length
[~, ~, Len] = oc.geomProp(X);
%Compute the double-well potential
fu = 0.25*conc1(:,:,k).^2.*(1-conc1(:,:,k)).^2;


%Compute the bending energy using the Trapezoid rule
Eb(count) = 0.5 * Len/N * sum(bu.*cur.^2);
%Compute the line energy using the Trapezoid rule
du = oc.diffFT(conc1(:,:,k))/Len;
Et1(count) = a/eps * Len/N * sum(fu + 0*(eps^2 /2)*du.^2);
Et2(count) = a/eps * Len/N * sum(0*fu + (eps^2 /2)*du.^2);
Etension(count) = sum(-ten1(:,:,count) - a/eps*(fu-(eps^2/2)*du.^2))*Len/N;
count = count + 1;
end
figure(3)
plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Eb)
hold on
plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Et1)
plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Et2)
plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Etension)
xlim([-24 28])
legend('Eb', 'Et1', 'Et2', 'Etension')
title(titlestr)




       