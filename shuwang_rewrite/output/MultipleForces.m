addpath ..

ax = [-50 50 -4 4];

% file1 = 'Chi2p5_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_rand.bin';
% [posx1,posy1,conc1,ea1,el1,time1,xvel1,yvel1,ten1] = loadFile(file1);
% file2 = 'Chi2p5_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_right.bin';
% [posx2,posy2,conc2,ea2,el2,time2,xvel2,yvel2,ten2] = loadFile(file2);
% file3 = 'Chi2p5_shax5p7_scL0p311_Conc0p5_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_left.bin';
% [posx3,posy3,conc3,ea3,el3,time3,xvel3,yvel3,ten3] = loadFile(file3);
% file4 = 'Chi2p5_shax5p7_scL0p311_Conc0_Beta0_n1024_nbd1024_dt1en4_bmax1_bmin0p1_eps0p04_a100_contracting_rand.bin';
% [posx4,posy4,conc4,ea4,el4,time4,xvel4,yvel4,ten4] = loadFile(file4);

%define a and epsilon from run
a = 100; % a is always 100
eps = 0.04;
%define bmin from file
bmin = 0.1;

irate = 1; 
istart = 1;
iend = min([numel(time1),numel(time2),numel(time3),numel(time4)]);
ntime = iend;
N = length(posx1(:,:,1));
oc = curve(N);
 
count = 1;
for k = istart:irate:iend
    %Compute b(u) in eq(10)
    bu1 = (bmin - 1)/2*tanh(3*(conc1(:,:,k) - 0.5)) + (1 + bmin)/2;
    bu2 = (bmin - 1)/2*tanh(3*(conc2(:,:,k) - 0.5)) + (1 + bmin)/2;
    bu3 = (bmin - 1)/2*tanh(3*(conc3(:,:,k) - 0.5)) + (1 + bmin)/2;
    bu4 = (bmin - 1)/2*tanh(3*(conc4(:,:,k) - 0.5)) + (1 + bmin)/2;

    %Compute the curvature
    X1 = [posx1(:,:,k);posy1(:,:,k)];
    X2 = [posx2(:,:,k);posy2(:,:,k)];
    X3 = [posx3(:,:,k);posy3(:,:,k)];
    X4 = [posx4(:,:,k);posy4(:,:,k)];
    [~,~,cur1] = oc.diffProp(X1);
    [~,~,cur2] = oc.diffProp(X2);
    [~,~,cur3] = oc.diffProp(X3);
    [~,~,cur4] = oc.diffProp(X4);
    %Compute the length
    [~, ~, Len1] = oc.geomProp(X1);
    [~, ~, Len2] = oc.geomProp(X2);
    [~, ~, Len3] = oc.geomProp(X3);
    [~, ~, Len4] = oc.geomProp(X4);
    %Compute the double-well potential
    fu1 = 0.25*conc1(:,:,k).^2.*(1-conc1(:,:,k)).^2;
    fu2 = 0.25*conc2(:,:,k).^2.*(1-conc2(:,:,k)).^2;
    fu3 = 0.25*conc3(:,:,k).^2.*(1-conc3(:,:,k)).^2;
    fu4 = 0.25*conc4(:,:,k).^2.*(1-conc4(:,:,k)).^2;

    %Compute the bending energy using the Trapezoid rule
    Eb1(count) = 0.5 * Len1/N * sum(bu1.*cur1.^2);
    Eb2(count) = 0.5 * Len2/N * sum(bu2.*cur2.^2);
    Eb3(count) = 0.5 * Len3/N * sum(bu3.*cur3.^2);
    Eb4(count) = 0.5 * Len4/N * sum(bu4.*cur4.^2);
    
    %Compute the line energy using the Trapezoid rule
    du1 = oc.diffFT(conc1(:,:,k))/Len1;
    du2 = oc.diffFT(conc2(:,:,k))/Len2;
    du3 = oc.diffFT(conc3(:,:,k))/Len3;
    du4 = oc.diffFT(conc4(:,:,k))/Len4;
    Et1(count) = a/eps * Len1/N * sum(fu1 + 0*(eps^2 /2)*du1.^2) + ...
                 a/eps * Len1/N * sum(0*fu1 + (eps^2 /2)*du1.^2);
    Et2(count) = a/eps * Len2/N * sum(fu2 + 0*(eps^2 /2)*du2.^2) + ...
                 a/eps * Len2/N * sum(0*fu2 + (eps^2 /2)*du2.^2);
    Et3(count) = a/eps * Len3/N * sum(fu3 + 0*(eps^2 /2)*du3.^2) + ...
                 a/eps * Len3/N * sum(0*fu3 + (eps^2 /2)*du3.^2);
    Et4(count) = a/eps * Len4/N * sum(fu4 + 0*(eps^2 /2)*du4.^2) + ...
                 a/eps * Len4/N * sum(0*fu4 + (eps^2 /2)*du4.^2);  
   
    %compute the tension
    Etens1(count) = sum(-ten1(:,:,count) - a/eps*(fu1-(eps^2/2)*du1.^2));
    Etens2(count) = sum(-ten2(:,:,count) - a/eps*(fu2-(eps^2/2)*du2.^2));
    Etens3(count) = sum(-ten3(:,:,count) - a/eps*(fu3-(eps^2/2)*du3.^2));
    Etens4(count) = sum(-ten4(:,:,count) - a/eps*(fu4-(eps^2/2)*du4.^2));
    
    %Update counter
    count = count + 1;
end

subplot(2,2,1)
semilogy(mean(squeeze(posx1(:,:,istart:irate:iend))),Eb1)
hold on
semilogy(mean(squeeze(posx2(:,:,istart:irate:iend))),Eb2)
semilogy(mean(squeeze(posx3(:,:,istart:irate:iend))),Eb3)
semilogy(mean(squeeze(posx4(:,:,istart:irate:iend))),Eb4)
xlim([0 19])
title("Bending energy")
legend('random', 'right', 'left', 'no conc')

subplot(2,2,2)
semilogy(squeeze(mean(posx1(:,:,istart:irate:iend))),Et1)
hold on
semilogy(squeeze(mean(posx2(:,:,istart:irate:iend))),Et2)
semilogy(squeeze(mean(posx3(:,:,istart:irate:iend))),Et3)
semilogy(squeeze(mean(posx4(:,:,istart:irate:iend))),Et4)
xlim([0 19])
title("Phase energy")
legend('random', 'right', 'left', 'no conc')

subplot(2,2,3)
semilogy(squeeze(mean(posx1(:,:,istart:irate:iend))),abs(Etens1))
hold on
semilogy(squeeze(mean(posx2(:,:,istart:irate:iend))),abs(Etens2))
semilogy(squeeze(mean(posx3(:,:,istart:irate:iend))),abs(Etens3))
semilogy(squeeze(mean(posx4(:,:,istart:irate:iend))),abs(Etens4))
xlim([0 19])
title("Tension")
legend('random', 'right', 'left', 'no conc')

% subplot(2,2,1)
% plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Eb1)
% hold on
% plot(squeeze(mean(posx2(:,:,istart:irate:iend))),Eb2)
% plot(squeeze(mean(posx3(:,:,istart:irate:iend))),Eb3)
% plot(squeeze(mean(posx4(:,:,istart:irate:iend))),Eb4)
% title("Bending energy")
% legend('random', 'right', 'left', 'no conc')
% 
% subplot(2,2,2)
% plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Et1)
% hold on
% plot(squeeze(mean(posx2(:,:,istart:irate:iend))),Et2)
% plot(squeeze(mean(posx3(:,:,istart:irate:iend))),Et3)
% plot(squeeze(mean(posx4(:,:,istart:irate:iend))),Et4)
% title("Phase energy")
% legend('random', 'right', 'left', 'no conc')
% 
% subplot(2,2,3)
% plot(squeeze(mean(posx1(:,:,istart:irate:iend))),Etens1)
% hold on
% plot(squeeze(mean(posx2(:,:,istart:irate:iend))),Etens2)
% plot(squeeze(mean(posx3(:,:,istart:irate:iend))),Etens3)
% plot(squeeze(mean(posx4(:,:,istart:irate:iend))),Etens4)
% title("Tension")
% legend('random', 'right', 'left', 'no conc')
