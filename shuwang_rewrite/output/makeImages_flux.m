%clc;clear;close all
addpath ..

 file = 'Chi200_RA0p75_Conc0p3_Beta0_y0p1_eps0p04_n20.bin';
%file = 'Test.bin';
farFieldSpeed = 200;
farFieldFlow = 'parabolic';

[posx,posy,conc,ea,el,time,xvel1,yvel1,lambTil] = loadFile(file);
N = numel(posx(:,:,end));
oc = curve(N);
op = poten(N);

istart = 2;
irate = 100; 
iend = numel(time);

for k = istart:irate:iend
    clf; 
    [~,A,~] = oc.geomProp([posx(:,:,k);posy(:,:,k)]);
    [len,theta,cur] = oc.computeOpeningAngle(N,[posx(:,:,k);posy(:,:,k)]);
    %Compute the curvature
    ves = capsules([posx(:,:,k);posy(:,:,k)],conc(:,:,k));
    ves.bendsti = 1;
    ves.bendratio = 0.1;
    % Compute Eu and Esigma, equations (13) and (14) 
    [Eu,Esigma] = ves.variationsNonStiff;
    
    Eu = -Eu;
    % Compute the force in eq (33)
    %   NOTE: this does not include u_s term.
    tau = [[-Esigma.*sin(theta) - Eu.*cos(theta)]; ...
           [Esigma.*cos(theta) - Eu.*sin(theta)]];
    % Calculate the Fourier derivative of lambdaTilde
    dlamTilds = oc.diffFT(lambTil(:,:,k))/len;
    % We can now compute the traction jump in first part of equation (39).
    % This comes from applying the product rule and using Frenet-Seret.
    tracJump = [(+lambTil(:,:,k).*cur.*sin(theta) - dlamTilds.*cos(theta));...
                (-lambTil(:,:,k).*cur.*cos(theta) - dlamTilds.*sin(theta))];
    % Now we add the jump conditions in eq (39) to (33) which is in the 
    % old variable tau
    tau = tau + tracJump;

    x = -0.8:0.05:1.2;
    y = -1.2:0.05:1.2;
    [xtar,ytar] = meshgrid(x,y);
    xtar = xtar(:); ytar = ytar(:);
    [vel] = op.StokesSLPtar(ves,tau,[xtar;ytar]);
    [velx,vely] = oc.getXY(vel);

    if strcmp(farFieldFlow,'shear')
        velx = velx + farFieldSpeed*ytar;
    elseif strcmp(farFieldFlow, 'parabolic')
        R0 = sqrt(A/pi);
        W = 10*R0; 
        velx = velx + farFieldSpeed*(1-(ytar/W).^2);
%         %clf; plot(velx);hold on;
%         dt = 1e-5;
%         dxdt = (posx(:,:,k+1)-posx(:,:,k))/dt + farFieldSpeed*(1-(posy(:,:,k)/W).^2);
%         velx = velx - (1/N)*sum(dxdt);
%        vesvel = xvel1(:,:,k) + farFieldSpeed*(1-(posy(:,:,k)/W).^2);
        velx = velx - mean(xvel1(:,:,k));
%         plot(velx,'r--')
%         pause
    elseif strcmp(farFieldFlow,'extensional')
        velx = velx + farFieldSpeed * (-xtar);
        vely = vely + farFieldSpeed * ytar;
    end 

    quiver(xtar,ytar,velx,vely)
    hold on
    px = posx(:,:,k);
    py = posy(:,:,k);
    plot(px,py)
    plot(px(1),py(1),'ko')
    pause(0.1)
end

