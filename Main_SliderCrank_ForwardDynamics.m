% Forward Dynamics

clear all
clc
close all
global localvec Iner mass FG Torq lengthA lengthB


mass.A = 5;
mass.B = 7.5;
mass.C = 6;
lengthA = 1;
lengthB = 1.5;

Iner.A(3,3) = mass.A*(lengthA^2)/12;
Iner.B(3,3) = mass.B*(lengthB^2)/12;

localvec.AO_AN = [-0.5*lengthA;0;0];
localvec.AO_AB = [0.5*lengthA;0;0];
localvec.BO_AB = [-0.5*lengthB;0;0];
localvec.BO_BC = [0.5*lengthB;0;0];

FG.A = [0;-9.81;0]*mass.A;
FG.B = [0;-9.81;0]*mass.B;
FG.C = [0;-9.81;0]*mass.C;

Torq.A = [0;0;20];

% Initial Values
angA0 = deg2rad(45); %rad
dangA0 = 12*pi; % rad/sec

initvalues = [angA0;dangA0];
[T,Y] = ode45(@SliderCrank_NE,[0 10],initvalues);
[lengthtime,~] = size(T);

for ii = 1:lengthtime

    angA = Y(ii,1);
    dangA = Y(ii,2);

    dY = SliderCrank_NE(T(ii),Y(ii,1:2)');
    d2angA = dY(2,1);
    Pos = SliderCrank_getPositions(angA);
    [Vel,omega] = SliderCrank_getVelocities(dangA,Pos);
    [Accl,alpha] = SliderCrank_getAccelerations(d2angA,Pos,omega);

r1A = -Pos.AO;
r2A = Pos.AB-Pos.AO;
r1B = Pos.AB-Pos.BO;
r2B = Pos.BC-Pos.BO;

r1Atilde = make_tilde(r1A);
r2Atilde = make_tilde(r2A);
r1Btilde = make_tilde(r1B);
r2Btilde = make_tilde(r2B);

    MAT = zeros(8,7);
    % FNA + FAB - massA*A_AO(unknown part) = -FGA + massA*A_AO(known part)
    MAT(1:2,1:2) = eye(2,2);
    MAT(1:2,3:4) = eye(2,2);
    VEC(1:2,1) = mass.A*Accl.AO(1:2,1) -FG.A(1:2,1);

    % -FAB + FBC - massB*A_BO(unknown part) = -FGB + massB*A_BO(known part)
    MAT(3:4,3:4) = -eye(2,2);
    MAT(3:4,5:6) = eye(2,2);
    VEC(3:4,1) = mass.B*Accl.BO(1:2,1) -FG.B(1:2,1);

    % -FBC + FNC - massC*A_CO(unknown part) = -FGC + massC*A_CO(known part)
    MAT(5:6,5:6) = -eye(2,2);
    MAT(6,7) = 1;
    VEC(5:6,1) = mass.C*Accl.CO(1:2,1) -FG.C(1:2,1);

    MAT(7,:) = [r1Atilde(3,1:2) r2Atilde(3,1:2) 0 0 0];
    VEC(7,1) = -Torq.A(3,1) + Iner.A(3,3)*alpha.A(3,1);

    MAT(8,:) = [0 0 -r1Btilde(3,1:2) r2Btilde(3,1:2) 0];
    VEC(8,1) =  Iner.B(3,3)*alpha.B(3,1);


    VE = MAT'*VEC;
    MA = (MAT'*MAT);

    sol = MA\VE;


    FNA(:,ii) = sol(1:2,1);
    FAB(:,ii) = sol(3:4,1);
    FBC(:,ii) = sol(5:6,1);
    FNC(:,ii) = [0;sol(7)];

    linexy = [Pos.AN Pos.AB Pos.BC];

    figure(1)
    clf
    plot(linexy(1,:),linexy(2,:))
    axis([-1 3 -2 2])

end
%

figure(2)
subplot(221)
plot(T,0.001*FNA)
title('time vs FNAx and FNAy (kN)')
subplot(222)
plot(T,0.001*FAB)
title('time vs FABx and FABy(kN)')
subplot(223)
plot(T,0.001*FBC)
title('time vs FBCx and FBCy(kN)')
subplot(224)
plot(T,0.001*FNC)
title('time vs FNCx and FNCy(kN)')
 
