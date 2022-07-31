function dyvar = SliderCrank_NE(~,yvar)

global Iner mass FG Torq



angA = yvar(1,1);
dangA = yvar(2,1);
dyvar(1,1) = dangA;

Pos = SliderCrank_getPositions(angA);
[Vel,omega] = SliderCrank_getVelocities(dangA,Pos);
dangA_temp = 1;
[CVel,Comg] = SliderCrank_getVelocities(dangA_temp,Pos);

d2angA_temp = 0;
[Accl_k,alpha_k] = SliderCrank_getAccelerations(d2angA_temp,Pos,omega);

r1A = -Pos.AO;
r2A = Pos.AB-Pos.AO;
r1B = Pos.AB-Pos.BO;
r2B = Pos.BC-Pos.BO;

r1Atilde = make_tilde(r1A);
r2Atilde = make_tilde(r2A);
r1Btilde = make_tilde(r1B);
r2Btilde = make_tilde(r2B);

    MAT = zeros(7,7);
    % FNA + FAB - massA*A_AO(unknown part) = -FGA + massA*A_AO(known part)
    MAT(1:2,1:2) = eye(2,2);
    MAT(1:2,3:4) = eye(2,2);
    MAT(1:2,8) = -mass.A*CVel.AO(1:2,1);
    VEC(1:2,1) = -FG.A(1:2,1) + mass.A*Accl_k.AO(1:2,1);

    % -FAB + FBC - massB*A_BO(unknown part) = -FGB + massB*A_BO(known part)
    MAT(3:4,3:4) = -eye(2,2);
    MAT(3:4,5:6) = eye(2,2);
    MAT(3:4,8) = -mass.B*CVel.BO(1:2,1);
    VEC(3:4,1) = -FG.B(1:2,1) + mass.B*Accl_k.BO(1:2,1);

    % -FBC + FNC - massC*A_CO(unknown part) = -FGC + massC*A_CO(known part)
    MAT(5:6,5:6) = -eye(2,2);
    MAT(6,7) = 1;
    MAT(5:6,8) = -mass.C*CVel.CO(1:2,1);
    VEC(5:6,1) = -FG.C(1:2,1) + mass.C*Accl_k.CO(1:2,1);

    MAT(7,:) = [r1Atilde(3,1:2) r2Atilde(3,1:2) 0 0 0 -Iner.A(3,3)*alpha_k.A(3,1)];
    VEC(7,1) = -Torq.A(3,1);

    MAT(8,:) = [0 0 -r1Btilde(3,1:2) r2Btilde(3,1:2) 0 -Iner.B(3,3)*alpha_k.B(3,1)];
    VEC(8,1) = 0;

    sol = MAT\VEC;
    d2angA = sol(8,1);

    dyvar(2,1) = d2angA;
end
