function [Vel,omega] = SliderCrank_getVelocities(dangA,Pos)

omega.A = [0;0;dangA];
Vel.AO = cross(omega.A,Pos.AO);
Vel.AB = cross(omega.A,Pos.AB);

rvec = Pos.BC-Pos.AB;
Mat = [1 rvec(2,1);0 -rvec(1,1)];
vec = [Vel.AB(1,1);Vel.AB(2,1)];

sol = Mat\vec;

Vbcx = sol(1,1);
dangB = sol(2,1);
omega.B = [0;0;dangB];
Vel.BC = [Vbcx;0;0];

Vel.CO = Vel.BC;
Vel.BO = 0.5*Vel.BC + 0.5*Vel.AB;
end

