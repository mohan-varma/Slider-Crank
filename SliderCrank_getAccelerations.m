function [Accl,alpha] = SliderCrank_getAccelerations(d2angA,Pos,omega)

alpha.A = [0;0;d2angA];
Accl.AO = cross(alpha.A,Pos.AO) + cross(omega.A,cross(omega.A,Pos.AO));
Accl.AB = cross(alpha.A,Pos.AB) + cross(omega.A,cross(omega.A,Pos.AB));


rvec = Pos.BC-Pos.AB;
Mat = [1 rvec(2,1);0 -rvec(1,1)];

vec = [Accl.AB(1,1);Accl.AB(2,1)];
Accl_w = cross(omega.B,cross(omega.B,rvec));
vec = vec - [Accl_w(1,1);Accl_w(2,1)];

sol = Mat\vec;


Abcx = sol(1,1);
d2angB = sol(2,1);
alpha.B = [0;0;d2angB];
Accl.BC = [Abcx;0;0];

Accl.CO = Accl.BC;
Accl.BO = 0.5*Accl.BC + 0.5*Accl.AB;

end
