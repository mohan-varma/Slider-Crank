function Pos = SliderCrank_getPositions(angA)
global localvec lengthB

    N_A = [cos(angA) -sin(angA) 0;
           sin(angA)  cos(angA) 0;
           0 0 1];


    Pos.AN = [0;0;0];
    Pos.AO = -N_A*localvec.AO_AN;
    Pos.AB = Pos.AO + N_A*localvec.AO_AB;

    [xc,yc] = linecirc(0,0,Pos.AB(1,1),Pos.AB(2,1),lengthB);
    pnt1 = [xc(1);yc(1);0];
    pnt2 = [xc(2);yc(2);0];


        if(xc(1)>0)
            Pos.BC = pnt1;
        elseif(xc(2)>0)
            Pos.BC = pnt2;
        else
            Pos.BC = [];
        end


        Pos.BO = 0.5*Pos.AB + 0.5*Pos.BC;
        Pos.CO = Pos.BC;

end
