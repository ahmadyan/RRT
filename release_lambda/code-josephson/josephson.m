function dxdt = josephson(t,x,p);
    C=1;    %pF
    G=0.25;    %uH
    I0=1;  %kOhm
    K=1;  %V
    IS=0;
    
    %X(1) is V_C
    %X(2) is Phi_L
    VC=x(1);
    PL=x(2);
    
    Is=p(1);
    PV=p(2);
    
    dxdt=[
        (1/C)*(-G*VC-I0*sin(K*(PL+PV))+Is);
        VC
        ]; 
end