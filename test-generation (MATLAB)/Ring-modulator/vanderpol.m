function dxdt = josephson(t,x,p);
    M=2;
    
    dxdt=[
         M*x(1) - x(2) - x(1)^3 + p(1);
         x(1)+p(2)
        ]; 
end