%solution of BVP problem y''=12x^2 with boundary conditions y(0)=0, y(1)=0
%exact solution is y(x)=x^4-x
%I approximate y''=12x^2 by taylor:
%Y_{j-1}-2Y_j+Y_{j+1}/\delta x^2 = 12x^2_j

J=10;
dx=1/J;
x=(0:dx:1);
b=zeros(J+1,1);
b(2:J)=12*dx^2*x(2:J).^2;
A   =  sparse(J+1,J+1);                          
A(1,1)       =  1.0;       
A(J+1,J+1)         =   1.0;
for j=2:J,
    A(j, [j-1, j, j+1])=[1,-2,1];
end


%solving the matrix problem looks like: 
Y  =   A  \   b;       %  solve      A  Y   =  b 
%  also     get    exact      soln     on   fine     grid: 
xf   =  0:1/1000:1;
yexact       =  xf.^4 -  xf; 
plot(x,Y,'o',xf,yexact) 
grid    on,    xlabel       x,   legend('finite             diff','exact') 