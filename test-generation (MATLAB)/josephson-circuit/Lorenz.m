function dxdt=Lorenz(t,x)
% The RHS of the Lorenz attractor
% Save this function in a separate file 'myLorenz.m'
sigma = 10;
r = 28;
b = 8/3;
dxdt=[ sigma*(x(2)-x(1));
(1+r)*x(1)-x(2)-x(1)*x(3);
x(1)*x(2)-b*x(3)];
end
