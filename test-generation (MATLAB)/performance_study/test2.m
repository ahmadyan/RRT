[X,Y] = meshgrid(0:.1:2.5, 0:.1:2.5);
U1=X.*(1.5-.5*X-Y);
V1=Y.*(2-Y-1.125*X);
U2=X.*(1-.5*Y);
V2=Y.*(-.25-.5*X);
warning off all
%Begin simulation
for alp=0:.01:0.01
    figure;hold on%Create new figure
    U=(1-alp)*U1+alp*U2;
    V=(1-alp)*V1+alp*V2;
    L=sqrt(U.^2+V.^2);
    quiver(X,Y, U./L,V./L,.5)%plot vector field
    axis equal 
    %Phase plot equation
    k=@(t,x)[(1-alp)*x(1)*(1.5-.5*x(1)-x(2))+alp*x(1)*(1-.5*x(2));(1-alp)*x(2)*(2-x(2)-1.125*x(1))+alp*x(2)*(-.25-.5*x(1))];
    for a=-2:.25:2%Phase plot loop
        for b=-2:.25:2
            [t,xa]=ode45(k,[0 10],[a b]);
            plot(xa(:,1),xa(:,2),'r');
            [t,xa]=ode45(k,[0 -5],[a b]);
            plot(xa(:,1),xa(:,2),'r');
         end
    end
    axis([0 2.5 0 2.5])
    title(alp)
end