sample = zeros(200,2);
count=1;

%tunnel-diode
%x'=0.5*(y-(17.76*x-103.79*(x^2)+229.62*(x^3)-226.31*(x^4)+83.72*(x^5)))
%y'=0.2*(-x-1.5*y+1.2)
%-0.4:1.6;-0.4:1.6
xmin=-0.4; 
xmax=1.6 ;
ymin=-0.4;
ymax=1.6;
resolution=10;
sim_time=0.7 ;
f = inline('[0.5*(y(2)-(17.76*y(1)-103.79*(y(1)^2)+229.62*(y(1)^3)-226.31*(y(1)^4)+83.72*(y(1)^5)));0.2*(-y(1)-1.5*y(2)+1.2)]','t','y');

hold on
for i=resolution*xmin:resolution*xmax,
    for j=resolution*ymin:resolution*ymax,
        [ts,ys] = ode45(f,[0,sim_time+3*rand],[i*(1/resolution);j*(1/resolution)]); 
        plot(ys(:,1),ys(:,2)); 
        sample(count, 1)=ys(size(ys,1),1);
        sample(count, 2)=ys(size(ys,1),2);
        count=count+1;
    end
end
vectfield(f,xmin:.05:xmax,ymin:.05:ymax); 
figure(2)
scatter(sample(:,1), sample(:,2))
hold off
