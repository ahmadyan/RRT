%http://www-users.math.umd.edu/~tvp/246/matlabode2.html
%f = inline('[ y(2); -sin(y(1))+sin(5*t) ]', 't', 'y');
%ode45(f,[0,20],[1;0]); %shows y(t) (blue) and y'(t) (green)


%[ts,ys] = ode45(f,[0,20],[1;0]); % find ts, ys, but don't show
%plot(ts,ys) % make plot of y1 and y2 vs. t
%[ts,ys] % show table with 3 columns for t, y1, y2

%[ts,ys] = ode45(f,[0,20],[1;0]); % find ts, ys, but don't show
%plot(ys(:,1),ys(:,2)) % make plot of y2 vs. y1


%f = inline('[y(2);-sin(y(1))]','t','y')
%vectfield(f,-2:.5:8,-2.5:.25:2.5)
%hold on
%for y20=0:0.3:2.7
%  [ts,ys] = ode45(f,[0,10],[0;y20]);
%  plot(ys(:,1),ys(:,2))
%end
%hold off


%f = inline('[(1-y(2))*y(1);(-1+y(1))*y(2)]','t','y')
%[ts,ys] = ode45(f,[0,20],[2;1]); 
%figure(1); plot(ts,ys)   
%figure(2); plot(ys(:,1),ys(:,2)); hold on
%vectfield(f,0:.2:2,0:.2:2); hold off

%competing species model
%x'=2x(1-x/2)-xy+u1
%y'=3y(1-y/3)-2xy+u2
%f = inline('[2*y(1)*(1-y(1)/2)-y(1)*y(2);3*y(2)*(1-y(2)/3)-2*y(1)*y(2)]','t','y')
%[ts,ys] = ode45(f,[0,20],[2;1]); 
%[ts,ys] = ode45(f,[0,20],[2;3]);
%figure(1); plot(ts,ys)   
%figure(2); plot(ys(:,1),ys(:,2)); hold on
%vectfield(f,0:.1:3,0:.1:3); hold off


%tunnel-diode
%x'=0.5*(y-(17.76*x-103.79*(x^2)+229.62*(x^3)-226.31*(x^4)+83.72*(x^5)))
%y'=0.2*(-x-1.5*y+1.2)
%-0.4:1.6;-0.4:1.6
f = inline('[0.5*(y(2)-(17.76*y(1)-103.79*(y(1)^2)+229.62*(y(1)^3)-226.31*(y(1)^4)+83.72*(y(1)^5)));0.2*(-y(1)-1.5*y(2)+1.2)]','t','y');
[ts,ys] = ode45(f,[0,20],[1;1]); 
figure(1); plot(ts,ys)   
figure(2); plot(ys(:,1),ys(:,2)); hold on
vectfield(f,-0.4:.05:1.6,-0.4:.05:1.6); hold off