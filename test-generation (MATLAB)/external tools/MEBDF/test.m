global variation0
global variation1

t_end=1e-3;
t_start=0;
h=-1e-3/1000;
x=zeros(1,15);
y=zeros(1,1000);
t=0;
for i=1:1000,
variation0 = 2*(rand-0.5)/10;
variation1 = 2*(rand-0.5)/20;
	i
    h
	x = ring(x, t, h);
	t=t+h;
	y(i)=x(8);
end