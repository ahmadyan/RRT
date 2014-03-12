function plotEnv(upper,lower,index)

hold on;
plot(upper(5,1:400),upper(index,1:400));
plot(lower(5,1:400),lower(index,1:400),'g');
hold off;


