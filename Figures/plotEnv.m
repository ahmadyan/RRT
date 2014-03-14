function plotEnv(upper,lower,index,t)

hold on;
plot(upper(t-2,:),upper(index,:));
plot(lower(t-2,:),lower(index,:),'g');
hold off;


