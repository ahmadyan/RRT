function plotEnv(upper,lower)


for i=1:4

figure(i);
hold on;
plot(upper(end,:),upper(i,:));
plot(lower(end,:),lower(i,:),'g');
hold off;
end



