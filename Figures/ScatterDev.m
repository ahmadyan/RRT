function Scatter(Dev,golden)


for i=4:6

figure(i);
hold on;
scatter(Dev(end,:),Dev(i,:));
plot(golden(end,:),golden(i,:),'g');
hold off;
end
