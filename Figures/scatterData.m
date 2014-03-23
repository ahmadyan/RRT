function Scatter(data,golden)


for i=3:6

figure(i);
hold on;
scatter(data(end,:),data(i,:));
plot(golden(end,:),golden(i,:),'g');
hold off;
end
