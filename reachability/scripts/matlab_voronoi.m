xmin=-10;
xmax=10;
ymin=-10;
ymax=10;

samples = 100;
x=zeros(1,100);
y=zeros(1,100);
z=zeros(1,100);

for i=1:samples,
    x(i)= rand*(xmax-xmin);
    y(i)= rand*(ymax-ymin);
    z(i)=sys(x(i), y(i));
end

scatter(x,y)



[v,c]=voronoi(x, y); 
for i = 1:length(c) 
if all(c{i}~=1)   % If at least one of the indices is 1, 
                  % then it is an open region and we can't 
                  % patch that.
patch(v(c{i},1),v(c{i},2),i); % use color i.
end
end