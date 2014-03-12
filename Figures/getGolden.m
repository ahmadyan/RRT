function Golden=getGolden(name)

r=rrt2mati(name);
Golden(:,1)=r(1:7,1)
for i=2:2037,
    Golden(:, i) = r(1:7, 2*(i-1));
end