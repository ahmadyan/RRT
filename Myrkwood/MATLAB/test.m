for i=1:num_nodes,
    if tree(i).x>1.3, 
        tree(i).x=1.3 ;
    end
    if tree(i).x<-0.3, 
        tree(i).x=-0.3;
    end
    if tree(i).y>1.3, tree(i).y=1.3 ;end
    if tree(i).y<-0.3, tree(i).y=-0.3; end
end

for i=1:num_nodes,
    vx(i) = tree(i).x ;
    vy(i) = tree(i).y ;
end
figure(2)

[vxp, vyp] = voronoi(vx,vy);
plot(vx, vy, 'r+', vxp, vyp, 'b-'); axis equal

figure(3)
vv = [vx;vy]' ;
[v,c] = voronoin(vv);
for i=1:length(c)
        if tree(i).direction==+1 %all(c{i}~= 1)
                patch(v(c{i},1), v(c{i},2),'y');
        else
                patch(v(c{i},1), v(c{i},2),'w');
        end
end