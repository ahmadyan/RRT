

mex discrepancy.cpp
d = 0.01 ;

discrepancy_result = [0,0];
discrepancy_index=1;
for x=0:d:1,
    for y=0:d:1,
        clear j;
        clear res ;
        clear up; clear lb;
        j=1;
        res = [0, 0];
        for i=1:3000,
            if(  (x <= rrt_sample(i, 1)) && (rrt_sample(i, 1)<= x+d) && (y <= rrt_sample(i, 2)) && (rrt_sample(i, 2)<= y+d)),
                res(j,:)= [(rrt_sample(i,1)-x)/d, (rrt_sample(i,2)-y)/d];
                j=j+1;
            end
        end
        if( size(res, 1)>20 ),
            x
            y
            [lb,  up] = discrepancy( res(:,1), res(:,2) );
            discrepancy_result(discrepancy_index, :) = [lb(1), up(1)];
            discrepancy_index = discrepancy_index+1;
        end
        
    end
end


