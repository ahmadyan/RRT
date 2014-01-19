m

mex discrepancy.cpp
d = 0.2 ;

discrepancy_result = [0,0];
discrepancy_index=1;
for x=-2:d:0,
    for y=1.5:d:3,
        mex discrepancy.cpp
        clear j;
        clear res ;
        clear up; clear lb;
        j=1;
        res = [0, 0];
        for i=1:3000,
           % x
           % y
           % rrt_sample(i, 1)
           % rrt_sample(i, 2)            
            
            if(  (x <= rrt_sample(i, 1)) && (rrt_sample(i, 1)<= x+d) && (y <= rrt_sample(i, 2)) && (rrt_sample(i, 2)<= y+d)),
                res(j,:)= [(rrt_sample(i,1)-x)/d, (rrt_sample(i,2)-y)/d];
                j=j+1;
            end
        end
        if( size(res, 1)>20 ),
            [lb,  up] = discrepancy( res(:,1), res(:,2) );
            discrepancy_result(discrepancy_index, :) = [lb(1), up(1)];
            discrepancy_index = discrepancy_index+1;
        end
        
    end
end


x=-0.1;
y=-0.1;
d = 0.2 ;
        mex discrepancy.cpp
        clear j;
        clear res ;
        clear up; clear lb;
        j=1;
        res = [0, 0];
        for i=1:3000,
           % x
           % y
           % rrt_sample(i, 1)
           % rrt_sample(i, 2)            
            
            if(  (x <= rrt_sample(i, 1)) && (rrt_sample(i, 1)<= x+d) && (y <= rrt_sample(i, 2)) && (rrt_sample(i, 2)<= y+d)),
                res(j,:)= [(rrt_sample(i,1)-x)/d, (rrt_sample(i,2)-y)/d];
                j=j+1;
            end
        end
        if( size(res, 1)>0 ),
            [lb,  up] = discrepancy( res(:,1), res(:,2) );
            discrepancy_result(discrepancy_index, :) = [lb(1), up(1)];
            discrepancy_index = discrepancy_index+1;
        end
        size(res)
        
        c=0;
        for i=1:3000,
            xb = rrt_sample(i, 1);
            yb = rrt_sample(i, 2);
            if( sqrt( xb^2 + yb^2 ) < 0.1 ),
                c = c+1; 
            end
        end
        c