function [upper,lower,sorted]=findEnv(data)

c=(sortrows(data',7))'; %%sort data w.r.t. time
k=1;
[x,y]=size(data)
upper=c(3:7,k);
lower=c(3:7,k);
sorted=c;
while (k<20001)
    temp_start=k;%%pointing at the start of a time 
    if (k<20000)
            while (c(end,k+1)==c(end,k))
            k=k+1
            end;
    end;
    array=c(3,temp_start:k)
    max_value=max(c(3,temp_start:k))
    min_value=min(c(3,temp_start:k))
    new_upper_column= [max(c(3,temp_start:k));
                       max(c(4,temp_start:k));
                       max(c(5,temp_start:k));
                       max(c(6,temp_start:k));
                       max(c(7,temp_start:k))];
    if(upper(end,end)~=new_upper_column(end,1))               
        upper=[upper,new_upper_column];
    end;
    new_lower_column= [min(c(3,temp_start:k));
                       min(c(4,temp_start:k));
                       min(c(5,temp_start:k));
                       min(c(6,temp_start:k));
                       min(c(7,temp_start:k))];
    if(lower(end,end)~=new_lower_column(end,1))               
        lower=[lower,new_lower_column];
    
    end
   t=new_upper_column;
   g=new_lower_column;
    
    k=k+1
    
end
end
    
    