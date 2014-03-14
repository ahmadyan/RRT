function [upper,lower,sorted]=getEnv(Dev_t)

c=(sortrows(Dev_t,16))'; %%sort data w.r.t. time
k=1;
[x,y]=size(c);
upper=c(3:16,k);
lower=c(3:16,k);
sorted=c;
while (k<3001)
    temp_start=k%%pointing at the start of a time 
    
    if k<y
   while (c(end,k+1)==c(end,k)) && (k<y-1)
            k=k+1;
   end;    
    end;
    if (c(end,k)==c(end,y))
            k=k+1;
    end;
        
    new_upper_column= [max(c(3,temp_start:k));
                       max(c(4,temp_start:k));
                       max(c(5,temp_start:k));
                       max(c(6,temp_start:k));
                       max(c(7,temp_start:k));
                       max(c(8,temp_start:k));
                       max(c(9,temp_start:k));
                       max(c(10,temp_start:k));
                       max(c(11,temp_start:k));
                       max(c(12,temp_start:k));
                       max(c(13,temp_start:k));
                       max(c(14,temp_start:k));
                       max(c(15,temp_start:k));
                       max(c(16,temp_start:k))];
                   
    if(upper(end,end)~=new_upper_column(end,1))               
        upper=[upper,new_upper_column];
    end;
    
    new_lower_column= [min(c(3,temp_start:k));
                       min(c(4,temp_start:k));
                       min(c(5,temp_start:k));
                       min(c(6,temp_start:k));
                       min(c(7,temp_start:k));
                       min(c(8,temp_start:k));
                       min(c(9,temp_start:k));
                       min(c(10,temp_start:k));
                       min(c(11,temp_start:k));
                       min(c(12,temp_start:k));
                       min(c(13,temp_start:k));
                       min(c(14,temp_start:k));
                       min(c(15,temp_start:k));
                       min(c(16,temp_start:k))];
    if(lower(end,end)~=new_lower_column(end,1))               
        lower=[lower,new_lower_column];
    
    end
   t=new_upper_column
   g=new_lower_column
    k
    k=k+1;
    
end
end
    
    