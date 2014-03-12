function [upper,lower,sorted]=getEnv(data)

c=(sortrows(data',7))'; %%sort data w.r.t. time
k=2;
[x,y]=size(data)
upper=c(3:7,1);
lower=c(3:7,1);

sorted=c;
while (k<20001)
    newarray=c(3:7,k);
    
    if(upper(end,end)~=newarray(end,1))               
        upper=[upper,newarray];
    else     
        upper(1,end)=max(upper(1,end),newarray(1,1));
        upper(2,end)=max(upper(2,end),newarray(2,1));
        upper(3,end)=max(upper(3,end),newarray(3,1));
        upper(4,end)=max(upper(4,end),newarray(4,1));
    end;
        
     if(lower(end,end)~=newarray(end,1))               
        lower=[lower,newarray];
    else
        lower(1,end)=min(upper(1,end),newarray(1,1));
        lower(2,end)=min(upper(2,end),newarray(2,1));
        lower(3,end)=min(upper(3,end),newarray(3,1));
        lower(4,end)=min(upper(4,end),newarray(4,1));
     end;
    k=k+1
    
end
end
    
    