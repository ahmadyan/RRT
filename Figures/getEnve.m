function [upper,lower,sorted]=findEnv(data)

a=(sortrows(data',[7,3]))'; %%sort data w.r.t. time
b=(sortrows(data',[7,4]))';
c=(sortrows(data',[7,5]))';
d=(sortrows(data',[7,6]))';


k=2;
[x,y]=size(data)
upper=c(3:7,1);
lower=c(3:7,1);
sorted=c;
while (k<y+1)
    t=k;%%pointing at the start of a time 
    if (k<y)
            while (c(end,k+1)==c(end,k))
            k=k+1
            end;
    end;
   maxarray=[a(3,k);b(4,k);c(5,k);d(6:7,k)];
   minarray=[a(3,t);b(4,t);c(5,t);d(6:7,t)];
    upper=[upper,maxarray];
    lower=[lower,minarray];
    
   
    k=k+1
    
end
end
    
    