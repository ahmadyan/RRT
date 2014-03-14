function numberOfsamplesAgainstTime=count_sample(dev)

NOSAT=[1;0]
c=(sortrows(dev',7))'; %%sort data w.r.t. time
k=2;
[x,y]=size(c);
dt=2e-7;
while (k<y+1)
    t=k;%%pointing at the start of a time 
    if (k<20000)
            while (c(end,k+1)==c(end,k))
            k=k+1
            end;
    end;
    
    NOSAT=[NOSAT,[k-t+1;(k-t+1)*dt]];
    k=k+1;
    
end
numberOfsamplesAgainstTime=NOSAT;
end
    
    