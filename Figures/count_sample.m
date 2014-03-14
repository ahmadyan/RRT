function numberOfsamplesAgainstTime=count_sample(dev)

NOSAT=[1;0];
c=(sortrows(dev,16))'; %%sort data w.r.t. time
k=2;
[x,y]=size(c);
dt=1e-5;
counttime=0;
while (k<y+1)
    t=k;%%pointing at the start of a time 

    if k<y
   while (c(end,k+1)==c(end,k)) && (k<y-1)
            k=k+1;
   end;    
    end;
    if (c(end,k)==c(end,y))
            k=k+1;
    end;
    NOSAT=[NOSAT,[k-t+1;counttime*dt]];
    k=k+1
    counttime=counttime+1;
    
end
numberOfsamplesAgainstTime=NOSAT;
end
    
    