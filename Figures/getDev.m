
function Deviation=getDev(theoticalVal, experimentalVal)
    

tV=theoticalVal;
eV=rrt2mat(experimentalVal);
[x,y]=size(eV);

time_interval=2e-7;

timemax=max(eV(7,:))
[x,y]=size(eV)
Dev=zeros(y);
for i=1:y
    i;
    t=eV(end,i);
    index=floor(t/time_interval)+1;
    if index>2037
        index=2037;
    end;
    index;
    newTimeDev=eV(3:6,i)-tV(3:6,index);
    
    if (i==1)
          Dev=newTimeDev;
    else 
          Dev=[Dev,newTimeDev];
    end
end
Dev;
nodecopy=eV(1:2,:);
timecopy=eV(end,:);
[a,b]=size(nodecopy)
[c,d]=size(timecopy)
[e,f]=size(Dev)

Deviation=[nodecopy;Dev;timecopy];


end