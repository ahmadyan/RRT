
function Deviation=getDev(theoticalVal, eV)
tV=theoticalVal;
%eV=rrt2mat(experimentalVal);   %Adel: I've changed this line, now the
%input is eV, not the text file
[x,y]=size(eV);

time_interval=2e-7;

timemax=max(eV(7,:))
[x,y]=size(eV)
Dev=zeros(y);
for i=1:y
    i;
    t=eV(end,i);
    index=floor(t/time_interval)+1;
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