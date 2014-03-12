function number=countExceptions(experimentalVal)
    

eV=rrt2mat(experimentalVal);

time_interval=2e-7;

timemax=max(eV(7,:))
[x,y]=size(eV);
number=0;
for i=1:y
    t=eV(end,i);
    
    if (t<6e-5) &&(t>4e-5)
        number=number+1;
    
    end
end


end