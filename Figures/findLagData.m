function [T3D,T2D]=findLagData(data,indexarray,nlags,currentTime)

k=10; %default number of data at each time.
timeinterval=2e-7;
A=indexarray;
T=k*round(currentTime/timeinterval)-k+2 %%%%%%%%round any floating point calculation, or the indexes cannot be read. 
N=nlags;
[x,y]=size(data);
target=indexarray(end,end);
[sizeX,sizeY]=size(indexarray);

if N>currentTime/timeinterval
    N=round(currentTime/timeinterval);
end;    


TT3D=zeros(N,sizeY-1,k);
for i=1:N  %%3
     for z=1:sizeY-1      %%2
         Array=zeros(1,k);
       for j=1:k  %%12
           %%2
            
        Array(j)=data(indexarray(z),(T-k*(i-1)+j-1));
             end;
    
    TT3D(i,z,:)=Array;
    end;
    
end;

%%%%%%change the demensions of T3D
T3D=zeros(N,k,sizeY-1);
for i=1:N     
   for j=1:k  
       for z=1:sizeY-1 
       T3D(i,j,z)=TT3D(i,z,j);
       end;
    end;
end;





T2D=zeros(N,k);
for i=1:N
    for j=1:k
        T2D(i,j)=data(target,T-k*(i-1)+j-1);
    end;
end;


