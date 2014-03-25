%T3D and T2D are output of 3D and 2D data.     
%data can be either the raw data or the deviation data. 
%nlags is the number of lagged data. 
%currentTime is the time T when we want to obeserve.
%indexarray is the array of indexes that need to be explored, and Array[end] is the node we are examining and the rest are the nodes surrounding this node.  i.e.  Array=[3,5,4], the output will have a T3D of 10 at each time, nlags lagged data for data row 3 and row 5, and T2D of 10 at each time, nlags lagged data for row 4.
function [T3D,T2D]=findLagData(data,indexarray,nlags,currentTime, k)
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


