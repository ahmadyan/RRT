%[X,Y]=meshgrid(-0.4:0.1:1.6, -0.4:0.1:1.6);
xmin=-0.4; 
xmax=1.6 ;
ymin=-0.4;
ymax=1.6;
resolution=10;
sim_time=0.7 ;

Z=zeros(20,20);

for i=1:size(sample,1),
    x=sample(i,1);
    y=sample(i,2);
   
    
    
    Z(4+floor(x*10), 4+floor(y*10)) = Z(4+floor(x*10), 4+floor(y*10))+1;
end