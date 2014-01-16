%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE OREGONATOR, BELUSOV-ZHABOTINSKII REACTION PROBLEM.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function OREGO
global MYtagSTIFF 

MYtagSTIFF = 1;
ND=3;
LWORK=(31+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);
N = 3;

MASBND = zeros(1,4);
MBND = zeros(1,4);
         

MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;

    
% C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% C     WORK() HOUSES THE FOLLOWING ARRAYS
% C
% C     Y(N,12)   , YHOLD(N,12) , YNHOLD(N,2) , YMAX(N)
% C     ERRORS(N) , SAVE1(N) , SAVE2(N) , SCALE(N) , ARH(N) , PW(MBND(4)*N)
% C     PWCOPY(MBND(4)*N) ,AM(MASBND(4)*N)
% C     IF THE SPARSE OPTION IS NOT BEING USED THEN MBND(4)=N=MASBND(4).
% C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
IWORK_15 = zeros(1,ND);
WORK_1   = 0; %WORK(1)
WORK_2   = 0; %WORK(2)
WORK_3   = zeros(ND,12); %WORK(3).....ETC
WORK_I1  = zeros(ND,12);
WORK_I2  = zeros(ND,2);
WORK_I3  = zeros(1,ND);
WORK_I4  = zeros(1,ND);
WORK_I5  = zeros(1,ND);
WORK_I6  = zeros(1,ND);
WORK_I7  = zeros(1,ND);
WORK_I8  = zeros(1,ND);
WORK_I9  = zeros(1,N*ND);
WORK_I10 = zeros(1,N*ND);
WORK_I11 = zeros(1,MASBND(4)*ND);
      
LOUT = fopen('oreg.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE OREG PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);


%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=360;
         XOUT=30;
         
%...  REQUIRED TOLERANCE

        RTOL=1E-7;
        ATOL=RTOL*1E-6;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
      Y(1)=1;
      Y(2)=2;
      Y(3)=3;

%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=6;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=4E-6;
        MASBND(1) = 0;
%...
%... CALL OF THE SUBROUTINE
%...


for I = 1:12
    MYtag = 0;

    while (MYtag == 0)
          
          [N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,...
              WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,...
              WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,...
              MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] =MEBDF...
              (N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,...
              WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,...
              WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,...
              IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR); 

           if (INDEX==1)
              INDEX = 0;
              
              continue
           else
             
              MYtag =34;
             
              
           end
           
    end
   XOUT = XOUT+30;
Y1(I) = Y(1);
Y2(I) = Y(2);
Y3(I) = Y(3);

end
X = 30;
for t = 1:12
    TRUE = SOLN(t);
    
    
         fprintf(LOUT,'\r\rAPROX Y1 = %E',Y1(t));
         fprintf(LOUT,'\rTRUE    Y1 = %E\r',TRUE(1));
         
         fprintf(LOUT,'\rAPROX Y2 = %E',Y2(t));
         fprintf(LOUT,'\rTRUE  Y2 = %E\r',TRUE(2));
         
         fprintf(LOUT,'\rAPROX Y3 = %E',Y3(t));
         fprintf(LOUT,'\rTRUE  Y3 = %E\r',TRUE(3));
         
         fprintf(LOUT,'\r AT X = %f\r\r',X);
X = X+30;         
end      

      fclose(LOUT)
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TRUE = SOLN(I)
    switch I
      case 1
      TRUE(1)=0.1000661467180497E+01;
      TRUE(2)=0.1512778937348249E+04;
      TRUE(3)=0.1035854312767229E+05;
      case 2
      TRUE(1)=0.1000874625199626E+01;
      TRUE(2)=0.1144336972384497E+04;
      TRUE(3)=0.8372149966624639E+02;
      case 3
      TRUE(1)=0.1001890368438751E+01;
      TRUE(2)=0.5299926232295553E+03;
      TRUE(3)=0.1662279579042420E+01;
      case 4
      TRUE(1)=0.1004118022612645E+01;
      TRUE(2)=0.2438326079910346E+03;
      TRUE(3)=0.1008822224048647E+01;
      case 5
      TRUE(1)=0.1008995416634061E+01;
      TRUE(2)=0.1121664388662539E+03;
      TRUE(3)=0.1007783229065319E+01;
      case 6
      TRUE(1)=0.1019763472537298E+01;
      TRUE(2)=0.5159761322947535E+02;
      TRUE(3)=0.1016985778956374E+01;
      case 7
      TRUE(1)= 0.1043985088527474E+01;
      TRUE(2)=0.2373442027531524E+02;
      TRUE(3)=0.1037691843544522E+01;
      case 8
      TRUE(1)=0.1100849071667922E+01;
      TRUE(2)=0.1091533805469020E+02;
      TRUE(3)=0.1085831969810860E+01;
      case 9
      TRUE(1)=0.1249102130020572E+01;
      TRUE(2)=0.5013945178605446E+01;
      TRUE(3)=0.1208326626237875E+01;
      case 10
      TRUE(1)=0.1779724751937019E+01;
      TRUE(2)=0.2281852385542403E+01;
      TRUE(3)=0.1613754023671725E+01;
      case 11
      TRUE(1)=0.1000889326903503E+01;
      TRUE(2)=0.1125438585746596E+04;
      TRUE(3)=0.1641049483777168E+05;
      case 12
      TRUE(1)=0.1000814870318523E+01;
      TRUE(2)=0.1228178521549889E+04;
      TRUE(3)=0.1320554942846513E+03;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NNN,X,Y,DY,ipar,rpar,ierr] = F(NNN,X,Y,DY,ipar,rpar,ierr)
%disp('now in F')
% C...  RIGHT-HAND SIDE OF OREGON EQUATION
     
      DY(1)=77.27*(Y(2)+Y(1)*(1-8.375E-6*Y(1)-Y(2)));
      DY(2)=(Y(3)-(1+Y(1))*Y(2))/77.27;
      DY(3)=0.161*(Y(1)-Y(3));
  
%disp('letf F')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MAS(n,am,m,IPAR,RPAR,ierr)

%DUMMY ROUTINE

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y,DFY,N,MEBND,ipar,rpar,ierr] = PDERV(X,Y,DFY,N,MEBND,ipar,rpar,ierr)
%disp('now in pderv')
 
  

      DFY(1,1)=77.27*(1-2*8.375E-6*Y(1)-Y(2));
      DFY(1,2)=77.27D0*(1.D0-Y(1));
      DFY(1,3)=0;
      DFY(2,1)=-Y(2)/77.27;
      DFY(2,2)=-(1+Y(1))/77.27;
      DFY(2,3)=1/77.27;
      DFY(3,1)=.161;
      DFY(3,2)=0;
      DFY(3,3)=-.161;


%disp('left pderv')
end



