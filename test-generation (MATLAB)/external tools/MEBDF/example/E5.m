%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE E5 PROBLEM.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function E5
global MYtagSTIFF 

MYtagSTIFF = 1;
ND=4;
LWORK=(34+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);
N = 4;

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
      
LOUT = fopen('E5.txt','w');


fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE E5 PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);



%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=10000000000000;
         XOUT=10;
         
%...  REQUIRED TOLERANCE

        RTOL=1E-7;
        ATOL=RTOL*1E-24;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        Y(1)=1.76E-3;
        Y(2)=0;
        Y(3)=0;
        Y(4)=0;
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=6;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-8;
        fprintf(LOUT,'Hstart = %E',H);
        MASBND(1) = 0;
%...
%... CALL OF THE SUBROUTINE
%...


for I = 1:7
    MYtag = 0;

    while (MYtag == 0)
          
         
           [N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] = MEBDF(N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR);  


           if (INDEX==1)
              INDEX = 0;
              
              continue
           else
             
              MYtag =34;
             
              
           end
           
    end
   XOUT = XOUT*100;
Y1(I) = Y(1);
Y2(I) = Y(2);
Y3(I) = Y(3);
Y4(I) = Y(4);


end
X = 10;
for t = 1:7
    TRUE = SOLN(t);
    
    
         fprintf(LOUT,'\r\rAPROX Y1 = %E',Y1(t));
         fprintf(LOUT,'\rTRUE  Y1 = %E\r',TRUE(1));
         
         fprintf(LOUT,'\rAPROX Y2 = %E',Y2(t));
         fprintf(LOUT,'\rTRUE  Y2 = %E\r',TRUE(2));
         
         fprintf(LOUT,'\rAPROX Y3 = %E',Y3(t));
         fprintf(LOUT,'\rTRUE  Y3 = %E\r',TRUE(3));
         
         fprintf(LOUT,'\rAPROX Y4 = %E',Y4(t));
         fprintf(LOUT,'\rTRUE  Y4 = %E\r',TRUE(4));
         
         fprintf(LOUT,'\r AT X = %f\r\r',X);
X = X*100;         
end      

      fclose(LOUT)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TRUE = SOLN(I)
    switch I
      case 1
      TRUE(1)=1.7599259497677897058E-003;
         TRUE(2)=1.3846281519376516449E-011;
         TRUE(3)=7.6370038530073911180E-013;
         TRUE(4)=1.3082581134075777338E-011;
      case 2
      TRUE(1)=1.6180769999072942552E-003;
         TRUE(2)=1.3822370304983735443E-010;
         TRUE(3)=8.2515735006838336088E-012;
         TRUE(4)=1.2997212954915352082E-010;
      case 3
      TRUE(1)=7.4813208224292220114E-006;
         TRUE(2)=2.3734781561205975019E-012;
         TRUE(3)=2.2123586689581663654E-012;
         TRUE(4)=1.6111948716243113653E-013;

      case 4
      TRUE(1)=4.7150333630401632232E-010;
         TRUE(2)=1.8188895860807021729E-014;
         TRUE(3)=1.8188812376786725407E-014;
         TRUE(4)=8.3484020296321693074E-020;

      case 5
      TRUE(1)=3.1317148329356996037E-014;
         TRUE(2)=1.4840957952870064294E-016;
         TRUE(3)=1.4840957948345691466E-016;
         TRUE(4)=4.5243728279782625194E-026;

      case 6
       TRUE(1)=3.8139035189787091771E-049;
         TRUE(2)=1.0192582567660293322E-020;
         TRUE(3)=1.0192582567660293322E-020;
         TRUE(4)=3.7844935507486221171E-065;
      case 7
       TRUE(1)=0.0000000000000000000E-000;
         TRUE(2)=8.8612334976263783420E-023;
         TRUE(3)=8.8612334976263783421E-023;
         TRUE(4)=0.0000000000000000000E-000;

    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,X,Y,DF,IPAR,RPAR,IER] = F(N,X,Y,DF,IPAR,RPAR,IER)
%disp('now in F')

% C...  RIGHT-HAND SIDE OF E5 EQUATION
     
      PROD1=7.89E-10*Y(1);
      PROD2=1.1E7*Y(1)*Y(3);
      PROD3=1.13E9*Y(2)*Y(3);
      PROD4=1.13E3*Y(4);
      DF(1)=-PROD1-PROD2;
      DF(2)=PROD1-PROD3;
      DF(4)=PROD2-PROD4;
      DF(3)=DF(2)-DF(4);
  
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
 
 
% C --- JACOBIAN OF E5 EQUATION
    
      A=7.89E-10;
      B=1.1E7;
      CM=1.13E9;
      C=1.13E3;
      DFY(1,1)=-A-B*Y(3);
      DFY(1,2)=0;
      DFY(1,3)=-B*Y(1);
      DFY(1,4)=0;
      DFY(2,1)=A;
      DFY(2,2)=-CM*Y(3);
      DFY(2,3)=-CM*Y(2);
      DFY(2,4)=0;
      DFY(3,1)=A-B*Y(3);
      DFY(3,2)=-CM*Y(3);
      DFY(3,3)=-B*Y(1)-CM*Y(2);
      DFY(3,4)=C;
      DFY(4,1)=B*Y(3);
      DFY(4,2)=0;
      DFY(4,3)=B*Y(1);
      DFY(4,4)=-C;

%disp('left pderv')
end




