%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE RING MODULATOR PROBLEM WITH cs=1e-9.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function [xf] = ring(x0, t0, dt)
global MYtagSTIFF 

MYtagSTIFF = 1;
ND=15;
LWORK=(41+2*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);
N = 15;

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
WORK_I9  = zeros(1,ND*ND);
WORK_I10 = zeros(1,ND*ND);
WORK_I11 = zeros(1,MASBND(4)*ND);
      
%LOUT = fopen('easyring.txt','w');
%fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE EASY RING PROBLEM USING MEBDF \r');
%fprintf(LOUT,'\r Number of Equations : %g \r',ND);


%... ENDPOINT OF INTEGRATION
         X=t0;
         XEND=t0+dt;
         
         XOUT=XEND;
         
%...  REQUIRED TOLERANCE

        RTOL=1E-7;
        ATOL=RTOL;
        ITOL=2;
        %fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        %fprintf(LOUT,'\r RTOL = %E',RTOL);
        %fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        Y = x0; %zeros(1,N);
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=300000;
        H=1E-8;
        %fprintf(LOUT,'Hstart = %E',H);
        MASBND(1) = 0;
%...
%... CALL OF THE SUBROUTINE
%...



MYtag = 0;
MYcase = 220;
    while (MYtag == 0)
        
        
        switch MYcase
            
            case 220
          
           [N,X,H,Y,XOUT,XEND,MF,INDEX,LWORK,WORK_1,WORK_2,WORK_3,...
              WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,...
              WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,...
              MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] =MEBDF...
              (N,X,H,Y,XOUT,XEND,MF,INDEX,LWORK,WORK_1,WORK_2,...
              WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,...
              WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,...
              IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR);  

           
           if (INDEX == 1)
               INDEX = 0;
               MYcase = 220;
               continue
           end
           
          MYtag = 34;
          continue
                
        end
    end
    
%      TRUE(1)=-0.170799033E-1;
%     TRUE(2)=-0.6660979d-2;
%    TRUE(3)=+0.27531919;
%      TRUE(4)=-0.39115732;
%      TRUE(5)=-0.38851731;
%      TRUE(6)=+0.27795920;
%      TRUE(7)=+0.1114600281;
%      TRUE(8)=+0.2979129627E-6;
%      TRUE(9)=-0.31427403E-7;
%      TRUE(10)=+0.70165883E-3;
%      TRUE(11)=+0.85207538E-3;
%      TRUE(12)=-0.77741454E-3;
%      TRUE(13)=-0.77631967E-3;
%      TRUE(14)=+0.78439426E-4;
%      TRUE(15)=+0.252322784E-4;


 %fprintf(LOUT,'\r\rX = %g',X);
%fprintf(LOUT,'\r\rXOUT = %g',XOUT);

xf=Y;

%for t = 1:N
    
 %        fprintf(LOUT,'\r\rAPROX Y%g = %E',t,Y(t));
  %       fprintf(LOUT,'\rTRUE  Y%g = %E\r\r',t,TRUE(t));

        
%end      

    %  fclose(LOUT)
    
end
