%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE PENDULUM PROBLEM: INDEX = 3.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function PENDUL3
global MYtagSTIFF
MYtagSTIFF = 1;
ND=5;
N=ND;
LWORK=(41+2*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);

MBND = zeros(1,4);
TRUE  = zeros(1,ND);
ERROR = zeros(1,ND);
MASBND = zeros(1,4);
MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;
        MASBND(1)=1;
        MASBND(2)=N;
        MASBND(3)=1;
        MASBND(4)=N;
        MBND(4)=N;

%...  DIMENSION OF THE SYSTEM
      
% C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% C     WORK() HOUSES THE FOLLOWING ARRAYS
% C
% C     Y(N,12)   , YHOLD(N,12) , YNHOLD(N,2) , YMAX(N)
% C     ERRORS(N) , SAVE1(N) , SAVE2(N) , SCALE(N) , ARH(N) , PW(MBND(4)*N)
% C     PWCOPY(MBND(4)*N) ,AM(MASBND(4)*N)
% C     IF THE SPARSE OPTION IS NOT BEING USED THEN MBND(4)=N=MASBND(4).
% C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
IWORK_15 = zeros(1,N);
WORK_1   = 0; %WORK(1)
WORK_2   = 0; %WORK(2)
WORK_3   = zeros(N,12); %WORK(3).....ETC
WORK_I1  = zeros(N,12);
WORK_I2  = zeros(N,2);
WORK_I3  = zeros(1,N);
WORK_I4  = zeros(1,N);
WORK_I5  = zeros(1,N);
WORK_I6  = zeros(1,N);
WORK_I7  = zeros(1,N);
WORK_I8  = zeros(1,N);
WORK_I9  = zeros(1,MBND(4)*N);
WORK_I10 = zeros(1,MBND(4)*N);
WORK_I11 = zeros(1,MASBND(4)*N);




      
LOUT = fopen('pendul3.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE PENDUL3 PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);


%... ENDPOINT OF INTEGRATION
    
         XEND=1;
         XOUT=XEND;
         X=0;
%...  REQUIRED TOLERANCE

        RTOL=10.0E-8;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
          Y(1)=1;
          Y(2)=0;
          Y(3)=0;
          Y(4)=1;
          Y(5)=1;
%...  SET DEFAULT VALUES
        MF=22;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-8;
        fprintf(LOUT,'Hstart = %E',H);
        IWORK(1)=N-3;
        IWORK(2)=2;
        IWORK(3)=1;
     

        

%...
%... CALL OF THE SUBROUTINE
%...

MYtag = 0

while (MYtag == 0)
             
          [N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] = MEBDF(N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR);  
 

           if (INDEX==1)
              INDEX = 0;
              
              continue
           else
             
              MYtag =34;
              
           end
           
end


      TRUE(1)=0.867348640604551746;
      TRUE(2)=0.497701050476373918;
      TRUE(3)=-0.337480180329866458E-1;
      TRUE(4)=0.588130114660715686E-1;
      TRUE(5)=-0.493103151528274442;


      
for t = 1:N
         fprintf(LOUT,'\r\r APROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE Y%g = %E\r\r',t,TRUE(t));
end      

      fclose(LOUT)
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,X,Y,DF,IPAR,RPAR,ierr] = F(N,X,Y,DF,IPAR,RPAR,ierr)

%disp('now in F')

      DF(1)=Y(3);
      DF(2)=Y(4);
      DF(3)=-Y(1)*Y(5);
      DF(4)=-Y(2)*Y(5)-1;
      DF(5)=Y(1)*Y(1)+Y(2)*Y(2)-1;



%disp('letf F')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n,am,m,IPAR,RPAR,ierr] = MAS(n,am,m,IPAR,RPAR,ierr)
%disp('IN MAS')
      am(1,1)=1;
      am(1,2)=0;
      am(2,1)=0;
      am(2,2)=1;
      am(3,3)=1;
      am(4,4)=1;
      am(5,5)=0;


%disp('LEFT MAS')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y,DFY,N,MEBAND,IPAR,RPAR,ierr] =  PDERV(X,Y,DFY,N,MEBAND,IPAR,RPAR,ierr)
%disp('now in pderv')

      DFY(1,3)=1;
      DFY(2,4)=1;
      DFY(3,1)=-Y(5);
      DFY(3,5)=-Y(1);
      DFY(4,2)=-Y(5);
      DFY(4,5)=-Y(2);
      DFY(5,1)=2*Y(1);
      DFY(5,2)=2*Y(2);

%disp('left pderv')
end
