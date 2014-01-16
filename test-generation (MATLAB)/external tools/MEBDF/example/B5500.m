%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE B5500 EQUATION.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function B5500
global MYtagSTIFF 

MYtagSTIFF = 1;
ND=6;
LWORK=(32+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);
N = 6;

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
      
LOUT = fopen('B5500.txt','w');


%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=20;
         XOUT=20;
         
%...  REQUIRED TOLERANCE

        RTOL=1E-10;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        Y = ones(1,6);
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-6;
        MASBND(1) = 0;
%...
%... CALL OF THE SUBROUTINE
%...



MYtag = 0;
MYcase = 220;
    while (MYtag == 0)
        
        
        switch MYcase
            
            case 220
          
         [N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,...
              WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,...
              WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,...
              MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] =MEBDF...
              (N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,...
              WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,...
              WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,...
              IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR);
          
           if (X > XOUT)
               MYcase = 995;
               continue
           end
           
           if (INDEX == 0)
               MYcase = 995;
               continue
           end
           
           INDEX = 3;
           
           MYcase = 220;
           continue
           
            case 995
        
                XOUT = X;
                TRUE(1)=exp(-10*XOUT)*(cos(500*XOUT)+sin(500*XOUT));
                TRUE(2)=exp(-10*XOUT)*(cos(500*XOUT)-sin(500*XOUT));
                TRUE(3)=exp(-4*XOUT);
                TRUE(4)=exp(-XOUT);
                TRUE(5)=exp(-0.5*XOUT);
                TRUE(6)=exp(-0.1*XOUT);  
                
                MYtag = 34;
                
                continue
                
        end
    end
                

 fprintf(LOUT,'\r\rX = %g',X);
fprintf(LOUT,'\r\rXOUT = %g',XOUT);


for t = 1:6
    
         fprintf(LOUT,'\r\rAPROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE  Y%g = %E\r\r',t,TRUE(t));

        
end      

      fclose(LOUT)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,X,Y,DF,IPAR,RPAR,IER] = F(N,X,Y,DF,IPAR,RPAR,IER)
%disp('now in F')

      DF(1)=-10*Y(1)+500*Y(2);
      DF(2)=-500*Y(1)-10*Y(2);
      DF(3)=-4*Y(3);
      DF(4)=-Y(4);
      DF(5)=-0.5*Y(5);
      DF(6)=-0.1*Y(6);
      
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
 
    
      DFY(1,1)=-10;
      DFY(1,2)=500;
      DFY(2,1)=-500;
      DFY(2,2)=-10;
      DFY(3,3)=-4;
      DFY(4,4)=-1;
      DFY(5,5)=-0.5;
      DFY(6,6)=-0.1;

%disp('left pderv')
end
