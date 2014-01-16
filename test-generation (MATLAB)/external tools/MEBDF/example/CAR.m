%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE CAR AXIS PROBLEM: INDEX = 3.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function CAR
global MYtagSTIFF 

%FOR USE IN RESMBS%


MYtagSTIFF = 1;

ND=10;
NEQN = 10;
N =10;

LWORK=(41+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);


%MASBND = zeros(1,4);
%MBND = zeros(1,4);
            
         MASBND(1) = 1;
         MASBND(2) = 0;
         MASBND(3) = 0;
         MASBND(4) = 1;
         MBND(1) = NEQN;
         MBND(2) = NEQN;
         MBND(3) = NEQN;
         MBND(4) = NEQN;
         

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
WORK_I9  = zeros(1,MBND(4)*ND);
WORK_I10 = zeros(1,MBND(4)*ND);
WORK_I11 = zeros(1,MASBND(4)*ND);
      
LOUT = fopen('car.txt','w');

fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE CAR PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);
%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=3;
         XOUT=XEND;
         
%...  REQUIRED TOLERANCE

        %RTOL=.316227766016837972E-06;
        RTOL=0.100000000000000009E-07;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        [N,X,Y,YPRIME] = INIT(N,X,Y,YPRIME);
%...  SET DEFAULT VALUES
        MF=22;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=200000;
        H=1E-8;
         fprintf(LOUT,'\rHstart = %E',H);
         IWORK(1)=4;
         IWORK(2)=4;
         IWORK(3)=2;    
      

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

           
           if (INDEX == 1)
               INDEX = 0;
               MYcase = 220;
               continue
           end
           
          MYtag = 34;
          continue
                
        end
    end
    
     TRUE = SOLN(1);


 fprintf(LOUT,'\r\rX = %g',X);
fprintf(LOUT,'\r\rXOUT = %g',XOUT);


for t = 1:N
    
         fprintf(LOUT,'\r\rAPROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE  Y%g = %E\r\r',t,TRUE(t));

        
end      

      fclose(LOUT)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  YEND = SOLN(dummy)
     
      YEND(  1) =  0.4934557842755629d-001;
      YEND(  2) =  0.4969894602303324d+000;
      YEND(  3) =  0.1041742524885400d+001;
      YEND(  4) =  0.3739110272652214d+000;
      YEND(  5) = -0.7705836840321485d-001;
      YEND(  6) =  0.7446866596327776d-002;
      YEND(  7) =  0.1755681574942899d-001;
      YEND(  8) =  0.7703410437794031d+000;
      YEND(  9) = -0.4736886750784630d-002;
      YEND( 10) = -0.1104680411345730d-002;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,T,y,YPRIME] = INIT(NEQN,T,y,YPRIME)
            

      
      M   = 10d0;
      eps = 1d-2;
      L   = 1d0;
      L0  = 0.5d0;

      xr   = L;
      xl   = 0;
      yr   = L0;
      yl   = yr;
      xra  = -L0/L;
      xla  = xra;
      yra  = 0d0;
      yla  = 0d0;
      lam1 = 0d0;
      lam2 = 0d0;

      y(1)  =  xl;
      y(2)  =  yl;
      y(3)  =  xr;
      y(4)  =  yr;
      y(5)  =  xla;
      y(6)  =  yla;
      y(7)  =  xra;
      y(8)  =  yra;
      y(9)  =  lam1;
      y(10) =  lam2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neqn,t,y,dy,ipar,rpar,ierr] = F(neqn,t,y,dy,ipar,rpar,ierr)
%disp('now in F')
   


      eps = 1d-2;
      M   = 10d0;
      L   = 1d0;
      L0  = 0.5d0;
      r   = 0.1d0;
      w   = 10d0;
      g   = 1d0;
      yb  = r*sin(w*t);
      xb  = sqrt(L*L-yb*yb);

      for i=1:4
         dy(i) = y(i+4);
      end
      

      xl   = y(1);
      yl   = y(2);
      xr   = y(3);
      yr   = y(4);
      lam1 = y(9);
      lam2 = y(10);

      Ll = sqrt(xl^2+yl^2);
      Lr = sqrt((xr-xb)^2+(yr-yb)^2);

      dy(5) =(L0-Ll)*xl/Ll +lam1*xb+2d0*lam2*(xl-xr);
      dy(6) =(L0-Ll)*yl/Ll +lam1*yb+2d0*lam2*(yl-yr)-M*eps*eps*g/2d0;
      dy(7) =(L0-Lr)*(xr-xb)/Lr -2d0*lam2*(xl-xr);
      dy(8) =(L0-Lr)*(yr-yb)/Lr -2d0*lam2*(yl-yr)-M*eps*eps*g/2d0;

      dy(9) = xb*xl+yb*yl;
      dy(10)= (xl-xr)^2+(yl-yr)^2-L*L;
      
     

%disp('letf F')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neqn,am,ldim,ipar,rpar,ier] = MAS(neqn,am,ldim,ipar,rpar,ier)

%disp('IN MAS')
 
      M   = 10d0;
      eps = 1d-2;
      
      %for i=1:4
      %   am(1,i) = 1d0;
      %end
      am(1,1:4) = 1;
      
      %for i=5:8
      %   am(1,i) = M*eps*eps/2d0;
      %end
      am(1,5:8) = M*eps*eps/2d0;
      
      %for i=9:10
      %   am(1,i) = 0d0;
      %end
      
      am(1,9:10) = 0;


%disp('LEFT MAS')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,y,jac,neqn,mn,ipar,rpar,ierr]...
    = PDERV(t,y,jac,neqn,mn,ipar,rpar,ierr)
%disp('now in pderv')

%DUMMY ROUTINE
     
%disp('left pderv')
end