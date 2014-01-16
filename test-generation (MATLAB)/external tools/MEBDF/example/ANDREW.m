%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE ANDREW SQUEEZING MECHANISM: INDEX = 3.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function ANDREW
global MYtagSTIFF 

%FOR USE IN RESMBS%

MYtagSTIFF = 1;

ND=27;
NEQN = 27;
N = 27;

LWORK=(38+3*ND)*ND+3;
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
      
LOUT = fopen('andrew.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE ANDREW PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);


%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=3E-2;
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
        MF=21;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-8;
        fprintf(LOUT,'\rHstart = %E',H);
         IWORK(1)=7;
         IWORK(2)=7;
         IWORK(3)=13  ;    
      

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

function  Y = SOLN(dummy)
   
  
% C
% c from CWI 
% c computed at Cray C90, using Cray double precision:
% c  Solving Andrews' squeezing mechanism using PSIDE
% c
% c  User input:
% c
% c  give relative error tolerance: 1d-14
% c  give absolute error tolerance: 1d-14
% c
% c  Integration characteristics:
% c
% c     number of integration steps        1112
% c     number of accepted steps           1022
% c     number of f evaluations           26000
% c     number of Jacobian evaluations      173
% c     number of LU decompositions        2056
% c
% c  CPU-time used:                          56.74 sec


      Y(  1) =  0.1581077119629904D+2;
      Y(  2) = -0.1575637105984298D+2;
      Y(  3) =  0.4082224013073101D-1;
      Y(  4) = -0.5347301163226948D+0;
      Y(  5) =  0.5244099658805304D+0;
      Y(  6) =  0.5347301163226948D+0;
      Y(  7) =  0.1048080741042263D+1;
      Y(  8) =  0.1139920302151208D+4;
      Y(  9) = -0.1424379294994111D+4;
      Y( 10) =  0.1103291221937134D+2;
      Y( 11) =  0.1929337464421385D+2;
      Y( 12) =  0.5735699284790808D+0;
      Y( 13) = -0.1929337464421385D+2;
      Y( 14) =  0.3231791658026955D+0;
      Y( 15) = -0.2463176316945196D+5;
      Y( 16) =  0.5185037701610329D+5;
      Y( 17) =  0.3241025686413781D+6;
      Y( 18) =  0.5667493645176213D+6;
      Y( 19) =  0.1674362929479361D+5;
      Y( 20) = -0.5667493645176222D+6;
      Y( 21) =  0.9826520791458422D+4;
      Y( 22) =  0.1991753333731910D+3;
      Y( 23) = -0.2975531228015052D+2;
      Y( 24) =  0.2306654119098399D+2;
      Y( 25) =  0.3145271365475927D+2;
      Y( 26) =  0.2264249232082739D+2;
      Y( 27) =  0.1161740700019673D+2;
 

      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,T,Y,YPRIME] = INIT(NEQN,T,Y,YPRIME)
            

      Y(1 ) = -0.0617138900142764496358948458001D0;
      Y(2 ) =  0D0;
      Y(3 ) =  0.455279819163070380255912382449D0;
      Y(4 ) =  0.222668390165885884674473185609D0;
      Y(5 ) =  0.487364979543842550225598953530D0;
      Y(6 ) = -0.222668390165885884674473185609D0;
      Y(7 ) =  1.23054744454982119249735015568D0;
      Y(8 ) =  0D0;
      Y(9 ) =  0D0;
      Y(10) =  0D0;
      Y(11) =  0D0;
      Y(12) =  0D0;
      Y(13) =  0D0;
      Y(14) =  0D0;
      Y(15) =  14222.4439199541138705911625887D0;
      Y(16) = -10666.8329399655854029433719415D0;
      Y(17) =  0D0;
      Y(18) =  0D0;
      Y(19) =  0D0;
      Y(20) =  0D0;
      Y(21) =  0D0;
      Y(22) =  98.5668703962410896057654982170D0;
      Y(23) = -6.12268834425566265503114393122D0;
      Y(24) =  0D0;
      Y(25) =  0D0;
      Y(26) =  0D0;
      Y(27) =  0D0;

      %for k=1:14
      %   yprime(k) = y(k+7);
      %end
      
      %for k=15:27
      %   yprime(k) = 0d0;
      %end
YPIRME(1:14) = Y(8:21);
YPRIME(15:27) = 0d0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neqn,t,y,dy,ipar,rpar,ierr] = F(neqn,t,y,dy,ipar,rpar,ierr)
%disp('now in F')
   
 
m1=.04325d0;m2=.00365d0;m3=.02373d0;m4=.00706d0;
m5=.07050d0;m6=.00706d0;m7=.05498d0;
xa=-.06934d0;ya=-.00227d0;
xb=-0.03635d0;yb=.03273d0;
xc=.014d0;yc=.072d0;c0=4530d0;
i1=2.194d-6;i2=4.410d-7;i3=5.255d-6;i4=5.667d-7;
i5=1.169d-5;i6=5.667d-7;i7=1.912d-5;
d=28d-3;da=115d-4;e=2d-2;ea=1421d-5;
rr=7d-3;ra=92d-5;l0=7785d-5;
ss=35d-3;sa=1874d-5;sb=1043d-5;sc=18d-3;sd=2d-2;
ta=2308d-5;tb=916d-5;u=4d-2;ua=1228d-5;ub=449d-5;
zf=2d-2;zt=4d-2;fa=1421d-5;mom=33d-3;
      
      sibe = sin(y(1));
      sith = sin(y(2));
      siga = sin(y(3));
      siph = sin(y(4));
      side = sin(y(5));
      siom = sin(y(6));
      siep = sin(y(7));
      
      cobe = cos(y(1));
      coth = cos(y(2));
      coga = cos(y(3));
      coph = cos(y(4));
      code = cos(y(5));
      coom = cos(y(6));
      coep = cos(y(7));
      
      sibeth = sin(y(1)+y(2));
      siphde = sin(y(4)+y(5));
      siomep = sin(y(6)+y(7));
      
      cobeth = cos(y(1)+y(2));
      cophde = cos(y(4)+y(5));
      coomep = cos(y(6)+y(7));
 
      bep = y(8);
      thp = y(9);
      php = y(11);
      dep = y(12);
      omp = y(13);
      epp = y(14);
 
      for j = 1:7
         for i = 1:7
            m(i,j) = 0d0;
         end
      end
 
      m(1,1) = m1*ra^2 + m2*(rr^2-2*da*rr*coth+da^2) + i1 + i2;
      m(2,1) = m2*(da^2-da*rr*coth) + i2;
      m(2,2) = m2*da^2 + i2;
      m(3,3) = m3*(sa^2+sb^2) + i3;
      m(4,4) = m4*(e-ea)^2 + i4;
      m(5,4) = m4*((e-ea)^2+zt*(e-ea)*siph) + i4;
      m(5,5) = m4*(zt^2+2*zt*(e-ea)*siph+(e-ea)^2) + m5*(ta^2+tb^2)+ i4 + i5;
      m(6,6) = m6*(zf-fa)^2 + i6;
      m(7,6) = m6*((zf-fa)^2-u*(zf-fa)*siom) + i6;
      m(7,7) = m6*((zf-fa)^2-2*u*(zf-fa)*siom+u^2) + m7*(ua^2+ub^2)+ i6 + i7;
      
      for j=2:7
         for i=1:j-1
            m(i,j) = m(j,i);
         end
      end
      
      xd = sd*coga + sc*siga + xb;
      yd = sd*siga - sc*coga + yb;
      lang  = sqrt ((xd-xc)^2 + (yd-yc)^2);
      force = - c0 * (lang - l0)/lang;
      fx = force * (xd-xc);
      fy = force * (yd-yc);
      ff(1) = mom - m2*da*rr*thp*(thp+2*bep)*sith;
      ff(2) = m2*da*rr*bep^2*sith;
      ff(3) = fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga);
      ff(4) = m4*zt*(e-ea)*dep^2*coph;
      ff(5) = - m4*zt*(e-ea)*php*(php+2*dep)*coph;
      ff(6) = - m6*u*(zf-fa)*epp^2*coom;
      ff(7) = m6*u*(zf-fa)*omp*(omp+2*epp)*coom;
 
      for j=1:7
         for i=1:6
            gp(i,j) = 0d0;
         end
      end
      
      gp(1,1) = - rr*sibe + d*sibeth;
      gp(1,2) = d*sibeth;
      gp(1,3) = - ss*coga;
      gp(2,1) = rr*cobe - d*cobeth;
      gp(2,2) = - d*cobeth;
      gp(2,3) = - ss*siga;
      gp(3,1) = - rr*sibe + d*sibeth;
      gp(3,2) = d*sibeth;
      gp(3,4) = - e*cophde;
      gp(3,5) = - e*cophde + zt*side;
      gp(4,1) = rr*cobe - d*cobeth;
      gp(4,2) = - d*cobeth;
      gp(4,4) = - e*siphde;
      gp(4,5) = - e*siphde - zt*code;
      gp(5,1) = - rr*sibe + d*sibeth;
      gp(5,2) = d*sibeth;
      gp(5,6) = zf*siomep;
      gp(5,7) = zf*siomep - u*coep;
      gp(6,1) = rr*cobe - d*cobeth;
      gp(6,2) = - d*cobeth;
      gp(6,6) = - zf*coomep;
      gp(6,7) = - zf*coomep - u*siep;

      g(1) = rr*cobe - d*cobeth - ss*siga - xb;
      g(2) = rr*sibe - d*sibeth + ss*coga - yb;
      g(3) = rr*cobe - d*cobeth - e*siphde - zt*code - xa;
      g(4) = rr*sibe - d*sibeth + e*cophde - zt*side - ya;
      g(5) = rr*cobe - d*cobeth - zf*coomep - u*siep - xa;
      g(6) = rr*sibe - d*sibeth - zf*siomep + u*coep - ya;

      for i=1:14
         dy(i) = y(i+7);
      end

      for i=15:21
         dy(i) = -ff(i-14);
         for j=1:7
            dy(i) = dy(i)+m(i-14,j)*y(j+14);
         end
         
         for j=1:6
            dy(i) = dy(i)+gp(j,i-14)*y(j+21);
         end
      end
         
      for i=22:27
         dy(i) = g(i-21);
      end


%disp('letf F')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,AM,LDIM,IPAR,RPAR,IERR] = MAS(NEQN,AM,LDIM,IPAR,RPAR,IERR)
%disp('IN MAS')
    
%       DO I=1,NEQN
%          AM(1,I) = 0.D0
%       ENDDO
%       DO J=1,14
%          AM(1,J) = 1.0D0        
%       ENDDO          
AM(1,1:NEQN) = 0;
AM(1,1:14) = 1;

     
%disp('LEFT MAS')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,y,jac,neqn,mn,ipar,rpar,ierr]...
    = PDERV(t,y,jac,neqn,mn,ipar,rpar,ierr)
%disp('now in pderv')

    
% c-----------------------------------------------------------------------
% c     the Jacobian computed here is an approximation, see p. 540 of
% c     Hairer & Wanner `solving ordinary differential equations II'
% c-----------------------------------------------------------------------     
  
m1=.04325d0;m2=.00365d0;m3=.02373d0;m4=.00706d0;
m5=.07050d0;m6=.00706d0;m7=.05498d0;
i1=2.194d-6;i2=4.410d-7;i3=5.255d-6;i4=5.667d-7;
i5=1.169d-5;i6=5.667d-7;i7=1.912d-5;
d=28d-3;da=115d-4;e=2d-2;ea=1421d-5;
rr=7d-3;ra=92d-5;
ss=35d-3;sa=1874d-5;sb=1043d-5;
ta=2308d-5;tb=916d-5;u=4d-2;ua=1228d-5;ub=449d-5;
zf=2d-2;zt=4d-2;fa=1421d-5;

      sibe = sin(y(1));
      siga = sin(y(3));
      siph = sin(y(4));
      side = sin(y(5));;
      siom = sin(y(6));
      siep = sin(y(7));
 
      cobe = cos(y(1));
      coth = cos(y(2));
      coga = cos(y(3));
      code = cos(y(5));
      coep = cos(y(7));
 
      sibeth = sin(y(1)+y(2));
      siphde = sin(y(4)+y(5));
      siomep = sin(y(6)+y(7));
 
      cobeth = cos(y(1)+y(2));
      cophde = cos(y(4)+y(5));
      coomep = cos(y(6)+y(7));
 
      for j = 1:7
         for i = 1:7
            m(i,j) = 0d0;
         end
      end
 
      m(1,1) = m1*ra^2 + m2*(rr^2-2*da*rr*coth+da^2) + i1 +i2;
      m(2,1) = m2*(da^2-da*rr*coth) + i2;
      m(2,2) = m2*da^2 + i2;
      m(3,3) = m3*(sa^2+sb^2) + i3;
      m(4,4) = m4*(e-ea)^2 + i4;
      m(5,4) = m4*((e-ea)^2+zt*(e-ea)*siph) + i4;
      m(5,5) = m4*(zt^2+2*zt*(e-ea)*siph+(e-ea)^2) + m5*(ta^2+tb^2)+ i4 + i5;
      m(6,6) = m6*(zf-fa)^2 + i6;
      m(7,6) = m6*((zf-fa)^2-u*(zf-fa)*siom) + i6;
      m(7,7) = m6*((zf-fa)^2-2*u*(zf-fa)*siom+u^2) + m7*(ua^2+ub^2)+ i6 + i7;

      for j=2:7
        for i=1:j-1
            m(i,j) = m(j,i);
        end
      end
 
      for j=1:7
         for i=1:6
            gp(i,j) = 0d0;
         end
      end
      
      gp(1,1) = - rr*sibe + d*sibeth;
      gp(1,2) = d*sibeth;
      gp(1,3) = - ss*coga;
      gp(2,1) = rr*cobe - d*cobeth;
      gp(2,2) = - d*cobeth;
      gp(2,3) = - ss*siga;
      gp(3,1) = - rr*sibe + d*sibeth;
      gp(3,2) = d*sibeth;
      gp(3,4) = - e*cophde;
      gp(3,5) = - e*cophde + zt*side;
      gp(4,1) = rr*cobe - d*cobeth;
      gp(4,2) = - d*cobeth;
      gp(4,4) = - e*siphde;
      gp(4,5) = - e*siphde - zt*code;
      gp(5,1) = - rr*sibe + d*sibeth;
      gp(5,2) = d*sibeth;
      gp(5,6) = zf*siomep;
      gp(5,7) = zf*siomep - u*coep;
      gp(6,1) = rr*cobe - d*cobeth;
      gp(6,2) = - d*cobeth;
      gp(6,6) = - zf*coomep;
      gp(6,7) = - zf*coomep - u*siep;

      for j=1:neqn
         for i=1:neqn
            jac(i,j) = 0d0;
         end
      end
      
      for i=1:14
         jac(i,i+7) = 1d0;
      end
      for i=1:7
         for j=1:7
            jac(14+j,14+i) = m(j,i);
         end
      end
      
      for i=1:6
         for j=1:7
            jac(14+j,21+i) = gp(i,j);
         end
      end
      
      for i=1:7
         for j=1:6
            jac(21+j,i) = gp(j,i);
         end
      end
      
     
%disp('left pderv')
end

