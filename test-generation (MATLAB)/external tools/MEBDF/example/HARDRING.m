%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE RING MODULATOR PROBLEM WITH cs=2e-12.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function HARDRING
global MYtagSTIFF  PI %MYcounter
%PI = 4*atan(1);
PI = 3.14159265358979324;
MYcounter = 0;
MYtagSTIFF = 1;
ND=15;
LWORK=(31+3*ND)*ND+3;
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
      
LOUT = fopen('hardring.txt','w');
fprintf(LOUT,'RESULTS FOR THE HARD RING PROBLEM');

%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=1E-3;
         XOUT=XEND;
         
%...  REQUIRED TOLERANCE

        RTOL=1E-7;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        Y = zeros(1,N);
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=5000000;
        H=1E-5;
        fprintf(LOUT,'\rHstart = %E',H);
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
           
           if (INDEX == 1)
               INDEX = 0;
               MYcase = 220;
               continue
           end
           
          MYtag = 34;
          continue
                
        end
    end
    
      TRUE(1)=-0.233905735842070174E-01 ;
      TRUE(2)= -0.736748548598215210E-02 ;
      TRUE(3)=0.258295670891733276;     
      TRUE(4)=-0.406446572164108511;     
      TRUE(5)=-0.403945566551075719 ;    
      TRUE(6)=0.260796676505321678;     
      TRUE(7)=0.110676186126208484;     
      TRUE(8)=0.293990434238965346E-06; 
      TRUE(9)=-0.284002993248608379E-07; 
      TRUE(10)=0.726719826722729682E-03 ;
      TRUE(11)=0.792948719688146885E-03 ;
      TRUE(12)=-0.725528349561929123E-03; 
      TRUE(13)=-0.794140196849026266E-03; 
      TRUE(14)=0.708849541688199994E-04; 
      TRUE(15)=0.239005907525775679E-04; 



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

function [n,tn,y,dy,ipar,rpar,ier] = F(n,tn,y,dy,ipar,rpar,ier)
%disp('now in F')
global PI
c=1.6E-8;
cs=2E-12;
cp=1E-8;
r=25E3;
rp=50;
lh=4.45;
ls1=2E-3;
ls2=5E-4;
ls3=5E-4;
rg1=36.3;
rg2=17.3;
rg3=17.3;
ri=5E1;
rc=6E2;
gamma=40.67286402E-9;
delta=17.7493332;


      uin1=0.5E0*sin(2.0E3*PI*tn);
      uin2=2.0E0*sin(2.0E4*PI*tn);
      ud1=+y(3)-y(5)-y(7)-uin2;
      ud2=-y(4)+y(6)-y(7)-uin2;
      ud3=+y(4)+y(5)+y(7)+uin2;
      ud4=-y(3)-y(6)+y(7)+uin2;
      qud1=gamma*(exp(delta*ud1)-1.0E0);
      qud2=gamma*(exp(delta*ud2)-1.0E0);
      qud3=gamma*(exp(delta*ud3)-1.0E0);
      qud4=gamma*(exp(delta*ud4)-1.0E0);
     
      dy(1)=(y(8)-0.5E0*y(10)+0.5E0*y(11)+y(14)-y(1)/r)/c;
      dy(2)=(y(9)-0.5E0*y(12)+0.5E0*y(13)+y(15)-y(2)/r)/c;
      dy(3)=(y(10)-qud1+qud4)/cs;
      dy(4)=(-y(11)+qud2-qud3)/cs;
      dy(5)=(y(12)+qud1-qud3)/cs;
      dy(6)=(-y(13)-qud2+qud4)/cs;
      dy(7)=(-y(7)/rp+qud1+qud2-qud3-qud4)/cp;
      dy(8)=-y(1)/lh;
      dy(9)=-y(2)/lh;
      dy(10)=(0.5E0*y(1)-y(3)-rg2*y(10))/ls2;
      dy(11)=(-0.5E0*y(1)+y(4)-rg3*y(11))/ls3;
      dy(12)=(0.5E0*y(2)-y(5)-rg2*y(12))/ls2;
      dy(13)=(-0.5E0*y(2)+y(6)-rg3*y(13))/ls3;
      dy(14)=(-y(1)+uin1-(ri+rg1)*y(14))/ls1;
      dy(15)=(-y(2)-(rc+rg1)*y(15))/ls1;
      
      
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

function [tn,y,pd,n,meband,ipar,rpar,ier] = PDERV(tn,y,pd,n,meband,ipar,rpar,ier)
%disp('now in pderv')


c=1.6E-8;cs=2E-12;cp=1E-8;r=25E3;rp=50;
lh=4.45;ls1=2E-3;ls2=5E-4;ls3=5E-4;
rg1=36.3;rg2=17.3;rg3=17.3;ri=5E1;rc=6E2;
gamma=40.67286402E-9;delta=17.7493332;

      
%       for i=1:n
%          do 177 j=1,n
%             pd(i,j)=0.0D+0
%       end      

pd = zeros(n,n);



      uin2=2.0E0*sin(2.0E4*pi*tn);  
      ud1=+y(3)-y(5)-y(7)-uin2;
      ud2=-y(4)+y(6)-y(7)-uin2;
      ud3=+y(4)+y(5)+y(7)+uin2;
      ud4=-y(3)-y(6)+y(7)+uin2;
      qpud1=gamma*delta*exp(delta*ud1);
      qpud2=gamma*delta*exp(delta*ud2);
      qpud3=gamma*delta*exp(delta*ud3);
      qpud4=gamma*delta*exp(delta*ud4);
      pd(1,1)=-1.0E0/(c*r);
      pd(1,8)=-1.0E0/c;
      pd(1,10)=-0.5E0/c;
      pd(1,11)=-pd(1,10);
      pd(1,14)=pd(1,8);
      pd(2,2)=pd(1,1);
      pd(2,9)=pd(1,8);
      pd(2,12)=pd(1,10);
      pd(2,13)=pd(1,11);
      pd(2,15)=pd(1,14);
      pd(3,3)=(-qpud1-qpud4)/cs;
      pd(3,5)=qpud1/cs;
      pd(3,6)=-qpud4/cs;
      pd(3,7)=(qpud1+qpud4)/cs;
      pd(3,10)=1.0E0/cs;
      pd(4,4)=(-qpud2-qpud3)/cs;
      pd(4,5)=-qpud3/cs;
      pd(4,6)=qpud2/cs;
      pd(4,7)=(-qpud2-qpud3)/cs;
      pd(4,11)=-1.0E0/cs;
      pd(5,3)=qpud1/cs;
      pd(5,4)=-qpud3/cs;
      pd(5,5)=(-qpud1-qpud3)/cs;
      pd(5,7)=(-qpud1-qpud3)/cs;
      pd(5,12)=1.0E0/cs;
      pd(6,3)=-qpud4/cs;
      pd(6,4)=qpud2/cs;
      pd(6,6)=(-qpud2-qpud4)/cs;
      pd(6,7)=(qpud2+qpud4)/cs;
      pd(6,13)=-1.0E0/cs;
      pd(7,3)=(qpud1+qpud4)/cp;
      pd(7,4)=(-qpud2-qpud3)/cp;
      pd(7,5)=(-qpud1-qpud3)/cp;
      pd(7,6)=(qpud2+qpud4)/cp;
      pd(7,7)=(-qpud1-qpud2-qpud3-qpud4-1.0E0/rp)/cp;
      pd(8,1)=-1.0E0/lh;
      pd(9,2)=pd(8,1);
      pd(10,1)=0.5E0/ls2;
      pd(10,3)=-1.0E0/ls2;
      pd(10,10)=-rg2/ls2;
      pd(11,1)=-0.5E0/ls3;
      pd(11,4)=1.0E0/ls3;
      pd(11,11)=-rg3/ls3;
      pd(12,2)=pd(10,1);
      pd(12,5)=pd(10,3);
      pd(12,12)=pd(10,10);
      pd(13,2)=pd(11,1);
      pd(13,6)=pd(11,4);
      pd(13,13)=pd(11,11);
      pd(14,1)=-1.0E0/ls1;
      pd(14,14)=-(ri+rg1)/ls1;
      pd(15,2)=pd(14,1);
      pd(15,15)=-(rc+rg1)/ls1;
%disp('left pderv')
end
