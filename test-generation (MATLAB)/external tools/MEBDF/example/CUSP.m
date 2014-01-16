%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE CUSP PROBLEM. WITH F, MAS AND PDERV TO BE 
% SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function CUSP
global MYtagSTIFF
MYtagSTIFF = 1;
%MAXN = 400;
ND=96;
LWORK =(33+3*ND)*ND+3;
LIWORK = ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
ML=3;
        MU=3;
        MBND(1)=ML;
        MBND(2)=MU;
        MBND(3)=ML+MU+1;
        MBND(4)=MBND(3)+ MU;

global NERVES NNERV DIFFCOEF DIFFUS


%MBND = zeros(1,4);
TRUE  = zeros(1,ND);
ERROR = zeros(1,ND);
MASBND = zeros(1,4);
%MYY = zeros(1,20);
MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;

%...  DIMENSION OF THE SYSTEM
NNERV = 32;
DIFFUS = NNERV*NNERV/144;
N=3*NNERV
IDID =1;
IERR = 0;



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




      
LOUT = fopen('CUSP.txt','w');

fprintf(LOUT,'\r COMPUTATIONAL STATISTICS OF THE CUSP PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r NUMBER OF EQUATIONS : %g \r',ND);

%... ENDPOINT OF INTEGRATION
    
         XEND=1.1;
         XOUT=XEND;
         X=0;
%...  REQUIRED TOLERANCE

        RTOL=1E-7;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
ANERV = NNERV;
DEL = 2*pi/ANERV;

for INERV=1:NNERV
           Y(3*INERV-2)=0;
           Y(3*INERV-1)=-2*cos(INERV*DEL);
        Y(3*INERV)=2*sin(INERV*DEL);
end


%...  SET DEFAULT VALUES
        MF=24;
        INDEX=1;
        MAXDER=7;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=200000;
        H=1E-8;
        fprintf(LOUT,'Hstart = %E',H);
        MASBND(1)=0;
        

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
           



           TRUE(1)=-1.335038235173363825;
      TRUE(2)=-0.141920661299976116;
      TRUE(3)=2.189999851122752954;
      TRUE(4)=-1.290165517136865728;
      TRUE(5)=0.292210513241939329;
      TRUE(6)=2.524498007953815718;
      TRUE(7)=-1.206268463248866837;
      TRUE(8)=0.702876002804259450;
      TRUE(9)=2.603037671957832137;
      TRUE(10)=-1.081173370722796200;
      TRUE(11)=1.054547339698463687;
      TRUE(12)=2.403900155664309457;
      TRUE(13)=-0.922551477213655165;
      TRUE(14)=1.326991956338080918;
      TRUE(15)=2.009305096775369514;
      TRUE(16)=-0.743049818521982534;
      TRUE(17)=1.516881284521938182;
      TRUE(18)=1.537256339189765457;
      TRUE(19)=-0.555201077072863929;
      TRUE(20)=1.632603197056899916;
      TRUE(21)=1.077437487481636876;
      TRUE(22)=-0.369158363066035655;
      TRUE(23)=1.687674223256960258;
      TRUE(24)=0.673204001909134905;
      TRUE(25)=-0.192671593795137051;
      TRUE(26)=1.695724385342398648;
      TRUE(27)=0.333758479535656796;
      TRUE(28)=-0.030615931836238017;
      TRUE(29)=1.667262708083148298;
      TRUE(30)=0.050978698243957534;
      TRUE(31)=0.117513584875619235;
      TRUE(32)=1.607508563419502269;
      TRUE(33)=-0.190604748914752998;
      TRUE(34)=0.259898961244445617;
      TRUE(35)=1.514823442340568479;
      TRUE(36)=-0.411323787318084589;
      TRUE(37)=0.411809029672400558;
      TRUE(38)=1.379804789392094260;
      TRUE(39)=-0.638121744946811883;
      TRUE(40)=0.590441346230457812;
      TRUE(41)=1.185589061664514451;
      TRUE(42)=-0.905945997151089418;
      TRUE(43)=0.803741778414404449;
      TRUE(44)=0.910756427168161885;
      TRUE(45)=-1.251345457775128863;
      TRUE(46)=1.037877442048335801;
      TRUE(47)=0.545036626743780133;
      TRUE(48)=-1.683821753687386274;
      TRUE(49)=1.239043542405442416;
      TRUE(50)=0.169981336507011738;
      TRUE(51)=-2.112958754094543822;
      TRUE(52)=1.406385681620871097;
      TRUE(53)=-0.235380986562835855;
      TRUE(54)=-2.450796096861493262;
      TRUE(55)=1.524334200774267799;
      TRUE(56)=-0.633461856049010260;
      TRUE(57)=-2.576413161518578492;
      TRUE(58)=1.588649099727842025;
      TRUE(59)=-0.986582203794960052;
      TRUE(60)=-2.442161394270367795;
      TRUE(61)=1.606022353430074771;
      TRUE(62)=-1.269240297385074114;
      TRUE(63)=-2.104018859237057304;
      TRUE(64)=1.588788794126354791;
      TRUE(65)=-1.473056296837721718;
      TRUE(66)=-1.670122571729852451;
      TRUE(67)=1.549115780473624505;
      TRUE(68)=-1.603417743271582846;
      TRUE(69)=-1.233609811984700428;
      TRUE(70)=1.495889929838369103;
      TRUE(71)=-1.672805947347039028;
      TRUE(72)=-0.844976238622025912;
      TRUE(73)=1.434154221021214460;
      TRUE(74)=-1.695067644864453216;
      TRUE(75)=-0.518751841694241591;
      TRUE(76)=1.365334914988092461;
      TRUE(77)=-1.681659890115195530;
      TRUE(78)=-0.249120054639391870;
      TRUE(79)=1.286403800980685275;
      TRUE(80)=-1.639285126097438457;
      TRUE(81)=-0.019980596158691364;
      TRUE(82)=1.184974025791693888;
      TRUE(83)=-1.567910985925968939;
      TRUE(84)=0.194039539947058944;
      TRUE(85)=1.011140518164431143;
      TRUE(86)=-1.455860565434940201;
      TRUE(87)=0.436843623543739620;
      TRUE(88)=-1.349821324547813729;
      TRUE(89)=-1.223845158570813537;
      TRUE(90)=0.809099908070360854;
      TRUE(91)=-1.355008974443540401;
      TRUE(92)=-0.926131110369012061;
      TRUE(93)=1.232945832067604962;
      TRUE(94)=-1.352261107347051711;
      TRUE(95)=-0.559070645046367311;
      TRUE(96)=1.716745798614099647;
      
for t = 1:96
         fprintf(LOUT,'\r\r APROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE Y%g = %E\r\r',t,TRUE(t));
end      

      fclose(LOUT)
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,T,Y,DF,IPAR,RPAR,IERR] = F(N,T,Y,DF,IPAR,RPAR,IERR)
  %F subroutine for the Cusp problem    
     global NERVES NNERV DIFFCOEF DIFFUS
  %disp('now in F')
  
      for INERV=1:NNERV
         X=Y(3*INERV-2);
         A=Y(3*INERV-1);
         B=Y(3*INERV);
         
         if (INERV == 1)
            XRIGHT=Y(3*NNERV-2);
            ARIGHT=Y(3*NNERV-1);
            BRIGHT=Y(3*NNERV);
         else
            XRIGHT=Y(3*INERV-5);
            ARIGHT=Y(3*INERV-4);
            BRIGHT=Y(3*INERV-3);
         end
         
         if (INERV == NNERV)
            XLEFT=Y(1);
            ALEFT=Y(2);
            BLEFT=Y(3);
         else
             
            XLEFT=Y(3*INERV+1);
            ALEFT=Y(3*INERV+2);
            BLEFT=Y(3*INERV+3);
            
         end
         
         XDOT=-10000*(B+X*(A+X*X));
         U=(X-0.7)*(X-1.3);
         V=U/(U+0.1);
         ADOT=B+0.07*V;
         BDOT=(1*(1-A*A)*B-A)-0.4*X+0.035*V;
         DF(3*INERV-2)=XDOT+DIFFUS*(XLEFT-2*X+XRIGHT);
         DF(3*INERV-1)=ADOT+DIFFUS*(ALEFT-2*A+ARIGHT);
         DF(3*INERV)  =BDOT+DIFFUS*(BLEFT-2*B+BRIGHT);
      end

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

function [X,Y,DFY,N,MEBAND,IPAR,RPAR,ierr] =  PDERV(X,Y,DFY,N,MEBAND,IPAR,RPAR,ierr)
%disp('now in pderv')

% DUMMY ROUTINE  

%disp('left pderv')
end

