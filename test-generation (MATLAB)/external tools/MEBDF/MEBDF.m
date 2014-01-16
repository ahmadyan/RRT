function [N, T0, HO, Y0, TOUT, TEND, MF, IDID, LWORK, WORK_1, ...
    WORK_2, WORK_3, WORK_I1, WORK_I2, WORK_I3, WORK_I4, WORK_I5, ...
    WORK_I6, WORK_I7, WORK_I8, WORK_I9, WORK_I10, WORK_I11, LIWORK, ...
    IWORK, IWORK_15, MBND, MASBND, MAXDER, ITOL, RTOL, ATOL, RPAR, ...
    IPAR, IERR] = MEBDF(N, T0, HO, Y0, TOUT, TEND, MF, IDID, ...
    LWORK, WORK_1, WORK_2, WORK_3, WORK_I1, ...
    WORK_I2, WORK_I3, WORK_I4, WORK_I5, WORK_I6, ...
    WORK_I7, WORK_I8, WORK_I9, WORK_I10, WORK_I11, ...
    LIWORK, IWORK, IWORK_15, MBND, MASBND, ...
    MAXDER, ITOL, RTOL, ATOL, RPAR, IPAR, ...
    IERR)

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$%
%                 A MATLAB version of MEBDF.f.                            %
%                                                                         %
%                   Translated by O. Selim,                               %
%                   Student at,                                           %
%                   Department of Mathematics,                            %
%                   Imperial College,                                     %
%                   London SW7 2AZ                                        %
%                   England                                               %
%                                                                         %
%       Please direct any questions and comments(good or bad)to:          %
%                            j.cash@ic.ac.uk                              %
%                                                                         %
%       The original comments from MEBDF.f have been left unchanged       %
%       and new comments have been added where necessary.                 %
%                                                                         %
%       Thank you to Prof. J.R. Cash and also to Dr J. Rogal-Salazar      %
%       for their unreserved help with this project. 08/10/2005.          %
%                                                                         %
%                                                                         %
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$%

% c
% c ***********************************************************************
% c ***********************************************************************
% c    Written by J.R. Cash,
% c    Department of Mathematics,
% c    Imperial College,
% c    London SW7 2AZ
% c    England
% c
% c          j.cash@ic.ac.uk
% c
% c    The author would be pleased to receive any comments
% c          good or bad!!
% c  This has been modified to take in comments from Walter Lioen and
% C  Jacques de Swart from Amsterdam.  Their many comments and suggestions
% C  have removed several bugs and transformed the code.  The author
% C  wishes to record his extreme gratitude to them.
% c ***********************************************************************
% c ***********************************************************************
% c
% C     THIS IS THE NOVEMBER 6th 1998 VERSION OF OVDRIV, A PACKAGE FOR
% C     THE SOLUTION OF THE INITIAL VALUE PROBLEM FOR SYSTEMS OF
% C     ORDINARY DIFFERENTIAL EQUATIONS
% C     DY/DT = F(Y,T),    Y=(Y(1),Y(2),Y(3), . . . ,Y(N))
% C     AND LINEARLY IMPLICIT DIFFERENTIAL ALGEBRAIC EQUATIONS
% C        M(DY/DT) = F(Y,T)
% C     SUBROUTINE OVDRIV IS A DRIVER ROUTINE FOR THIS PACKAGE
% C
% C                    REFERENCES
% C
% C     1.  J. R. CASH, THE INTEGRATION OF STIFF INITIAL VALUE PROBLEMS
% C         IN O.D.E.S USING MODIFIED EXTENDED BACKWARD DIFFERENTIATION
% C         FORMULAE, COMP. AND MATHS. WITH APPLICS., 9, 645-657, (1983).
% C     2.  J.R. CASH AND S. CONSIDINE, AN MEBDF CODE FOR STIFF
% C         INITIAL VALUE PROBLEMS, ACM TRANS MATH SOFTWARE, 142-158,
% C         (1992).
% C     3.  J.R. CASH, STABLE RECURSIONS WITH APPLICATIONS TO THE
% C         NUMERICAL SOLUTION OF STIFF SYSTEMS, ACADEMIC PRESS,(1979).
% C     4.  A.C. HINDMARSH, ODEPACK, A SYSTEMISED COLLECTION OF ODE
% C         SOLVERS, in SCIENTIFIC COMPUTING, R.S. STEPLEMAN et. al.
% C         (eds) North-Holland, AMSTERDAM, pp55-64 , (1983).
% C     5.  E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
% C         EQUATIONS II, STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS,
% C         SPRINGER 1996, page 267.
% C
% C     ----------------------------------------------------------------
% C     OVDRIV IS TO BE CALLED ONCE FOR EACH OUTPUT VALUE OF T, AND
% C     IN TURN MAKES REPEATED CALLS TO THE CORE INTEGRATOR STIFF.
% C
% C     THE INPUT PARAMETERS ARE ..
% C     N     =  THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS.
% C     T0    =  THE INITIAL VALUE OF T, THE INDEPENDENT VARIABLE
% C              (USED ONLY ON THE FIRST CALL)
% C     HO    =  THE NEXT STEP SIZE IN T (USED FOR INPUT ONLY ON THE
% C              FIRST CALL)
% C     Y0    =  A VECTOR OF LENGTH N CONTAINING THE INITIAL VALUES OF Y
% C              (USED FOR INPUT ONLY ON FIRST CALL)
% C     TOUT  =  THE VALUE OF T AT WHICH OUTPUT IS DESIRED NEXT.
% C              INTEGRATION WILL NORMALLY GO SLIGHTLY BEYOND TOUT
% C              AND THE PACKAGE WILL INTERPOLATE TO T = TOUT
% C     TEND  =  END OF THE RANGE OF INTEGRATION.
% C     MF    =  THE METHOD FLAG.  AT PRESENT MF=21,22,23 OR 24 IS
% C              ALLOWED. THESE ARE EXTENDED BACKWARD DIFFERENTIATION
% C              FORMULAE USING THE CHORD METHOD WITH ANALYTIC OR NUMERICAL
% C              JACOBIAN FOR MF=21,22 RESPECTIVELY. MF=23/24 ARE  THE SAME
% C              AS FOR 21/22 BUT THE JACOBIAN IS NOW BANDED.   THE USER
% C              NEEDS TO SPECIFY SUBROUTINE PDERVIF MF=21 OR 23.
% C     IDID   = THE INTEGER USED ON INPUT TO INDICATE THE TYPE OF CALL.
% C              1   THIS IS THE FIRST CALL FOR THE PROBLEM.
% C              0   THIS IS NOT THE FIRST CALL FOR THIS PROBLEM
% C                  AND INTEGRATION IS TO CONTINUE.
% C             -1   THIS IS NOT THE FIRST CALL FOR THE PROBLEM,
% C                  AND THE USER HAS RESET N, RTOL, ATOL,H  AND/OR MF.
% C              2   SAME AS 0 EXCEPT THAT TOUT HAS TO BE HIT
% C                  EXACTLY (NO INTERPOLATION IS DONE).
% C                  ASSUMES TOUT .GE. THE CURRENT T.
% C              3   SAME AS 0 EXCEPT CONTROL RETURNS TO CALLING
% C                  PROGRAM AFTER ONE STEP. TOUT IS IGNORED, UNTIL THE
% C                  INTEGRATION REACHES TOUT OR BEYOND. IF IT PASSES TOUT
% C                  THE PROGRAM INTERPOLATES THE SOLUTION VALUES AND
% C                  RETURNS THE SOLUTION VALUE AT TOUT.
% C              SINCE THE NORMAL OUTPUT VALUE OF IDID IS 0,
% C              IT NEED NOT BE RESET FOR NORMAL CONTINUATION.
% C              THE FIRST CALL TO THE DRIVER IS WITH IDID=1 AND FOR
% C              A SUCCESSFUL STEP THE DRIVER RETURNS WITH IDID=1.THUS
% C              THE CALL WITH IDID = 1 IS SIMPLY THE FIRST
% C              INITIALISING STEP FOR THE CODE.  THE USER
% C              THEN NEEDS TO CONTINUE WITH IDID=0,-1,2 OR 3 AS ABOVE.
% C     LOUT   = THE LOGICAL OUTPUT CHANNEL FOR MESSAGE PASSING.
% C     MBND   = AN ARRAY OF DIMENSION 4 FOR USE WHEN THE JACOBIAN
% C              MATRIX IS BANDED.  IF THE JACOBIAN HAS ML DIAGONALS
% C              BELOW THE MAIN DIAGONAL AND MU DIAGONALS ABOVE THE
% C              MAIN DIAGONAL THEN:
% C              MBND(1) = ML
% C              MBND(2) = MU
% C              MBND(3) = MU + ML + 1
% C              MBND(4) = 2*ML + MU + 1
% C     MASBND = AN ARRAY OF DIMENSION 4 DESCRIBING THE MASS MATRIX.
% C              MASBND(1) = IMAS IS ZERO IF THE MASS MATRIX IS THE
% C              IDENTITY. OTHERWISE IMAS = 1.
% C              MASBND(2) = N IF THE MASS MATRIX IS FULL. OTHERWISE
% C             = MLMAS THE LOWER BANDWIDTH.
% C              MASBND(3) = MUMAS THE UPPER BANDWIDTH.  IF THE MASS MATRIX
% C              IS FULL THIS NEED NOT BE SPECIFIED.
% C              MASBND(4) = LMAS = N IF THE MASS MATRIX IS FULL.
% C                                = 0 IF IMAS = 0
% C                                = MLMAS + MUMAS + 1 OTHERWISE
% C              THE RESTRICTION MLMAS .LE. ML AND MUMAS .LE. MU IS
% C              ALSO IMPOSED.
% C     MAXDER=  THE MAXIMUM ORDER IS MAXDER + 1.
% C              THE VALUE OF MAXDER CANNOT EXCEED 7.  THIS IS THE
% C              VALUE RECOMMENDED UNLESS IT IS BELIEVED THAT THERE
% C              ARE SEVERE STABILITY PROBLEMS IN WHICH CASE MAXDER=3
% C              OR 4 SHOULD BE TRIED INSTEAD.
% C     ITOL  =  AN INDICATOR OF THE TYPE OF ERROR CONTROL. SEE
% C              DESCRIPTION BELOW UNDER ATOL.
% C     RTOL  =  A RELATIVE ERROR TOLERANCE PARAMETER. CAN BE EITHER A
% C              SCALAR OR AN ARRAY OF LENGTH N.  SEE DESCRIPTION
% C              BELOW UNDER ATOL.
% C     ATOL  =  THE ABSOLUTE ERROR BOUND.
% C              THE INPUT PARAMETERS ITOL, RTOL AND ATOL DETERMINE
% C              THE ERROR CONTROL PERFORMED BY THE SOLVER.  THE
% C              SOLVER WILL CONTROL THE VECTOR e = (e(i)) OF ESTIMATED
% C              LOCAL ERRORS IN y ACCORDING TO AN INEQUALITY OF THE FORM
% C                  RMS-NORM OF (e(i)/ewt(i)) .LE. 1
% C              THE ROOT MEAN SQUARE NORM IS
% C                   RMS-NORM(V) = SQRT((SUM v(i)**2)/N).  HERE
% C                ewt = (ewt(i)) IS A VECTOR OF WEIGHTS WHICH MUST
% C              ALWAYS BE POSITIVE, AND THE VALUES OF RTOL AND ATOL
% C              SHOULD BE NON-NEGATIVE. IF ITOL = 1 THEN SINGLE STEP ERROR
% C              ESTIMATES DIVIDED BY YMAX(I) WILL BE KEPT LESS THAN 1
% C              IN ROOT-MEAN-SQUARE NORM.  THE VECTOR YMAX OF WEIGHTS IS
% C              COMPUTED IN OVDRIV. INITIALLY YMAX(I) IS SET AS
% C              THE MAXIMUM OF 1 AND ABS(Y(I)).  THEREAFTER YMAX(I) IS
% C              THE LARGEST VALUE OF ABS(Y(I)) SEEN SO FAR, OR THE
% C              INITIAL VALUE YMAX(I) IF THAT IS LARGER.
% C              IF ITOL = 1 THE USER NEEDS TO SET ATOL = RTOL =
% C              THE PRECISION REQUIRED.  THEN
% C                         ewt(i) = RTOL(1)*YMAX(i)
% C              IF ITOL IS GREATER THAN 1 THEN
% C                 ewt(i) = rtol(i)*abs(y(i)) + atol(i)
% C              THE FOLLOWING TABLE GIVES THE TYPES (SCALAR/ARRAY)
% C              OF RTOL AND ATOL, AND THE CORRESPONDING FORM OF ewt(i)
% C                  ITOL   RTOL      ATOL       ewt(i)
% C                   2    SCALAR    SCALAR   rtol*abs(y(i))   + atol
% C                   3    SCALAR    ARRAY    rtol*abs(y(i))   + atol(i)
% C                   4    ARRAY     SCALAR   rtol(i)*abs(y(i))+ atol
% C                   5    ARRAY     ARRAY    rtol(i)*abs(y(i))+ atol(i)
% C              IF EITHER OF THESE PARAMETERS IS A SCALAR, IT NEED
% C              NOT BE DIMENSIONED IN THE USER'S CALLING PROGRAM.
% C     NIND1    = THE NUMBER OF VARIABLES OF INDEX 1,2,3 RESPECTIVELY.
% C     NIND2,NIND3  THESE ARE SET IN IWORK(1),(2),(3).
% C              THE EQUATIONS MUST BE DEFINED SO THAT THE INDEX 1
% C              VARIABLES PRECEDE THE INDEX 2 VARIABLES WHICH IN
% C              TURN PRECEDE THE INDEX 3 VARIABLES.
% C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS)
% C              WHICH CAN BE USED FOR COMMUNICATION BETWEEN  THE USER'S
% C              CALLING PROGRAM AND THE F, PDERV AND MAS SUBROUTINES.
% C     IERR     IERR IS AN INTEGER FLAG WHICH IS ALWAYS EQUAL TO ZERO
% C              ON INPUT.  SUBROUTINES F, PDERV AND MAS SHOULD ALTER
% C              IERR ONLY IF ONE OF THEM ENCOUNTERS AN ILLEGAL OPERATION SUCH
% C              AS THE SQUARE ROOT OF A NEGATIVE NUMBER OR EXPONENT
% C              OVERFLOW. THE USER CAN THEN ALTER H AND CALL THE
% C              SUBROUTINE AGAIN WITH IDID=-1 IF HE WISHES.
% C
% C     AFTER THE INITIAL CALL, IF A NORMAL RETURN OCCURED AND A NORMAL
% C     CONTINUATION IS DESIRED, SIMPLY RESET TOUT AND CALL AGAIN.
% C     ALL OTHER PARAMETERS WILL BE READY FOR THE NEXT CALL.
% C     A CHANGE OF PARAMETERS WITH IDID = -1 CAN BE MADE AFTER
% C     EITHER A SUCCESSFUL OR AN UNSUCCESSFUL RETURN.
% C
% C     THE OUTPUT PARAMETERS ARE..
% C     T0    =  THE VALUE OF T WHICH RELATES TO THE CURRENT SOLUTION
% C              POINT Y0()
% C     HO    =  THE STEPSIZE H USED LAST, WHETHER SUCCESSFULLY OR NOT.
% C     Y0    =  THE COMPUTED VALUES OF Y AT T = TOUT
% C     TOUT  =  UNCHANGED FROM ITS INPUT VALUE.
% C     IDID  =  INTEGER USED ON OUTPUT TO INDICATE RESULTS, WITH
% C              THE FOLLOWING VALUES AND MEANINGS..
% C
% C      0   INTEGRATION WAS COMPLETED TO TOUT OR BEYOND.
% C
% C     -1   THE INTEGRATION WAS HALTED AFTER FAILING TO PASS THE
% C          ERROR TEST EVEN AFTER REDUCING H BY A FACTOR OF
% C          1.E10 FROM ITS INITIAL VALUE.
% C
% C     -2   AFTER SOME INITIAL SUCCESS, THE INTEGRATION WAS
% C          HALTED EITHER BY REPEATED ERROR TEST FAILURES OR BY
% C          A TEST ON RTOL/ATOL.  TOO MUCH ACCURACY HAS BEEN REQUESTED.
% C
% C     -3   THE INTEGRATION WAS HALTED AFTER FAILING TO ACHIEVE
% C          CORRECTOR CONVERGENCE EVEN AFTER REDUCING H BY A
% C          FACTOR OF 1.E10 FROM ITS INITIAL VALUE.
% C
% C     -4   IMMEDIATE HALT BECAUSE OF ILLEGAL VALUES OF INPUT
% C          PARAMETERS.  SEE PRINTED MESSAGE.
% C
% C     -5   IDID WAS -1 ON INPUT, BUT THE DESIRED CHANGES OF
% C          PARAMETERS WERE NOT IMPLEMENTED BECAUSE TOUT
% C          WAS NOT BEYOND T.  INTERPOLATION AT T = TOUT WAS
% C          PERFORMED AS ON A NORMAL RETURN.  TO TRY AGAIN,
% C          SIMPLY CALL AGAIN WITH IDID = -1 AND A NEW TOUT.
% C
% C     -6   MAXIMUM ALLOWABLE NUMBER OF INTEGRATION STEPS EXCEEDED.
% C          TO CONTINUE THE USER SHOULD RESET IWORK(14).
% C
% C
% C     -7   STEPSIE IS TOO SMALL (LESS THAN SQRT(UROUND)/100)
% C
% C
% C     -11   INSUFFICIENT REAL WORKSPACE FOR THE INTEGRATION
% C
% C     -12   INSUFFICIENT INTEGER WORKSPACE FOR THE INTEGRATION
% C
% C
% C     IN ADDITION TO OVDRIVE, THE FOLLOWING ROUTINES ARE PROVIDED
% C     IN THE PACKAGE..
% C
% C     INTERP( - )   INTERPOLATES TO GET THE OUTPUT VALUES
% C                   AT T=TOUT FROM THE DATA IN THE Y ARRAY.
% C     STIFF( - )    IS THE CORE INTEGRATOR ROUTINE.  IT PERFORMS A
% C                   SINGLE STEP AND ASSOCIATED ERROR CONTROL.
% C     COSET( - )    SETS COEFFICIENTS FOR BACKWARD DIFFERENTIATION
% C                   SCHEMES FOR USE IN THE CORE INTEGRATOR.
% C     PSET( - )     COMPUTES AND PROCESSES THE JACOBIAN
% C                   MATRIX J = DF/DY
% C     DEC( - )      PERFORMS AN LU DECOMPOSITION ON A MATRIX.
% C     SOL( - )      SOLVES LINEAR SYSTEMS A*X = B AFTER DEC
% C                   HAS BEEN CALLED FOR THE MATRIX A
% C     DGBFA ( - )   FACTORS A DOUBLE PRECISION BAND MATRIX BY
% C                   ELIMINATION.
% C     DGBSL ( - )   SOLVES A BANDED LINEAR SYSTEM A*x=b
% C
% C                   ALSO SUPPLIED ARE THE BLAS ROUTINES
% C
% C                   daxpy, dscal, idamax, ddot.
% C
% C
% C     THE FOLLOWING ROUTINES ARE TO BE SUPPLIED BY THE USER AND
% C                   SHOULD BE DECLARED AS EXTERNAL.
% C
% C     F(N,T,Y,YDOT,IPAR,RPAR)   COMPUTES THE FUNCTION YDOT = F(Y,T),
% C                         THE RIGHT HAND SIDE OF THE O.D.E.
% C                         HERE Y AND YDOT ARE VECTORS OF LENGTH N
% C     PDERV(T,Y,PD,N,MEBAND,IPAR,RPAR)  COMPUTES THE N*N JACOBIAN MATRIX
% C                         OF PARTIAL DERIVATIVES AND STORES IT IN PD
% C                         AS AN N BY N ARRAY IF THE JACOBIAN IS FULL.
% C                         IF THE JACOBIAN IS BANDED THE ARRAY PD IS OF
% C                         SIZE MBND(4)*N.
% C                         IF THE JACOBIAN IS FULL PD(I,J) IS TO BE SET
% C                         TO THE PARTIAL DERIVATIVE OF YDOT(I) WITH
% C                         RESPECT TO Y(J).  IF THE JACOBIAN IS BANDED
% C                         WITH mu DIAGONALS ABOVE THE MAIN DIAGONAL
% C                         THE PARTIAL DERIVATIVE DF(I)/DY(J) SHOULD BE
% C                         PUT IN PD(i-j+mu+1,j). PDERV IS CALLED ONLY IF
% C                         MITER = 1 OR 3.  OTHERWISE A DUMMY ROUTINE CAN
% C                         BE SUBSTITUTED.
% C
% C     MAS(N,AM,MASBND(4),IPAR,RPAR) COMPUTES THE MASS MATRIX M.  IF IMAS=0
% C                      THIS MATRIX IS THE IDENTITY MATRIX AND THE
% C                      SUBROUTINE MAS CAN BE DUMMY.  IF IMAS=1
% C                      THE SUBROUTINE MAS IS OF THE FORM
% C                      MAS(N,AM,MASBND(4))
% C                      IF MASBND(4) = N THE MASS MATRIX IS STORED AS
% C                      A FULL MATRIX WITH AM(I,J) = M(I,J). OTHERWISE
% C                      THE MATRIX IS BANDED AND IS STORED
% C                      COMPONENTWISE AS AM(I-J+MUMAS+1,J)=M(I,J).
% C
% C     THE DIMENSION OF PD/ MAS  MUST BE AT LEAST N**2 FOR A FULL JACOBIAN
% C     /MASS MATRIX AND MBND(4)*N/ MASBND(4)*N FOR A BANDED JACOBIAN
% C     / MASS MATRIX(MF=23 OR 24) .  THE DIMENSIONS
% C     OF YMAX,ERROR,SAVE1,SAVE2,IPIV AND THE FIRST DIMENSION
% C     OF Y SHOULD ALL BE AT LEAST N.
% C
% C
% C     UROUND   THIS IS THE UNIT ROUNDOFF AND HAS TO BE SET AS
% C                        UROUND = DLAMCH('Epsilon')
% C     EPSJAC   = sqrt(UROUND).
% C
% C     HUSED  (=WORK(2))    LAST STEPSIZE SUCCESSFULLY USED BY THE INTEGRATOR
% C     NQUSED (=IWORK(4))   LAST ORDER SUCCESSFULLY USED
% C     NSTEP  (=IWORK(5))   NUMBER OF SUCCESSFUL STEPS TAKEN SO FAR
% C     NFAIL  (=IWORK(6))   NUMBER OF FAILED STEPS
% C     NFE    (=IWORK(7))   NUMBER OF FUNCTION EVALUATIONS  SO FAR
% C     NJE    (=IWORK(8))   NUMBER OF JACOBIAN EVALUATIONS  SO FAR
% C     NDEC   (=IWORK(9))   NUMBER OF LU DECOMPOSITIONS  SO FAR
% C     NBSOL  (=IWORK(10))   NUMBER OF 'BACKSOLVES'  SO FAR
% C     NPSET  (=IWORK(11))   NUMBER OF TIMES A NEW COEFFICIENT MATRIX HAS BEEN
% C               FORMED SO FAR
% C     NCOSET (=IWORK(12))   NUMBER OF TIMES THE ORDER OF THE METHOD USED HAS
% C               BEEN CHANGED SO FAR
% C     MAXORD (=IWORK(13))   THE MAXIMUM ORDER USED SO FAR IN THE INTEGRATION
% c     MAXSTP (=IWORK(14))   THE MAXIMUM ALLOWED NUMBER OF STEPS SET BY THE USER
% C
% C    IF IT IS ONLY REQUIRED TO CONTROL THE ACCURACY IN THE
% C    DIFFERENTIAL VARIABLES THEN THE USER SHOULD FIND THE
% C    STRING 'AMMEND' AND MAKE CHANGES THERE
% C
% c

% C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% C     THIS SUBROUTINE IS FOR THE PURPOSE               *
% C     OF SPLITTING UP THE WORK ARRAYS WORK AND IWORK   *
% C     FOR USE INSIDE THE INTEGRATOR STIFF              *
% C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


persistent  EPSJAC

if (IDID == 1)

    if (N <= 0)

        disp(['*** ERROR ***** INTEGRATION HALTED IN DRIVER \n >>>' ...
            ' ILLEGAL VALUE FOR NUMBER OF EQUATIONS <<< \n WITH' ...
            ' N = %d'],N)
        IDID = -4;
    else

        if (MF<23)
            MBND(4) = N;
        end

        UROUND = eps(1/2);

        WORK_1 = UROUND;
        EPSJAC = sqrt(WORK_1);
    end

    if (IDID < 0)
        return
    end

end

% c
% c    THE DIMENSION OF THE REAL WORKSPACE, WORK, HAS TO BE AT LEAST
% c    (32 + 2*MBND(4) + MASBND(4))*N+2 WHILE THE DIMENSION OF THE INTEGER
% c    WORKSPACE HAS TO BE AT LEAST N+14.
% c


IWORK_1  = IWORK(1) ;
IWORK_2  = IWORK(2) ;
IWORK_3  = IWORK(3) ;
IWORK_4  = IWORK(4) ;
IWORK_5  = IWORK(5) ;
IWORK_6  = IWORK(6) ;
IWORK_7  = IWORK(7) ;
IWORK_8  = IWORK(8) ;
IWORK_9  = IWORK(9) ;
IWORK_10 = IWORK(10) ;
IWORK_11 = IWORK(11) ;
IWORK_12 = IWORK(12) ;
IWORK_13 = IWORK(13) ;
IWORK_14 = IWORK(14) ;

%ALL INPUTS PASSED ON TO OVDRIV() ARE INITIALISED AT THE DRIVER STAGE.

[N, T0, HO, Y0, TOUT, TEND, MF, IDID, WORK_3, WORK_I1, WORK_I2, ...
    WORK_I3, WORK_I4, WORK_I5, WORK_I6, WORK_I7, WORK_I8, WORK_I9, WORK_I10, ...
    WORK_I11, IWORK_15, MBND,MASBND, IWORK_1 ,IWORK_2, IWORK_3, MAXDER, ITOL, ...
    RTOL, ATOL, RPAR, IPAR, IWORK_4, IWORK_5, IWORK_6, IWORK_7, IWORK_8, ...
    IWORK_9, IWORK_10, IWORK_11, IWORK_12, IWORK_13, IWORK_14, WORK_1, WORK_2, ...
    EPSJAC, IERR] = OVDRIV(N, T0, HO, Y0, TOUT, TEND, MF, IDID, WORK_3, ...
    WORK_I1, WORK_I2, WORK_I3, WORK_I4, WORK_I5, ...
    WORK_I6, WORK_I7, WORK_I8, WORK_I9, WORK_I10, ...
    WORK_I11, IWORK_15, MBND, MASBND, IWORK_1, IWORK_2, ...
    IWORK_3, MAXDER, ITOL, RTOL, ATOL, RPAR, IPAR, ...
    IWORK_4, IWORK_5, IWORK_6, IWORK_7, IWORK_8, ...
    IWORK_9, IWORK_10, IWORK_11, IWORK_12, IWORK_13, ...
    IWORK_14, WORK_1,WORK_2,EPSJAC,IERR);


IWORK(1) = IWORK_1;
IWORK(2) = IWORK_2;
IWORK(3) = IWORK_3;
IWORK(4) = IWORK_4;
IWORK(5) = IWORK_5;
IWORK(6) = IWORK_6;
IWORK(7) = IWORK_7;
IWORK(8) = IWORK_8;
IWORK(9) = IWORK_9;
IWORK(10) = IWORK_10;
IWORK(11) = IWORK_11;
IWORK(12) = IWORK_12;
IWORK(13) = IWORK_13;
IWORK(14) = IWORK_14;

% C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% C     WORK() HOUSES THE FOLLOWING ARRAYS
% C
% C     Y(N,12)   , YHOLD(N,12) , YNHOLD(N,2) , YMAX(N)
% C     ERRORS(N) , SAVE1(N) , SAVE2(N) , SCALE(N) , ARH(N) , PW(MBND(4)*N)
% C     PWCOPY(MBND(4)*N) ,AM(MASBND(4)*N)
% C     IF THE SPARSE OPTION IS NOT BEING USED THEN MBND(4)=N=MASBND(4).
% C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%========================================================================
%========================================================================



function [N, T0, HO, Y0, TOUT, TEND, MF, IDID, Y, YHOLD, YNHOLD, ...
    YMAX, ERRORS, SAVE1, SAVE2, SCALE, ARH, PW, PWCOPY, AM, IPIV, ...
    MBND, MASBND, NIND1, NIND2, NIND3, MAXDER, ITOL, RTOL, ATOL, ...
    RPAR, IPAR, NQUSED, NSTEP, NFAIL, NFE, NJE, NDEC, NBSOL, NPSET, ...
    NCOSET, MAXORD, MAXSTP, UROUND, HUSED, EPSJAC, IERR] = OVDRIV(N, ...
    T0, HO, Y0, TOUT, TEND, ...
    MF,IDID, Y, YHOLD, ...
    YNHOLD, YMAX, ERRORS, ...
    SAVE1, SAVE2, SCALE, ...
    ARH, PW, PWCOPY, AM, ...
    IPIV, MBND, MASBND, ...
    NIND1, NIND2, NIND3, ...
    MAXDER, ITOL, RTOL, ATOL, ...
    RPAR, IPAR, NQUSED, ...
    NSTEP, NFAIL, NFE, NJE, ...
    NDEC, NBSOL, NPSET, ...
    NCOSET, MAXORD, MAXSTP, ...
    UROUND,HUSED,EPSJAC,IERR)

% C
% C     START OF PROGRAM PROPER
% C

persistent T H HMIN HMAX KFLAG JSTART NHCUT

if (IDID==0)
    %I.E. NORMAL CONTINUATION OF INTEGRATION
    T0=T;
    HMAX = abs(TEND-T0)*10;

    if ((T-TOUT)*H >= 0)
        % HAVE OVERSHOT THE OUTPUT POINT, SO INTERPOLATE
        [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0);
        IDID = KFLAG;
        T0 = TOUT;
        HO = H;
        return
    end

elseif (IDID == 2)
    % C         I.E. CONTINUING INTEGRATION BUT WISH TO HIT TOUT
    T0 = T;
    HMAX = abs(TEND-T0)*10;

    if (((T+H)-TOUT)*H > 0)
        % C             WE HAVE ALREADY OVERSHOT THE OUTPUT POINT OR WE WILL
        % C             DO SO ON THE NEXT STEP

        if (((T-TOUT)*H >= 0)||(abs(T-TOUT) <= 100*UROUND*HMAX))

            % C                 HAVE OVERSHOT THE OUTPUT POINT, SO INTERPOLATE
            [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0);
            T0 = TOUT;
            HO = H;
            IDID = KFLAG;
            return

        else
            % C                 WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE
            % C                 SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY'
            H = (TOUT-T)* (1-4*UROUND);
            JSTART = -1;
        end

    end

elseif (IDID==-1)
    % C         NOT FIRST CALL BUT PARAMETERS RESET
    H = HO;
    if (H < EPSJAC/100)
        disp('\n STEPSIZE IS TOO SMALL \n');         %(6,9160)
        IDID = -7;
        return
    end
    T0 = T;

    if ((T-TOUT)*H >= 0)
        % C             HAVE OVERSHOT TOUT
        disp(['\n IDID = -1 ON INPUT WITH (T-TOUT)*H .GE. 0. \n T =' ...
            ' %f  TOUT = %f  H = %f . \n INTERPOLATION WAS DONE AS' ...
            ' ON NORMAL RETURN. \n DESIRED PARAMETER CHANGES WERE' ...
            ' NOT MADE. \n'], T,TOUT,H);%WRITE(LOUT,9080) T,TOUT,H
        
        [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0);
        HO = H;
        T0 = TOUT;
        IDID = -5;     
        return
    else
        JSTART = -1;
    end

elseif (IDID==3)
    T0 = T;

    if ((T-TOUT)*H >= 0)
        % C             HAVE OVERSHOT,SO INTERPOLATE
        [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0);
        IDID = KFLAG;
        T0 = TOUT;
        HO = H;             
        return
    end

else
    %IDID SHOULD BE 1 AND THIS IS THE FIRST CALL FOR THIS PROBLEM
    %CHECK THE ARGUMENTS THAT WERE PASSED FOR CORRECTNESS
    
    
    if (IDID~=1)
        % C             VALUE OF IDID NOT ALLOWED
        disp('\n ILLEGAL INPUT.. IDID = %d \n',IDID);
        %WRITE (LOUT,9070) IDID
        IDID = -4;
    end

    NN=N;
    
    if (ITOL <= 3)
        NN = 1;
    end

    for I = 1:NN
        if (RTOL(I)<0)
            % C                 ILLEGAL VALUE FOR RELATIVE ERROR TOLERENCE
            disp('\n ILLEGAL INPUT.. RTOL .LE. 0. \n');
             %WRITE (LOUT,9040)
            IDID = -4;
        end
    end

    NN=N;

    if (ITOL==1||ITOL==2||ITOL==4)
        NN=1;
    end

    for I=1:NN
    if (ATOL(1:NN) < 0)
        % C                 ILLEGAL ABSOLUTE ERROR TOLERANCE
        disp('\n ILLEGAL INPUT.. ATOL .LE. 0. \n'); %WRITE(LOUT,9045)
        IDID=-4;

          end
    end

    if ((ITOL==1) && (RTOL(1)== 0))
        % C             ILLEGAL ERROR TOLERANCE
        disp(' ILLEGAL INPUT.. RTOL .LE. 0. \n');
        %WRITE(LOUT,9040)
        IDID = -4;
    end

    if (ITOL ~= 1)
        VHOLD = 0;

        for I = 1:N

            if (ITOL==2)
                VHOLD = max([RTOL(1),ATOL(1)]);
            elseif (ITOL==3)
                VHOLD = max([RTOL(1),ATOL(I)]);
            elseif (ITOL==4)
                VHOLD = max([RTOL(I),ATOL(1)]);
            elseif (ITOL==5)
                VHOLD = max([RTOL(I),ATOL(I)]);
            end

            if(VHOLD <= 0)
                disp('\n ILLEGAL INPUT.. RTOL .LE. 0. \n')
                IDID = -4;
            end
        end
    end

    if (N <= 0)
        % C             ILLEGAL VALUE FOR THE NUMBER OF EQUATIONS
        disp(' ILLEGAL INPUT.. N .LE. 0 \n');
        %WRITE (LOUT,9050)
        IDID = -4;

    end

    if ((T0-TOUT)*HO >= 0)
        % C             PARAMETERS FOR INTEGRATION ARE ILLEGAL
        disp('\n ILLEGAL INPUT.. (T0-TOUT)*H .GE. 0. \n');
        %WRITE (LOUT,9060)
       IDID = -4;
    end

    if ((MF ~= 21)&&(MF ~= 22)&&(MF ~= 23)&&(MF~=24))
        % C             ILLEGAL VALUE FOR METHOD FLAG
        disp(['\n ILLEGAL INPUT.. METHOD FLAG, MF, = %f, ALLOWED' ...
            ' VALUES ARE 21 OR 22 \n'],MF);       %WRITE (LOUT,9090) MF
        IDID = -4;

    end

    if (ITOL < 1 || ITOL > 5)
        %C              ILLEGAL VALUE FOR ERROR CONTROL PARAMETER
        disp('\n ILLEGAL VALUE FOR ITOL \n');      %WRITE (LOUT,9110)
        IDID=-4;

    end

    if (MAXDER < 1 || MAXDER > 7)
        %C              ILLEGAL VALUE FOR MAXIMUM ORDER
        disp('\n ILLEGAL VALUE FOR MAXDER \n');     %WRITE(LOUT,9120)
        IDID = -4;
    end

    if (NIND1 == 0)
        NIND1=N;
    end

    if(NIND1 + NIND2 + NIND3 ~= N)
        %C              SUM OF VARIABLES OF DIFFERENT INDEX SHOULD BE N.
        disp('\n BAD INPUT FOR NUMBER OF VARIABLES OF INDEX 1,2,3 \n');
        %WRITE(LOUT,9140)
        IDID = -4;
    end
   
    if (IDID ~= 1)
 
        return
    else

        % C         THE INITIAL PARAMETERS ARE O.K. SO INITIALISE EVERYTHING
        % C         ELSE NECESSARY FOR THE INTEGRATION.
        % C         IF VALUES OF YMAX OTHER THAN THOSE SET BELOW ARE DESIRED,
        % C         THEY SHOULD BE SET HERE. ALL YMAX(I) MUST BE POSITIVE. IF
        % C         VALUES FOR HMIN OR HMAX, THE BOUNDS ON DABS(H), OTHER THAN
        % C         THOSE BELOW ARE DESIRED, THEY SHOULD BE SET BELOW.

        if(ITOL == 1)
            YMAX(1:N)=max([abs(Y0(1:N)),1]);
        end
      
        Y(1:N,1)=Y0(1:N);
        
        T = T0;
        H = HO;
        HMIN = abs(HO);
        HMAX = abs(T0-TEND)*10;
        JSTART = 0;
        NHCUT = 0;
    end

end

% C     <<<<<<<<<<<<<<<<<
% C     <  TAKE A STEP  >
% C     <<<<<<<<<<<<<<<<<

MYstop = 0;
MYcase = 20;
while (MYstop == 0)

    switch MYcase

        case 20 %**********************************************************
             
            if ((T+H)==T)
                disp('\n WARNING..  T + H = T ON NEXT STEP. \n');
                %WRITE (LOUT,9000)
            end

            [H, HMAX, HMIN, JSTART, KFLAG, MF, MBND, MASBND, NIND1, NIND2, NIND3, ...
                T, TOUT, TEND, Y, N, YMAX, ERRORS, SAVE1, SAVE2, SCALE, PW, PWCOPY, ...
                AM, YHOLD, YNHOLD, ARH, IPIV, MAXDER, ITOL, RTOL, ATOL, RPAR, ...
                IPAR, NQUSED, NSTEP, NFAIL, NFE, NJE, NDEC, NBSOL, NPSET, NCOSET, ...
                MAXORD, MAXSTP, UROUND, EPSJAC, HUSED, IERR] = STIFF(H, HMAX, HMIN, ...
                JSTART, KFLAG, MF, MBND, ...
                MASBND, NIND1, NIND2, ...
                NIND3, T, TOUT, TEND,Y,N, ...
                YMAX,ERRORS,SAVE1,SAVE2, ...
                SCALE,PW,PWCOPY,AM,YHOLD, ...
                YNHOLD,ARH,IPIV, ...
                MAXDER, ITOL,RTOL,ATOL,RPAR, ...
                IPAR,NQUSED,NSTEP,NFAIL,NFE, ...
                NJE,NDEC,NBSOL,NPSET, ...
                NCOSET, MAXORD,MAXSTP, ...
                UROUND,EPSJAC,HUSED,IERR);


            KGO = 1 - KFLAG;

            if (KGO == 1)
                % C        NORMAL RETURN FROM STIFF
                %GOTO 30
                MYcase = 30;
                continue

            elseif (KGO == 2)

                % C        COULD NOT ACHIEVE REQUIRED PRECISION WITH HMIN
                % C        SO CHOP HMIN IF WE HAVEN'T DONE SO 10 TIMES

                %GOTO 60
                MYcase = 60;
                continue

            elseif (KGO == 3)

                % C    ERROR REQUIREMENT SMALLER THAN CAN BE HANDLED FOR THIS PROBLEM
                disp(['\r\r KFLAG = -2 FROM INTEGRATOR AT T = %f  H = %f, \n THE REQUESTED ERROR IS SMALLER THAN CAN BE HANDLED\r \n'],T,H);
                % WRITE (LOUT,9010) T,H

                %GOTO 70
                MYcase = 70;
                continue

            elseif (KGO == 4)

                % C        COULD NOT ACHIEVE CONVERGENCE WITH HMIN
                disp([' \r\rKFLAG = -3 FROM INTEGRATOR AT T = %f \n  CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED \n',T]);
                %WRITE (LOUT,9030)
                %GOTO 60
                MYcase = 60;
                continue
            end
            error('KGO is not 1,2,3 or 4')

        case 30
            %**********************************************************
                 %disp('label 30')
            % C -------------------------------------------------------------------
            % C    NORMAL RETURN FROM THE INTEGRATOR.
            % C
            % C     THE WEIGHTS YMAX(I) ARE UPDATED IF ITOL=1.
            % C    IF DIFFERENT VALUES ARE DESIRED, THEY SHOULD BE SET HERE.
            % C
            % C     ANY OTHER TESTS OR CALCULATIONS THAT ARE REQUIRED AFTER EVERY
            % C    STEP SHOULD BE INSERTED HERE.
            % C
            % C     IF IDID = 3, Y0 IS SET TO THE CURRENT Y VALUES ON RETURN.
            % C    IF IDID = 2, H IS CONTROLLED TO HIT TOUT (WITHIN ROUNDOFF
            % C    ERROR), AND THEN THE CURRENT Y VALUES ARE PUT IN Y0 ON RETURN.
            % C    FOR ANY OTHER VALUE OF IDID, CONTROL RETURNS TO THE INTEGRATOR
            % C    UNLESS TOUT HAS BEEN REACHED.  THEN INTERPOLATED VALUES OF Y ARE
            % C    COMPUTED AND STORED IN Y0 ON RETURN.
            % C    IF INTERPOLATION IS NOT DESIRED, THE CALL TO INTERP SHOULD BE
            % C    REMOVED AND CONTROL TRANSFERRED TO STATEMENT 500 INSTEAD OF 520.
            % C --------------------------------------------------------------------

            if (NSTEP > MAXSTP)
                KGO = 5;
                KFLAG = 4;
                % C             TOO MUCH WORK
                disp('\n NUMBER OF STEPS EXCEEDS MAXIMUM \n'); %(LOUT,9130)
                IDID = -6;
                %GOTO 70
                MYcase = 70;
                continue
            end

            if (ITOL == 1)
                D = 0;
                AYI=abs(Y(1:N,1));
                YMAX(1:N)=max([YMAX(1:N),AYI]);            
            end


            if  (IDID==3 ||IDID==1)
                %GOTO 70
                MYcase = 70;
                continue
            end

            if (abs(T-TOUT) <= abs(15*UROUND*TOUT))
                % C             EFFECTIVELY WE HAVE HIT TOUT
                IDID = KFLAG;
                T0 = TOUT;
                Y0(1:N)=Y(1:N,1);
                HO = H;
                return
            end

            if (IDID == 2)
                % C        CONTINUING INTEGRATION BUT MUST HIT TOUT EXACTLY
                if (((T+H)-TOUT)*H > 0)
                    % C      WE HAVE ALREADY OVERSHOT THE OUTPUT POINT OR WE WILL DO
                    % C      SO ON THE NEXT STEP
                    if (((T-TOUT)*H >= 0) || (abs(T-TOUT) <= 100*UROUND*HMAX))
                        % C                 HAVE OVERSHOT, SO INTERPOLATE
                        [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0);
                        T0 = TOUT;
                        HO = H;
                        IDID = KFLAG;                         
                        return

                    else
                        % C         WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE
                        % C         SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY'
                        H = (TOUT-T)* (1-4*UROUND);
                        JSTART = -1;
                    end

                end

            elseif ((T-TOUT)*H >= 0)
                % C             HAVE OVERSHOT, SO INTERPOLATE
                [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0);
                IDID = KFLAG;
                HO = H;
                T0 = TOUT;
                return

            end

            %GOTO 20
            MYcase = 20;
            continue

            % C -------------------------------------------------------------------
            % C  ON AN ERROR RETURN FROM THE INTEGRATOR, AN IMMEDIATE RETURN OCCURS
            % C  IF KFLAG = -2, AND RECOVERY ATTEMPTS ARE MADE OTHERWISE.
            % C  H AND HMIN ARE REDUCED BY A FACTOR OF .1 UP TO 10 TIMES
            % C  BEFORE GIVING UP.
            % C --------------------------------------------------------------------


        case 60
            %*********************************************************

            if (NHCUT==10)
                % C             HAVE REDUCED H TEN TIMES
                disp('\n PROBLEM APPEARS UNSOLVABLE WITH GIVEN INPUT \n  HMIN REDUCED BY A FACTOR OF 1.0E10 \n');
                %WRITE (LOUT,9100)
                %GOTO 70
                MYcase = 70;
                continue
            end

            NHCUT = NHCUT + 1;
            HMIN = 0.1*HMIN;
            H = 0.1*H;
            JSTART = -1;

            %GOTO 20
            MYcase = 20;
            continue

        case 70
            %**************************************************************
            if (abs(T-TOUT) > 1000*UROUND)
                Y0(1:N)=Y(1:N,1);           
                T0 = T;
            else
                %               HAVE PASSED TOUT SO INTERPOLATE
                [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0);
                T0 = TOUT;
                IDID = KFLAG;
            end

            HO = H;
            if (KFLAG ~= 0)
                IDID = KFLAG;
            end
            return
    end
end

%========================================================================
%========================================================================


function [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,Y,TOUT,Y0)


Y0(1:N)=Y(1:N,1);

L = JSTART + 2;
S = (TOUT-T)/H;
S1 = 1;
for J = 2:L
    S1 = S1* (S+(J-2))/(J-1);
    for I = 1:N
        Y0(I) = Y0(I) + S1*Y(I,J);
    end
end

%========================================================================
%========================================================================


function [NQ,EL,ELST,TQ,NCOSET,MAXORD] = COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD)
 
% C --------------------------------------------------------------------
% C     COSET IS CALLED BY THE INTEGRATOR AND SETS THE COEFFICIENTS USED
% C     BY THE CONVENTIONAL BACKWARD DIFFERENTIATION SCHEME AND THE
% C     MODIFIED EXTENDED BACKWARD DIFFERENTIATION SCHEME.  THE VECTOR
% C     EL OF LENGTH NQ+1 DETERMINES THE BASIC BDF METHOD WHILE THE VECTOR
% C     ELST OF LENGTH NQ+2 DETERMINES THE MEBDF.  THE VECTOR TQ OF
% C     LENGTH 4 IS INVOLVED IN ADJUSTING THE STEPSIZE IN RELATION TO THE
% C     TRUNCATION ERROR.  ITS VALUES ARE GIVEN BY THE PERTST ARRAY.  THE
% C     VECTORS EL AND TQ BOTH DEPEND ON METH AND NQ.  THE
% C     COEFFICIENTS IN PERTST NEED TO BE GIVEN TO ONLY ABOUT ONE PERCENT
% C     ACCURACY.  THE ORDER IN WHICH THE GROUPS APPEAR BELOW IS:
% C     COEFFICIENTS FOR ORDER NQ-1, COEFFICIENTS FOR ORDER NQ,
% C     COEFFICIENTS FOR ORDER NQ+1.
% C -------------------------------------------------------------------

PERTST(1,1)=1.; PERTST(2,1)=2.; PERTST(3,1)=4.5;
PERTST(4,1)=7.333; PERTST(5,1)=10.42; PERTST(6,1)=13.7;
PERTST(7,1)=17.15; PERTST(8,1)=20.74;
PERTST(1,2)=2.; PERTST(2,2)=4.5; PERTST(3,2)=7.333;
PERTST(4,2)=10.42; PERTST(5,2)=13.7; PERTST(6,2)=17.15;
PERTST(7,2)=20.74; PERTST(8,2)=24.46;
PERTST(1,3)=4.5; PERTST(2,3)=7.333; PERTST(3,3)=10.42;
PERTST(4,3)=13.7; PERTST(5,3)=17.15; PERTST(6,3)=20.74; 
PERTST(7,3)=24.46; PERTST(8,3)=1.;

% C     ..
% C -------------------------------------------------------------------
% C     THE FOLLOWING COEFFICIENTS SHOULD BE DEFINED TO MACHINE ACCURACY.
% C     THEIR DERIVATION IS GIVEN IN REFERENCE 1.
% C -------------------------------------------------------------------

if (NQ > MAXORD)
    MAXORD = NQ;
end

NCOSET = NCOSET + 1;


switch NQ

    case 1
        EL(1) = 1.0E+0;
      ELST(1) = 1.5E+0;
      ELST(3) = -0.5E+0;

    case 2
        EL(1) = 6.6666666666667E-01;
      EL(3) = 3.3333333333333E-01;
      ELST(1) = 9.5652173913043E-01;
      ELST(3) = 2.1739130434782E-01;
      ELST(4) = -1.7391304347826E-01;

    case 3
        EL(1) = 5.4545454545455E-01;
      EL(3) = 4.5454545454545E-01;
      EL(4) = 1.8181818181818E-01;
      ELST(1) = 7.6142131979695E-01;
      ELST(3) = 3.2994923857868E-01;
      ELST(4) = 8.6294416243654E-02;
      ELST(5) = -9.1370558375634E-02;
  
    case 4
        EL(1) = 0.48E+0;
      EL(3) = 0.52E+0;
      EL(4) = 0.28E+0;
      EL(5) = 0.12E+0;
      ELST(1) = 6.5733706517393E-01;
      ELST(3) = 4.0023990403838E-01;
      ELST(4) = 1.5793682526989E-01;
      ELST(5) = 4.4382247101159E-02;
      ELST(6) = -5.7576969212315E-02;

    case 5
        EL(1) = 4.3795620437956E-01;
      EL(3) = 5.62043795620436E-01;
      EL(4) = 3.43065693430656E-01;
      EL(5) = 1.97080291970802E-01;
      EL(6) = 8.75912408759123E-02;
      ELST(1) = 5.9119243917152E-01;
      ELST(3) = 4.4902473356122E-01;
      ELST(4) = 2.1375427307460E-01;
      ELST(5) = 9.0421610027481503E-02;
      ELST(6) = 2.6409276761177E-02;
      ELST(7) = -4.0217172732757E-02;
        
    case 6
        EL(1) = 4.08163265306120E-01;
      EL(3) = 5.91836734693874E-01;
      EL(4) = 3.87755102040813E-01;
      EL(5) = 2.51700680272107E-01;
      EL(6) = 1.49659863945577E-01;
      EL(7) = 6.80272108843534E-02;
      ELST(1) = 5.4475876041119E-01;
      ELST(3) = 4.8525549636077E-01;
      ELST(4) = 2.5789750131312E-01;
      ELST(5) = 1.3133738525800E-01;
      ELST(6) = 5.7677396763462E-02;
      ELST(7) = 1.7258197643881E-02;
      ELST(8) = -3.0014256771967E-02;

    case 7
        EL(1) = 3.85674931129476E-01;
      EL(3) = 6.14325068870521E-01;
      EL(4) = 4.21487603305783E-01;
      EL(5) = 2.9292929292929E-01;
      EL(6) = 1.96510560146923E-01;
      EL(7) = 1.19375573921028E-01;
      EL(8) = 5.50964187327820E-02;
      ELST(1) = 5.0999746293734E-01;
      ELST(3) = 5.1345839935281E-01;
      ELST(4) = 2.9364346131937E-01;
      ELST(5) = 1.6664672120553E-01;
      ELST(6) = 8.8013735242353E-02;
      ELST(7) = 3.9571794884069E-02;
      ELST(8) = 1.2039080338722E-02;
      ELST(9) = -2.3455862290154E-02;

end

TQ(1:3)=PERTST(NQ,1:3);

TQ(4) = 0.5*TQ(2)/NQ;

if(NQ~=1)
    TQ(5)=PERTST(NQ-1,1);
end
%========================================================================
%========================================================================

%PAGE 16
function [Y, N, H, T, UROUND, EPSJAC, CON, MITER, MBND, MASBND, NIND1, ...
    NIND2, NIND3, IER, NRENEW, YMAX, SAVE1, SAVE2, PW, PWCOPY, AM, ...
    WRKSPC, IPIV, ITOL, RTOL, ATOL, NPSET, NJE, NFE, NDEC, IPAR, ...
    RPAR, IERR] = PSET(Y, N, H, T, UROUND, EPSJAC, CON, MITER, MBND, ...
    MASBND, NIND1, NIND2, NIND3, IER, NRENEW, ...
    YMAX, SAVE1, SAVE2, PW, PWCOPY, AM, WRKSPC, ...
    IPIV, ITOL, RTOL, ATOL, NPSET, NJE, NFE, NDEC, ...
    IPAR, RPAR, IERR)
% C -------------------------------------------------------------------
% C   PSET IS CALLED BY STIFF TO COMPUTE AND PROCESS THE MATRIX
% C   M/(H*EL(1)) - J  WHERE J IS AN APPROXIMATION TO THE RELEVANT JACOBIAN
% C   AND M IS THE MASS MATRIX.  THIS MATRIX IS THEN SUBJECTED TO LU
% C   DECOMPOSITION IN PREPARATION FOR LATER SOLUTION OF LINEAR SYSTEMS
% C   OF ALGEBRAIC EQUATIONS WITH LU AS THE COEFFICIENT MATRIX.  THE
% C   MATRIX J IS FOUND BY THE USER-SUPPLIED ROUTINE PDERV IF MITER=1
% C   OR 3 OR BY FINITE DIFFERENCING IF MITER = 2 OR 4.
% C   IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION WITH
% C   PSET USES THE FOLLOWING ..
% C   EPSJAC = DSQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS.
% C *******************************************************************
% C     THE ARGUMENT NRENEW IS USED TO SIGNAL WHETHER OR NOT
% C     WE REQUIRE A NEW JACOBIAN TO BE CALCULATED.
% C
% C        IF NRENEW > 0 THEN WE REQUIRE A NEW J TO BE COMPUTED
% C                  = 0 THEN USE A COPY OF THE LAST J COMPUTED
% C *******************************************************************

NPSET = NPSET + 1;
IMAS = MASBND(1);
MLMAS = MASBND(2);
MUMAS = MASBND(3);

ML = MBND(1);
MU = MBND(2);

MYstop = 0;
MYcase = 1;

while MYstop == 0

    switch MYcase

        case 1 
            %**************************************************************            
            if (NRENEW==0)
                if (MITER<3)
                    PW(1:N^2)=-PWCOPY(1:N^2);                   
                else
                    for I=0:N-1
                        for J = (MBND(1)+1):MBND(4)
                            PW(I*MBND(4)+J) = -PWCOPY(I*MBND(4)+J);
                        end
                    end
                end

                if (IMAS == 0)
                    %GOTO 70                    
                    MYcase = 70;
                    continue
                end


                if (MASBND(4)==N)
                    if(MLMAS==N)
                        PW(1:N^2)=PW(1:N^2)+AM(1:N^2)/CON;                      
                    else
                        FOUR = MUMAS + MLMAS + 1;
                        FIVE = FOUR + MLMAS;
                        for I=1:N
                            for J=(I-MUMAS):(I+MLMAS)
                                if( J<1 || J>N )
                                else
                                    PW((I-1)*N+J) = PW((I-1)*N+J)+ ...
                                        AM((J-1)*FIVE+I-J+FOUR)/CON;
                                end
                            end
                        end
                    end
                else
                    for I=0:N-1
                        for J=(ML+1):(ML+MLMAS+MUMAS+1)
                            PW(I*MBND(4)+J)=PW(I*MBND(4)+J)+ ...
                                AM(I*MASBND(4)+J-MBND(1))/CON;
                        end
                    end
                end

                %GOTO 70                
                MYcase = 70;
                continue

            end

            if (MITER==2||MITER==4)
                %GOTO 30
                MYcase = 30;
                continue
            end

            NJE = NJE + 1;
            
            if ( MITER~=3 )
                 PWCOPY = reshape(PWCOPY,N,N);
                [T,Y,PWCOPY,N,N,IPAR,RPAR,IERR] = PDERV(T,Y,PWCOPY,N, ...
                    N,IPAR,RPAR,IERR);
                 PWCOPY = reshape(PWCOPY,1,N*N);
                PW(1:N^2)=-PWCOPY(1:N^2);              

                if (IMAS==0)
                    %GOTO 70                   
                    MYcase = 70;
                    continue
                end

                if (MLMAS==N)
                    PW(1:N^2)=PW(1:N^2)+AM(1:N^2)/CON;     
                else
                    FOUR = MUMAS + MLMAS + 1;
                    FIVE = FOUR + MLMAS;
                    for I=1:N
                        for J=(I-MUMAS):(I+MLMAS)
                            if (J < 1 || J > N)
                            else
                                PW((I-1)*N+J) = PW((I-1)*N+J) + ...
                                    AM((J-1)*FIVE+I-J+FOUR)/CON;
                            end
                        end
                    end
                end
            else
                
                 PWCOPY = reshape(PWCOPY,MBND(4),N);
                [T,Y,PWCOPY,N,MBND(4),IPAR,RPAR,IERR] = PDERV(T, Y, ...
                    PWCOPY, N, ...
                    MBND(4), IPAR,RPAR,IERR);
                 PWCOPY = reshape(PWCOPY,1,MBND(4)*N);
                     
                for I=0:N-1
                    ITMB=I*MBND(4);
                    for J=(MBND(1)+1):(MBND(4))
                        PW(ITMB+J) = -PWCOPY(ITMB+J-MBND(1));
                      
                    end
                end
               
                if (IMAS == 0)
                    %GOTO 70                    
                    MYcase = 70;
                    continue
                end

                for J = 1:N
                    I1 = max([1, MUMAS+1+1-J]);
                    I2 = min([MLMAS+MUMAS+1, MUMAS+1+N-J]);
                    for I = I1:I2
                        K = (J-1)*MBND(4)+ML+MU-MUMAS;
                        PW(K+I) = PW(K+I)+AM((J-1)*MASBND(4)+I)/CON;
                    end
                end
            end

            %GOTO 70            
            MYcase = 70;
            continue

        case 30
            %**********************************************************

            NJE = NJE + 1;
            for I=1:N
                if (ITOL==2)
                    YMAX(I)=abs(Y(I,1))*RTOL(1)+ATOL(1);
                elseif (ITOL==3)
                    YMAX(I)=zbs(Y(I,1))*RTOL(1)+ATOL(I);
                elseif (ITOL==4)
                    YMAX(I)=abs(Y(I,1))*RTOL(I)+ATOL(1);
                elseif (ITOL==5)
                    YMAX(I)=abs(Y(I,1))*RTOL(I)+ATOL(I);
                end
            end

            D = 0;
                 for I = 1:N
                   D = D + (SAVE2(I)*(YMAX(I)))^2;
                 end
	    
            if (ITOL==1)
                D=D*RTOL(1)^2;
            end

            D = sqrt(D)/N;
            R0 = abs(H)*D*N*(1.0E+03)*UROUND;
            if (R0==0)
                R0=1;
            end

            J1 = 0;

            if (MITER==4)
                %GOTO 51
                MYcase = 51;
                continue
            end
            

            for J = 1:N
                YJ = Y(J,1);
                if(ITOL==1)
                    R=max([EPSJAC*abs(YJ),R0/(YMAX(J)*RTOL(1))]);
                else
                    R=max([EPSJAC*abs(YJ),R0/YMAX(J)]);
                end

                Y(J,1) = Y(J,1) + R;
                D = CON/R;
                [N,T,Y,WRKSPC,IPAR,RPAR,IERR] = F(N,T,Y,WRKSPC,IPAR,...
                    RPAR,IERR);

                for I = 1:N
                    JJKK = I + J1;
                    TEMPRY = (WRKSPC(I)-SAVE2(I));
                    PWCOPY(JJKK) = TEMPRY/R;
                    PW(JJKK) =  -PWCOPY(JJKK);
   
                end

                Y(J,1) = YJ;
                J1 = J1 + N;
            end

            NFE = NFE + N;
            
            if (IMAS==0)
                %GOTO 70
                MYcase = 70;
                continue
            end
           

            if (MLMAS==N)
                PW(1:N^2)=PW(1:N^2)+AM(1:N^2)/CON;
                %GOTO 70
                MYcase = 70;
                continue
            else
                
                FOUR = MUMAS+MLMAS+1;
                FIVE = FOUR + MLMAS;
                for I=1:N
                    for J=(I-MUMAS):(I+MLMAS)
                        if (J<1 || J>N)
                        else
                            PW((I-1)*N+J) = PW((I-1)*N+J) + ...
                                AM((J-1)*FIVE+I-J+FOUR)/CON;
                        end
                    end
                end
            end

            %GOTO 70
            MYcase = 70;
            continue

        case 51
            %****************************************************

            MBA = min([MBND(3),N]);
            for J=1:MBA
                for I=J:MBND(3):N 
                    SAVE1(I) = Y(I,1);
                    YI=Y(I,1);
                    if (ITOL==1)
                        R=max([EPSJAC*abs(YI),R0/(YMAX(I)*RTOL(1))]);
                    else
                        R=max([EPSJAC*abs(YI),R0/YMAX(I)]);
                    end
                    Y(I,1)=Y(I,1)+R;
                end

                [N,T,Y,WRKSPC,IPAR,RPAR,IERR] = F(N,T,Y,WRKSPC,...
                    IPAR,RPAR,IERR);
                
                for JJ=J:MBND(3):N
                    Y(JJ,1)=SAVE1(JJ);
                    YJJ=Y(JJ,1);
                    if (ITOL==1)
                        R=max([EPSJAC*abs(YJJ),R0/(YMAX(JJ)*RTOL(1))]);
                    else
                        R=max([EPSJAC*abs(YJJ),R0/YMAX(JJ)]);
                    end
                    D=CON/R;
                    
                    I1 = max([JJ-MU,1]);
                    I2 = min([JJ+ML,N]);
                    II = JJ*(MBND(4)-1)-ML;

                    for I = I1:I2
                        TEMPRY = WRKSPC(I) - SAVE2(I);
                        PWCOPY(II+I)=TEMPRY/R;
                        PW(II+I) = -PWCOPY(II+I);
                    end
                end
            end

            NFE=NFE+MBND(3);

            if (IMAS==0)
                %GOTO 70
                MYcase = 70;
                continue
            end

            for J = 1:N
                I1 = max([1, MUMAS+1+1-J]);
                I2 = min([MLMAS + MUMAS + 1, MUMAS + 1 + N - J]);
                for I = I1:I2
                    PW((J-1)*MBND(4)+ ML + MU - MUMAS + I) = ...
                        PW((J-1)*MBND(4) + ML + MU - MUMAS + I)+ ...
                        AM((J-1) * MASBND(4)+I)/CON;
                end
            end
    %GOTO 70
    MYcase = 70;
    continue
   

        case 70
            %*************************************************
            if (MITER > 2)
                if(MASBND(1)==0)
                    ML=MBND(1);
                    MU=MBND(2);
                    II=ML+MU+1;
                    for I=1:N
                        PW(II) = PW(II) + 1/CON;
                        II = II + MBND(4);
                    end
                end
                 %RESHAPE IS NEEDED TO PASS PW TO DGBFA.
                 PW = reshape(PW,MBND(4),N);
                        [PW,N,ML,MU,IPIV,IER] = DGBFA(PW,MBND(4),N,ML,MU,IPIV,IER);
                 PW = reshape(PW,1,MBND(4)*N);                
                NDEC = NDEC + 1;
            else
                if(MASBND(1)==0)
                    J = 1;
                    NP1 = N + 1;
                    for I = 1:N
                        PW(J) = PW(J) + 1/CON;
                        J = J + NP1;
                    end
                end                              
                PW = reshape(PW,N,N);                            
                [N,N,PW,IPIV,IER] = DEC(N,N,PW,IPIV,IER);              
                PW = reshape(PW,1,N*N);               
                NDEC = NDEC + 1;
            end            
            return
    end
end

%========================================================================
%========================================================================


function [N,NDIM,A,IP,IER] = DEC(N,NDIM,A,IP,IER)

% C -------------------------------------------------------------------
% C     MATRIX TRIANGULISATION BY GAUSSIAN ELIMINATION
% C     INPUT..
% C     N = ORDER OF MATRIX.
% C     NDIM = DECLARED DIMENSION OF ARRAY A.
% C     A = MATRIX TO BE TRIANGULARISED.
% C     OUTPUT..
% C     A(I,J),  I.LE.J = UPPER TRIANGULAR FACTOR, U.
% C     A(I,J),  I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I-L.
% C     IP(K), K.LT.N = INDEX OF KTH PIVOT ROW.
% C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0.
% C     IER = 0 IF MATRIX IS NON-SINGULAR, OR K IF FOUND TO BE SINGULAR
% C                  AT STAGE K.
% C     USE SOL TO OBTAIN SOLUTION OF LINEAR SYSTEM.
% C     DETERM(A) = IP(N)*A(1,1)*A(2,2)* . . . *A(N,N).
% C     IF IP(N) = 0, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
% C
% C     REFERENCE.
% C     C.B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, C.A.C.M
% C     15 (1972), P.274.
% C     ------------------------------------------------------------------

IER = 0;
IP(N) = 1;

if (N==1)
    %GO TO 70
    %%%%%%%%%%%%%%%%%%%%%%%%%LABEL 70%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = N;
    if (A(N,N)==0)
        %GOTO 80
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LABEL 80%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IER = K;
        IP(N) = 0;
        return
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END LABEL 80%%%%%%%%%%%%%%%%%%%%%%%
    end
    return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END LABEL 70%%%%%%%%%%%%%%%%%%%%%%%
end

NM1 = N - 1;

for K = 1:NM1
    KP1 = K + 1;
    M = K;
    for I = KP1:N
        if (abs(A(I,K)) > abs(A(M,K)))
            M = I;
        end
    end

    IP(K) = M;
    T = A(M,K);

    if (M==K)
    else

        IP(N) = -IP(N);
        A(M,K) = A(K,K);
        A(K,K) = T;
    end

    if (T==0)
        %GO TO 80
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LABEL 80%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IER = K;
        IP(N) = 0;
        return
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END LABEL 80%%%%%%%%%%%%%%%%%%%%%%%
    end

    T = 1/T;

    A(KP1:N,K)=-A(KP1:N,K)*T;

    for J = KP1:N
        T = A(M,J);
        A(M,J) = A(K,J);
        A(K,J) = T;
        if (T==0)
        else
            for I = KP1:N;
                A(I,J) = A(I,J) + A(I,K)*T;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%LABEL 70%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = N;
if (A(N,N)==0)
    %GO TO 80
    %%%%%%%%%%%%%%%%%%%%%%%%%LABEL 80%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IER = K;
    IP(N) = 0;
    return
    %%%%%%%%%%%%%%%%%%%%%%%END LABEL 80%%%%%%%%%%%%%%%%%%%%%%%
end
return
%%%%%%%%%%%%%%%%%%%%%%%END LABEL 70%%%%%%%%%%%%%%%%%%%%%%%

%========================================================================
%========================================================================


function [N,NDIM,A,B,IP] = SOL(N,NDIM,A,B,IP)

% C     ------------------------------------------------------------------
% C     SOLUTION OF LINEAR SYSTEM, A*X = B.
% C     INPUT ..
% C     N = ORDER OF MATRIX.
% C     NDIM = DECLARED DIMENSION OF MATRIX A.
% C     A = TRIANGULARISED MATRIX OBTAINED FROM DEC.
% C     B = RIGHT HAND SIDE VECTOR.
% C     IP = PIVOT VECTOR OBTAINED FROM DEC.
% C     DO NOT USE IF DEC HAS SET IER .NE. 0
% C     OUTPUT..
% C     B = SOLUTION VECTOR, X.
% C     ------------------------------------------------------------------

if (N==1)
    %GO TO 50
    %%%%%%%%%%%%%%%%%%%%LABEL 50%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B(1) = B(1)/A(1,1);
    return
    %%%%%%%%%%%%%%%%%%%%END LABEL 50%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

NM1 = N - 1;

for K = 1:NM1
    KP1 = K + 1;
    M = IP(K);
    T = B(M);
    B(M) = B(K);
    B(K) = T;
 for I = KP1:N
        B(I) = B(I) + A(I,K)*T;
    end
end

for KB = 1:NM1
    KM1 = N - KB;
    K = KM1 + 1;
    B(K) = B(K)/A(K,K);
    T = -B(K);
    for I = 1:KM1
        B(I) = B(I) + A(I,K)*T;
    end
end

%%%%%%%%%%%%%%%%%%%%LABEL 50%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B(1) = B(1)/A(1,1);
return
%%%%%%%%%%%%%%%%%%%%END LABEL 50%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%========================================================================
%========================================================================

function [abd,n,ml,mu,ipvt,info] = DGBFA(abd,lda,n,ml,mu,ipvt,info)

%       integer lda,n,ml,mu,ipvt(1),info
%       double precision abd(lda,1)
% c
% c     dgbfa factors a double precision band matrix by elimination.
% c
% c     dgbfa is usually called by dgbco, but it can be called
% c     directly with a saving in time if  rcond  is not needed.
% c
% c     on entry
% c
% c        abd     double precision(lda, n)
% c                contains the matrix in band storage.  the columns
% c                of the matrix are stored in the columns of  abd  and
% c                the diagonals of the matrix are stored in rows
% c                ml+1 through 2*ml+mu+1 of  abd .
% c                see the comments below for details.
% c
% c        lda     integer
% c                the leading dimension of the array  abd .
% c                lda must be .ge. 2*ml + mu + 1 .
% c
% c        n       integer
% c                the order of the original matrix.
% c
% c        ml      integer
% c                number of diagonals below the main diagonal.
% c                0 .le. ml .lt. n .
% c
% c        mu      integer
% c                number of diagonals above the main diagonal.
% c                0 .le. mu .lt. n .
% c                more efficient if  ml .le. mu .
% c     on return
% c
% c        abd     an upper triangular matrix in band storage and
% c                the multipliers which were used to obtain it.
% c                the factorization can be written  a = l*u  where
% c                l  is a product of permutation and unit lower
% c                triangular matrices and  u  is upper triangular.
% c
% c        ipvt    integer(n)
% c                an integer vector of pivot indices.
% c
% c        info    integer
% c                = 0  normal value.
% c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
% c                     condition for this subroutine, but it does
% c                     indicate that  will divide by zero if
% c                     called.  use  rcond  in dgbco for a reliable
% c                     indication of singularity.
% c
% c     band storage
% c
% c           if  a  is a band matrix, the following program segment
% c           will set up the input.
% c
% c                   ml = (band width below the diagonal)
% c                   mu = (band width above the diagonal)
% c                   m = ml + mu + 1
% c                   do 20 j = 1, n
% c                      i1 = max0(1, j-mu)
% c                      i2 = min0(n, j+ml)
% c                      do 10 i = i1, i2
% c                         k = i - j + m
% c                         abd(k,j)  = a(i,j)
% c                10    continue
% c                20 continue
% c
% c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
% c           in addition, the first  ml  rows in  abd  are used for
% c           elements generated during the triangularization.
% c           the total number of rows needed in  abd  is  2*ml+mu+1 .
% c           the  ml+mu by ml+mu  upper left triangle and the
% c           ml by ml  lower right triangle are not referenced.
% c
% c     linpack. this version dated 08/14/78 .
% c     cleve moler, university of new mexico, argonne national lab.
% c
% c     subroutines and functions
% c
% c     blas daxpy,dscal,idamax
% c     fortran max0,min0
% c

m = ml + mu + 1;

info = 0;

% c zero initial fill-in columns
% c
j0 = mu + 2;
j1 = min([n,m]) - 1;


if (j1 < j0)
else
    for jz = j0:j1
        i0 = m + 1 - jz;
        abd(i0:ml,jz)=0;       
    end
end

jz = j1;
ju = 0;
% c
% c gaussian elimination with partial pivoting
% c
nm1 = n - 1;

if (nm1 < 1)
    %GOTO 130
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LABEL 130%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ipvt(n) = n;
    if (abd(m,n) == 0)
        info = n;
    end
    return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END LABEL 130%%%%%%%%%%%%%%%%%%%%%%%%%%
end

for k = 1: nm1
  
    kp1 = k + 1;
    % c
    % c     zero next fill-in column
    % c
    jz = jz + 1;
   
    if (jz > n)
    else

        if (ml < 1)
        else
            abd(1:ml,jz)=0; 
        end
    end
    
    % c
    % c     find l = pivot index
    % c

    lm = min([ml,n-k]);  
        SIZE_ABD = size(abd);
        abd = reshape(abd,1,SIZE_ABD(1)*SIZE_ABD(2));
        
        for IDAMAX_COUNT = 0:lm
            ABD2(IDAMAX_COUNT+1) = abd((k-1)*SIZE_ABD(1)+m+IDAMAX_COUNT);
        end
       
        l = IDAMAX(lm+1,ABD2,1) + m - 1;

        abd = reshape(abd,SIZE_ABD(1),SIZE_ABD(2));

    ipvt(k) = l + k - m;
  
    % c
    % c     zero pivot implies this column already triangularized
    % c

    if (abd(l,k) == 0)
        %go to 100
    else
        % c
        % c     interchange if necessary
        % c
        
        if (l==m)
            %go to 60
        else
            t = abd(l,k);
            abd(l,k) = abd(m,k);
            abd(m,k) = t;
            %60       continue
        end
        % c
        % c         compute multipliers
        % c
        t = -1/abd(m,k);

        SIZE_ABD = size(abd);
        abd = reshape(abd,1,SIZE_ABD(1)*SIZE_ABD(2));


        for DRUN = 0:(lm-1)
        abd((k-1)*SIZE_ABD(1)+(m+1)+DRUN) = t* abd((k-1)*SIZE_ABD(1)+(m+1)+DRUN);
        end

        abd = reshape(abd,SIZE_ABD(1),SIZE_ABD(2));
        % c
        % c         row elimination with column indexing
        % c
        ju = min([max([ju,mu+ipvt(k)]),n]);
        mm = m;

        if (ju < kp1)
          % go to 90
        else
          for j = kp1:ju
            l = l - 1;
            mm = mm - 1;
            t = abd(l,j);
            if (l==mm)
              %go to 70
            else
              abd(l,j) = abd(mm,j);
              abd(mm,j) = t;
              %70              continue
            end
            %CALL DAXPY
                        SIZE_ABD = size(abd);
                        abd =reshape(abd,1,SIZE_ABD(1)*SIZE_ABD(2));

                        %" constant times a vector plus a vector "

                        for DAXPY_COUNT = 0:(lm-1)
                           abd((j-1)*SIZE_ABD(1)+mm+1+ DAXPY_COUNT)  = abd((j-1)*SIZE_ABD(1)+mm+1+ DAXPY_COUNT)+t*abd((k-1)*SIZE_ABD(1)+ m+1 + DAXPY_COUNT);
                        end

                        abd = reshape(abd,SIZE_ABD(1),SIZE_ABD(2));       
          end
          %90          continue
        end
	
        %GOTO 110
        continue %110 LABEL THAT SENDS YOU BACK TO THE START OF THE LOOP.
		%100      continue
    end
    info = k;
end

%LABEL 130
ipvt(n) = n;
if (abd(m,n) == 0)
  info = n;
end
return

%========================================================================
%========================================================================
function idamax = IDAMAX(n,dx,incx)
% c
% c     finds the index of element having max. absolute value.
% c     jack dongarra, linpack, 3/11/78.
% c     modified to correct problem with negative increment, 8/21/90.
% c

idamax = 0;

if (n<1)
    return
end

idamax = 1;

if (n==1)
    return
end

if (incx==1)
    %go to 20
else
    % c
    % c     code for increment not equal to 1
    % c
    ix = 1;

    if(incx<0)
        ix = (-n+1)*incx + 1;
    end

    dmax = abs(dx(ix));
    ix = ix + incx;

    for i = 2:n
        if(abs(dx(ix))<=dmax)
            %go to 5
        else
            idamax = i;
            dmax = abs(dx(ix));
        end
        ix = ix + incx;
    end
    return
    % c
    % c     code for increment equal to 1
    % c
end
dmax = abs(dx(1)); %20
for i = 2:n
    if (abs(dx(i))<=dmax)
        %go to 30
    else
        idamax = i;
        dmax = abs(dx(i));
    end
end
 %disp('LEFT IDAMAX')
return

%========================================================================
%========================================================================


function [abd,lda,n,ml,mu,ipvt,b] = DGBSL(abd,lda,n,ml,mu,ipvt,b,job)

% c     dgbsl solves the double precision band system
% c     a * x = b  or  trans(a) * x = b
% c     using the factors computed by dgbco or dgbfa.
% c
% c     on entry
% c
% c        abd     double precision(lda, n)
% c                the output from dgbco or dgbfa.
% c
% c        lda     integer
% c                the leading dimension of the array  abd .
% c
% c        n       integer
% c                the order of the original matrix.
% c
% c        ml      integer
% c                number of diagonals below the main diagonal.
% c
% c        mu      integer
% c                number of diagonals above the main diagonal.
% c
% c        ipvt    integer(n)
% c                the pivot vector from dgbco or dgbfa.
% c
% c        b       double precision(n)
% c                the right hand side vector.
% c
% c        job     integer
% c                = 0         to solve  a*x = b ,
% c                = nonzero   to solve  trans(a)*x = b , where
% c                            trans(a)  is the transpose.
% c
% c     on return
% c
% c        b       the solution vector  x .
% c
% c     error condition
% c
% c        a division by zero will occur if the input factor contains a
% c        zero on the diagonal.  technically this indicates singularity
% c        but it is often caused by improper arguments or improper
% c        setting of lda .  it will not occur if the subroutines are
% c        called correctly and if dgbco has set rcond .gt. 0.0
% c        or dgbfa has set info .eq. 0 .
% c
% c     to compute  inverse(a) * c  where  c  is a matrix
% c     with  p  columns
% c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
% c           if (rcond is too small) go to ...
% c           do 10 j = 1, p
% c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
% c        10 continue
% c
% c     linpack. this version dated 08/14/78 .
% c     cleve moler, university of new mexico, argonne national lab.
% c
% c     subroutines and functions
% c
% c     blas daxpy,ddot
% c     fortran min0
% c
% c     internal variables
% c
%       double precision ddot,t
%       integer k,kb,l,la,lb,lm,m,nm1
% c
m = mu + ml + 1;
nm1 = n - 1;

if (job ~= 0)
else

    % c
    % c     job = 0 , solve  a * x = b
    % c     first solve l*y = b
    % c

    if (ml==0)
    else
        if (nm1<1)
        else
            for k = 1:nm1
                lm = min([ml,n-k]);
                l = ipvt(k);
                t = b(l);
                if (l==k)
                else
                    b(l) = b(k);
                    b(k) = t;
                end
                    %CALL DAXPY
                    SIZE_ABD = size(abd);
                    
                    abd =reshape(abd,1,SIZE_ABD(1)*SIZE_ABD(2));
                    
                    %" Vecter times constant times vecter "
                    
                    for DAXPY_COUNT = 0:(lm-1)
                       b(k+1+DAXPY_COUNT) = b(k+1+DAXPY_COUNT)+t*abd((k-1)*SIZE_ABD(1)+m+1 + DAXPY_COUNT);
                    end
                  
                    abd = reshape(abd,SIZE_ABD(1),SIZE_ABD(2));
                    
            end

        end
    end

    % c
    % c     now solve  u*x = y
    % c

    for kb = 1:n
        k = n + 1 - kb;
        b(k) = b(k)/abd(m,k);
        lm = min([k,m]) - 1;
        la = m - lm;
        lb = k - lm;
        t = -b(k);
        
        %CALL DAXPY
                    
                    SIZE_ABD = size(abd);
                    
                    abd =reshape(abd,1,SIZE_ABD(1)*SIZE_ABD(2));
                    
                    %" Vecter times constant times vecter "
                    
                    for DAXPY_COUNT = 0:(lm-1)
                       b(lb+DAXPY_COUNT) =b(lb+DAXPY_COUNT)+ t*abd((k-1)*SIZE_ABD(1)+la + DAXPY_COUNT);
                    end
                
                    abd = reshape(abd,SIZE_ABD(1),SIZE_ABD(2));    
    end
    return
end
% c
% c job = nonzero, solve  trans(a) * x = b
% c first solve  trans(u)*y = b
% c
for k = 1:n
    lm = min([k,m]) - 1;
    la = m - lm;
    lb = k - lm;
    %WARNING: THIS SECTION OF THE CODE HAS NOT BEEN TESTED.
    
    %t = DDOT(lm,abd(la,k),1,b(lb),1);
    SIZE_ABD = size(abd);
    abd =reshape(abd,1,SIZE_ABD(1)*SIZE_ABD(2));
    abd_temp(1:lm) =abd((la-1)*SIZE_ABD(1)+k:(la-1)*SIZE_ABD(1)+k+lm-1);
    abd = reshape(abd,SIZE_ABD(1),SIZE_ABD(2));
    b_temp(1:lm) = b(lb:lb+lm-1);
    
    t = dot(b_temp,abd_temp);
    
    b(k) = (b(k) - t)/abd(m,k);
end
% c
% c now solve trans(l)*x = y
% c
if (ml==0)
    return
end

if (nm1<1)
    return
end

for kb = 1:nm1
    k = n - kb;
    lm = min([ml,n-k]);
    %WARNING: THIS SECTION OF THE CODE HAS NOT BEEN TESTED.
    %b(k) = b(k) + DDOT(lm,abd(m+1,k),1,b(k+1),1);
    SIZE_ABD = size(abd);
    abd =reshape(abd,1,SIZE_ABD(1)*SIZE_ABD(2));
    abd_temp(1:lm) =abd((m)*SIZE_ABD(1)+k:(m)*SIZE_ABD(1)+k+m);
    abd = reshape(abd,SIZE_ABD(1),SIZE_ABD(2));
    b_temp(1:k+1) = b(k+1:k+lm);
    
    t = dot(b_temp,abd_temp);
    
    l = ipvt(k);
    if (l==k)
    else
        t = b(l);
        b(l) = b(k);
        b(k) = t;
    end
end
return

%========================================================================
%========================================================================

function ddot = DDOT(n,dx,incx,dy,incy)
  %disp('now in DDOT')
% c
% c     forms the dot product of two vectors.
% c     uses unrolled loops for increments equal to one.
% c     jack dongarra, linpack, 3/11/78.
% c
%       double precision dx(1),dy(1),dtemp
%       integer i,incx,incy,ix,iy,m,mp1,n
% c
ddot = 0;
dtemp = 0;
if(n<=0)
    return
end

MYstop = 0;
MYcase = 1;

while (MYstop == 0)

    switch MYcase

        case 1

            if ((incx==1)&&(incy==1))
                %go to 20
                MYcase = 20;
                continue
            end

            % c
            % c         code for unequal increments or equal increments
            % c         not equal to 1
            % c
            ix = 1;
            iy = 1;
            if(incx<0)
                ix = (-n+1)*incx + 1;
            end

            if(incy<0)
                iy = (-n+1)*incy + 1;
            end

            for i = 1:n
                dtemp = dtemp + dx(ix)*dy(iy);
                ix = ix + incx;
                iy = iy + incy;
            end
            ddot = dtemp;
            return
            % c
            % c         code for both increments equal to 1
            % c
            % c
            % c         clean-up loop
            % c

        case 20
            % HERE THE FORTRAN MOD IS USED, WHERE MOD(A,B) = A-floor(A/B)*B.
            m = n-floor(n/5)*5;
            %m = mod(n,5)
            if(m==0)
                %go to 40
                MYcase = 40;
                continue
            end

	    dtemp=dtemp+dx(1:m)*dy(1:m);
            if( n<5 )
                %go to 60
                MYcase = 60;
                continue
            end

            %go to 40
            MYcase = 40;
            continue

        case 40

            mp1 = m + 1;

            for i = mp1:5:n
                dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + ...
                    dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + ...
                    dx(i + 4)*dy(i + 4);
            end

            %go to 60
            MYcase = 60;
            continue

        case 60

            ddot = dtemp;
            return

    end
end

%========================================================================
%========================================================================

function [N,TQ,EDN,E,EUP,BND,EDDN] = FERRORS(N,TQ,EDN,E,EUP,BND,EDDN)
%PREVIOUSLY NAMED ERRORS()

% C     ***************************************************
% C
% C     THIS ROUTINE CALCULATES ERRORS USED IN TESTS
% C     IN STIFF .
% C
% C     ***************************************************

SQHOL = N;
EDN = TQ(1)*TQ(1)*SQHOL;

% C
% C     ** ERROR ASSOCIATED WITH  METHOD OF ORDER ONE LOWER.
% C

E = TQ(2)*TQ(2)*SQHOL;

% C
% C     ** ERROR ASSOCIATED WITH PRESENT ORDER
% C

EUP = TQ(3)*TQ(3)*SQHOL;
      
% C
% C     ** ERROR ASSOCIATED WITH HIGHER ORDER METHOD
% C

BND = TQ(4)*TQ(4)*SQHOL*0.5;

EDDN=TQ(5)*TQ(5)*SQHOL;

% C
% C     ** ERROR ASSOCIATED WITH METHOD OF ORDER TWO LOWER.

%========================================================================
%========================================================================

function [T, H, Y, L, N, YPRIME, NFE, IPAR, RPAR, IERR] = PRDICT( T, H, ...
    Y, L, N, YPRIME, NFE, IPAR, RPAR, IERR)

% C **********************************************************************
% C     PREDICTS A VALUE FOR Y AT (T+H) GIVEN THE HISTORY ARRAY AT T
% C     THEN EVALUATES THE DERIVATIVE AT THIS POINT, THE RESULT OF THIS
% C     EVALUATION BEING STORED IN YPRIME()
% C **********************************************************************

for J2 = 2:L
    for I = 1:N
        Y(I,1) = Y(I,1) + Y(I,J2);
    end
end

T = T + H;

[N,T,Y,YPRIME,IPAR,RPAR,IERR] = F(N,T,Y,YPRIME,IPAR,RPAR,IERR);

NFE = NFE + 1;

%========================================================================
%========================================================================

function [QQQ, Y, N, T, HBETA, ERRBND, ARH, CRATE, TCRATE, M, WORKED, YMAX, ...
    ERROR, SAVE1, SAVE2, SCALE, PW, MF, MBND, AM, MASBND, NIND1, ...
    NIND2, NIND3, IPIV, ITOL, RTOL, ATOL, IPAR, RPAR, HUSED, NBSOL, ...
    NFE, NQUSED,IERR] = ITRAT2(QQQ, Y, N, T, HBETA, ERRBND, ARH, ...
    CRATE, TCRATE, M, WORKED, YMAX, ERROR, ...
    SAVE1, SAVE2, SCALE, PW, MF, MBND, ...
    AM, MASBND, NIND1, NIND2, NIND3, IPIV, ...
    LMB, ITOL, RTOL, ATOL, IPAR, RPAR, ...
    HUSED, NBSOL, NFE, NQUSED, IERR)

ZERO = 0;
% C     ..


IMAS = MASBND(1);
MLMAS = MASBND(2);
MUMAS = MASBND(3);

for I=1:N
    AYI = abs(Y(I,1));
    if (ITOL==1)
        SCALE(I) = YMAX(I);
    elseif (ITOL==2)
        SCALE(I) = RTOL(1)*AYI + ATOL(1);
    elseif (ITOL==3)
        SCALE(I) = RTOL(1)*AYI + ATOL(I);
    elseif (ITOL == 4)
        SCALE(I) = RTOL(I)*AYI + ATOL(1);
    elseif (ITOL == 5)
        SCALE(I) = RTOL(I)*AYI + ATOL(I);
    end
end


if (NIND2~=0)
    for I = (NIND1+1):(NIND2+NIND1)
        SCALE(I)=SCALE(I)/HUSED;
    end
end

if (NIND3~=0)
    for I = (NIND1 +NIND2 + 1):(NIND3+NIND2+NIND1)
        SCALE(I)=SCALE(I)/(HUSED^2);
    end
end

if (LMB==1)
    %GOTO 25
else
    if (IMAS==0)
        SAVE1(1:N)=(-SAVE1(1:N) + HBETA*SAVE2(1:N) +ARH(1:N))/QQQ;
    else

        SAVE1(1:N)=(-SAVE1(1:N)+ARH(1:N));
        if ((MF <= 22) && (MLMAS == N))
            for I=1:N
                SAVE2(I) = SAVE2(I)*HBETA;
                for J=1:N
                    SAVE2(I) = SAVE2(I) + AM((J-1)*N+I)*SAVE1(J);
                end
            end
        else
            C1 = MLMAS + MUMAS;
            C2 = C1 + 1;
            for I=1:N
                SAVE2(I) = SAVE2(I)*HBETA;
                for J = -MLMAS:MUMAS
                    ISUM = I + J;
                    if ((ISUM<1) ||(ISUM > N))
                    else
                        KK = C2*I+C1*J-MLMAS;
                        SAVE2(I) = SAVE2(I) + AM(KK)*SAVE1(I+J);
                    end
                end
            end
        end

        SAVE1(1:N)=SAVE2(1:N)/QQQ;
    end

    if (MF >= 23)
        PW = reshape(PW,MBND(4),N);
        [PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE1] = DGBSL(PW, MBND(4), ...
            N, MBND(1), MBND(2), IPIV, SAVE1, 0);
        PW = reshape(PW,1,MBND(4)*N);
        NBSOL = NBSOL + 1;
    else

        PW = reshape(PW,N,N);
        [N,N,PW,SAVE1,IPIV] = SOL(N,N,PW,SAVE1,IPIV);
        PW = reshape(PW,1,N*N);

        NBSOL = NBSOL + 1;
    end

    D = ZERO;

    for I = 1:N
        ERROR(I) = ERROR(I) + SAVE1(I);
        D = D + (SAVE1(I)/(SCALE(I)))^2;
        SAVE1(I) = Y(I,1) + ERROR(I);
    end
    
    if (ITOL==1)
        D=D/(RTOL(1)^2);
    end

    TCRATE = TCRATE + CRATE;
    D1 = D;
    M = 1;
    [N,T,SAVE1,SAVE2,IPAR,RPAR,IERR] = F(N,T,SAVE1,SAVE2,IPAR,RPAR,IERR);

    NFE = NFE + 1;

end       %25   CONTINUE

%   WORKED = .TRUE.
WORKED = 1;
MYstop = 0;

while (MYstop  == 0)

    if (IMAS == 0)
        
        SAVE1(1:N)=(-SAVE1(1:N)+HBETA*SAVE2(1:N)+ARH(1:N))/QQQ;
 
    else
        SAVE1(1:N) = (-SAVE1(1:N) + ARH(1:N));


        if ((MF <= 22)&&(MLMAS==N))
            for I=1:N
                SAVE2(I) = SAVE2(I)*HBETA;
                for J=1:N
                    SAVE2(I) = SAVE2(I) + AM((J-1)*N+I)*SAVE1(J);
                end
            end
        else
            C1 = MLMAS + MUMAS;
            C2 = C1 + 1;
            for I=1:N
                SAVE2(I) = SAVE2(I)*HBETA;
                for J=(-MLMAS):(MUMAS)
                    ISUM = I+J;
                    if ((ISUM<1)||(ISUM>N))
                    else
                        KK = C2*I + J*C1-MLMAS;
                        SAVE2(I) = SAVE2(I) + AM(KK)*SAVE1(I+J);
                    end
                end
            end
        end
        SAVE1(1:N)=SAVE2(1:N)/QQQ;

    end
    % C
    % C IF WE ARE HERE THEN PARTIALS ARE O.K.
    % C

    if ( MF>=23)

        PW = reshape(PW,MBND(4),N);
        [PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE1] = DGBSL(PW,MBND(4),N, ...
            MBND(1), MBND(2), IPIV, SAVE1, 0);
        PW = reshape(PW,1,MBND(4)*N);

        NBSOL=NBSOL + 1;
    else
 

               PW = reshape(PW,N,N);   
               [N,N,PW,SAVE1,IPIV] = SOL(N,N,PW,SAVE1,IPIV);  
               PW = reshape(PW,1,N*N); 
        NBSOL = NBSOL + 1;
    end
    % C
    % C WE NOW CALCULATE A WEIGHTED RMS TYPE NORM
    % C
    D = ZERO;
 
    for I = 1:N
        ERROR(I) = ERROR(I) + SAVE1(I);
        D = D + (SAVE1(I)/(SCALE(I)))^2;
        SAVE1(I) = Y(I,1) + ERROR(I);
    end


    if (ITOL== 1)
        D=D/(RTOL(1)^2);
    end
    % C -------------------------------------------------------------------
    % C     TEST FOR CONVERGENCE.  IF M.GT.0 , AN ESTIMATE OF THE CONVERGENCE
    % C     RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
    % C -------------------------------------------------------------------

    if (M~=0)
        if (D1~=ZERO)
            CRATE = max([0.9*CRATE,D/D1]);
        end
    end

    TCRATE = TCRATE + CRATE;

    if ((D*min([1,2*CRATE])) < ERRBND/NQUSED)
        return
    end

    if (M~=0)
        if (D > D1)
            %WORKED = .FALSE.
            WORKED = 0;
            return
        end
    end
    D1 = D;

    if (M==4)
        %WORKED = .FALSE.
        WORKED = 0;
        return
    end

    M = M + 1;
    [N,T,SAVE1,SAVE2,IPAR,RPAR,IERR] = F(N,T,SAVE1,SAVE2,IPAR,RPAR,IERR);
    NFE = NFE + 1;

end

%========================================================================
%========================================================================

%PAGE 34
function [H, HMAX, HMIN, JSTART, KFLAG, MF, MBND, MASBND, NIND1, NIND2, ...
    NIND3, T, TOUT, TEND, Y, N, YMAX, ERROR, SAVE1, SAVE2, SCALE, ...
    PW, PWCOPY, AM, YHOLD, YNHOLD, ARH, IPIV, MAXDER, ITOL, ...
    RTOL, ATOL, RPAR, IPAR, NQUSED, NSTEP, NFAIL, NFE, NJE, NDEC, ...
    NBSOL, NPSET, NCOSET, MAXORD, MAXSTP, UROUND, EPSJAC, HUSED, ...
    IERR] = STIFF(H, HMAX, HMIN, JSTART, KFLAG, MF, MBND, MASBND, ...
    NIND1, NIND2, NIND3, T, TOUT, TEND, Y, N, YMAX, ...
    ERROR, SAVE1, SAVE2, SCALE, PW, PWCOPY, AM, YHOLD, ...
    YNHOLD, ARH, IPIV, MAXDER, ITOL, RTOL, ATOL, ...
    RPAR, IPAR, NQUSED, NSTEP, NFAIL, NFE, NJE, NDEC, ...
    NBSOL, NPSET, NCOSET, MAXORD, MAXSTP, UROUND, ...
    EPSJAC, HUSED, IERR)


% C     ------------------------------------------------------------------
% C     THE SUBROUTINE STIFF PERFORMS ONE STEP OF THE INTEGRATION OF AN
% C     INITIAL VALUE PROBLEM FOR A SYSTEM OF ORDINARY DIFFERENTIAL
% C     EQUATIONS OR LINEARLY IMPLICIT DIFFERENTIAL ALGEBRAIC EQUATIONS.
% C     COMMUNICATION WITH STIFF IS DONE WITH THE FOLLOWING VARIABLES..
% C     Y      AN N BY LMAX+3 ARRAY CONTAINING THE DEPENDENT VARIABLES
% C              AND THEIR BACKWARD DIFFERENCES.  MAXDER (=LMAX-1) IS THE
% C              MAXIMUM ORDER AVAILABLE.  SEE SUBROUTINE COSET.
% C              Y(I,J+1) CONTAINS THE JTH BACKWARD DIFFERENCE OF Y(I)
% C     T      THE INDEPENDENT VARIABLE. T IS UPDATED ON EACH STEP TAKEN.
% C     H      THE STEPSIZE TO BE ATTEMPTED ON THE NEXT STEP.
% C              H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING
% C              THE PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE BUT
% C              ITS SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM.
% C     HMIN   THE MINIMUM AND MAXIMUM ABSOLUTE VALUE OF THE STEPSIZE
% C     HMAX   TO BE USED FOR THE STEP.  THESE MAY BE CHANGED AT ANY
% C              TIME BUT WILL NOT TAKE EFFECT UNTIL THE NEXT H CHANGE.
% C     RTOL,ATOL  THE ERROR BOUNDS. SEE DESCRIPTION IN OVDRIV.
% C     N      THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS.
% C     MF     THE METHOD FLAG.  MUST BE SET TO 21,22,23 OR 24 AT PRESENT
% C     KFLAG  A COMPLETION FLAG WITH THE FOLLOWING MEANINGS..
% C                  0  THE STEP WAS SUCCESSFUL
% C                 -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED
% C                       WITH ABS(H) = HMIN.
% C                 -2  THE REQUESTED ERROR IS SMALLER THAN CAN
% C                       BE HANDLED FOR THIS PROBLEM.
% C                 -3  CORRECTOR CONVERGENCE COULD NOT BE
% C                       ACHIEVED FOR ABS(H)=HMIN.
% C            ON A RETURN WITH KFLAG NEGATIVE, THE VALUES OF T AND
% C            THE Y ARRAY ARE AS AT THE BEGINNING OF THE LAST
% C            STEP ATTEMPTED, AND H IS THE LAST STEP SIZE ATTEMPTED.
% C     JSTART  AN INTEGER USED ON INPUT AND OUTPUT.
% C          ON INPUT IT HAS THE FOLLOWING VALUES AND MEANINGS..
% C              0  PERFORM THE FIRST STEP.
% C          .GT.0  TAKE A NEW STEP CONTINUING FROM THE LAST
% C          .LT.0  TAKE THE NEXT STEP WITH A NEW VALUE OF H OR N.
% C          ON EXIT, JSTART IS NQUSED, THE ORDER OF THE METHOD LAST USED.
% C     YMAX     AN ARRAY OF N ELEMENTS WITH WHICH THE ESTIMATED LOCAL
% C              ERRORS IN Y ARE COMPARED
% C     ERROR    AN ARRAY OF N ELEMENTS.
% C     SAVE1,2  TWO ARRAYS OF WORKING SPACE BOTH OF LENGTH N.
% C     PW       A BLOCK OF LOCATIONS USED FOR PARTIAL DERIVATIVES
% C     IPIV     AN INTEGER ARRAY OF LENGTH N USED FOR PIVOT INFORMATION.
% C     JNEWIM   IS TO INDICATE IF PRESENT ITERATION MATRIX
% C                WAS FORMED USING A NEW J OR OLD J.
% C     JSNOLD   KEEPS TRACK OF NO. OF STEPS TAKEN WITH
% C                PRESENT ITERATION MATRIX (BE IT FORMED BY
% C                A NEW J OR NOT).
% C     AVNEWJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION
% C                MATRIX WAS FORMED BY A NEW J.
% C     AVOLDJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION
% C                MATRIX WAS FORMED BY AN OLD J.
% C     NRENEW   FLAG THAT IS USED IN COMMUNICATION WITH SUBROUTINE PSET.
% C                IF  NRENEW > 0  THEN FORM A NEW JACOBIAN BEFORE
% C                                COMPUTING THE COEFFICIENT MATRIX FOR
% C                                THE NEWTON-RAPHSON ITERATION
% C                           = 0  FORM THE COEFFICIENT MATRIX USING A
% C                                COPY OF AN OLD JACOBIAN
% C     NEWPAR   FLAG USED IN THIS SUBROUTINE TO INDICATE IF A JACOBIAN
% C              HAS BEEN EVALUATED FOR THE CURRENT STEP
% C **********************************************************************
global MYtagSTIFF;
global STPSZE HSTPSZ;
persistent MFOLD HOLD LMAX EDN EUP BND EDDN EL TQ ELST E IWEVAL;
persistent NEWPAR NRENEW JSINUP JSNOLD IDOUB JCHANG L NQ MEQC1 QI QQQ ;
persistent MEQC2 MQ1TMP MQ2TMP ISAMP RH RMAX TOLD CRATE1 CRATE2 IER ;
persistent TCRAT1 TCRAT2 AVNEW2 AVOLD2 AVNEWJ AVOLDJ MITER IBND UPBND;
persistent CFAIL KFAIL RC JNEWIM VTOL SAMPLE IEMB OLDLO IREDO OVRIDE;
persistent FFAIL PLLFAL

%SOME INITIALISATION.
WORKED = 0;
FINISH = 0;
if (MYtagSTIFF == 1)
    FFAIL = 0;
    PLLFAL = 0;
    NQ = 0;
    EL = zeros(1,10);
    ELST = zeros(1,10);
    TQ = zeros(1,5);
    OVRIDE = 0;
    MYtagSTIFF = 0;
    EL(2)=1;ELST(2)=1;OLDLO=1;
end

% C     ..
% C     .. DATA STATEMENTS ..
ZERO=0;ONE=1;

% C     ..

MYstop = 0;
MYcase = 6000;
while (MYstop == 0)

    switch MYcase

        case 6000 %********************************************************
           
            TOLD = T;
            KFLAG = 0;
            IMAS = MASBND(1);
            MUMAS = MASBND(3);
            MLMAS = MASBND(2);
            LDMAS  = MASBND(4);
            
            if (JSTART > 0)
                %GO TO 60
                MYcase = 60;
                continue
            end
            
            if (JSTART ~= 0)
                %GO TO 30
                MYcase = 30;
                continue
            end
            
            % C  ------------------------------------------------------------------
            % C  ON THE FIRST CALL, THE ORDER IS SET TO 1 AND THE INITIAL YDOT
            % C  IS CALCULATED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE
            % C  INCREASED IN A SINGLE STEP.  RMAX IS SET EQUAL TO 1.D4 INITIALLY
            % C  TO COMPENSATE FOR THE SMALL INITIAL H, BUT THEN IS NORMALLY = 10.
            % C  IF A FAILURE OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST),
            % C  RMAX IS SET AT 2. FOR THE NEXT INCREASE.
            % C  ------------------------------------------------------------------

            [N,T,Y,SAVE1,IPAR,RPAR,IERR] = F(N,T,Y,SAVE1,IPAR,RPAR,IERR);

            if (IMAS~=0 )
                AM = reshape(AM,LDMAS,N);
                [N,AM,LDMAS,IPAR,RPAR,IERR] = MAS(N,AM,LDMAS,IPAR,RPAR,IERR);
                AM = reshape(AM,1,LDMAS*N);
            end

            if(IERR~=0)
                %GOTO 8000
                MYcase = 8000;
                continue
            end

            Y(1:N,2)=H*SAVE1(1:N);
            
            METH = 2;
            MITER = MF - 10*METH;
          
            IBND=5;
            UPBND=0.2;

            if(MF>22)
                UPBND=0.1;
            end

            NQ = 1;
            NQUSED = NQ;
            L = 2;
            IDOUB = 3;
            KFAIL = 0;
            RMAX = 10000;
            IER=0;
            RC = ZERO;
            CRATE1 = 0.1;
            CRATE2 = 0.1;
            JSNOLD = 0;
            %JNEWIM = .TRUE.
            JNEWIM = 1;
            TCRAT1 = ZERO;
            TCRAT2 = ZERO;
            VTOL=max([RTOL(1),ATOL(1)])/10;

            for I=1:12
                HSTPSZ(1,I)=1;
                HSTPSZ(2,I)=VTOL;
            end

            HOLD = H;
            MFOLD = MF;
            NSTEP = 0;
            NFE = 1;
            NJE = 0;
            NDEC = 0;
            NPSET = 0;
            NCOSET = 0;
            MAXORD = 1;
            NFAIL = 0;
            %CFAIL = .TRUE.
            CFAIL = 1;
            AVNEWJ = ZERO;
            AVOLDJ = ZERO;
            AVNEW2 = ZERO;
            AVOLD2 = ZERO;
            %SAMPLE = .FALSE.
            SAMPLE = 0;
            ISAMP = 0;
            IEMB=0;
            % C     **************************************************
            % C     CFAIL=.TRUE. ENSURES THAT WE CALCULATE A NEW
            % C     J ON THE FIRST CALL.
            % C     **************************************************
            MEQC1 = 0;
            MEQC2 = 0;
            MQ1TMP = 0;
            MQ2TMP = 0;
            NBSOL = 0;
            HUSED = H;
            % C  -----------------------------------------------------------------
            % C  IF THE CALLER HAS CHANGED N , THE CONSTANTS E, EDN, EUP
            % C  AND BND MUST BE RESET.  E IS A COMPARISON FOR ERRORS AT THE
            % C  CURRENT ORDER NQ.  EUP IS TO TEST FOR INCREASING THE ORDER,
            % C  EDN FOR DECREASING THE ORDER.  BND IS USED TO TEST FOR CONVERGENCE
            % C  OF THE CORRECTOR ITERATES.   IF THE CALLER HAS CHANGED H, Y MUST
            % C  BE RE-SCALED.  IF H IS CHANGED, IDOUB IS SET TO L+1 TO PREVENT
            % C  FURTHER CHANGES IN H FOR THAT MANY STEPS.
            % C  -----------------------------------------------------------------
            [NQ,EL,ELST,TQ,NCOSET,MAXORD] = COSET(NQ,EL,ELST,TQ,...
                NCOSET,MAXORD);
            LMAX = MAXDER + 1;
            RC = RC*EL(1)/OLDLO;
            OLDLO = EL(1);
            IWEVAL = MITER;
            NRENEW = 1;
            NEWPAR = 0;
            % C     *****************************************************
            % C     NRENEW AND NEWPAR ARE TO INSTRUCT ROUTINE THAT
            % C     WE WISH A NEW J TO BE CALCULATED FOR THIS STEP.
            % C     *****************************************************
            [N,TQ,EDN,E,EUP,BND,EDDN] = FERRORS(N,TQ,EDN,E,EUP,BND,EDDN);

            ARH(1:N)=EL(2)*Y(1:N,1);
            [Y,YHOLD] = CPYARY(N*L,Y,YHOLD);

            QI = H*EL(1);
            QQ = ONE/QI;
            
            

            [T,H,Y,L,N,SAVE2,NFE,IPAR,RPAR,IERR] = PRDICT(T, H, Y, L, N,...
                SAVE2,NFE, IPAR,RPAR,IERR);
            

            if (IERR~=0)
                H = H/2;
                IERR = 0;
                %goto 6000
                MYcase = 6000;
                continue
            end
            %GO TO 110
            MYcase = 110;
            continue

        case 30 
            %***********************************************************             
            % C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            % C     DIFFERENT PARAMETERS ON THIS CALL        <
            % C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            [YHOLD,Y] = CPYARY(N*L,YHOLD,Y);

            if (MF~=MFOLD)
                METH = MF/10;
                MITER = MF - 10*METH;
                MFOLD = MF;
                IWEVAL = MITER;
            end

            if (NSTEP > 0)
                %GOTO 35
            else

                NJE = 0;
                NFE = 1;
                %CFAIL = .TRUE.
                CFAIL = 1;
                NEWPAR = 0;
                MQ1TMP = 0;
                MQ2TMP = 0;
                MEQC1 = 0;
                MEQC2 = 0;
                TCRAT1 = 0;
                TCRAT2 = 0;
                CRATE1 = 1E-1;
                CRATE2 = 1E-1;
                NSTEP = 0;
                NBSOL = 0;
                NPSET = 0;
                NCOSET = 0;
                NDEC = 0;
                %35   CONTINUE
            end

            if (H~=HOLD)
                RH = H/HOLD;
                H = HOLD;
                IREDO = 3;
                %GO TO 50
                MYcase = 50;
                continue

            else
                %GO TO 60
                MYcase = 60;
                continue
            end

            error('should have gone to 50 or 60 or somewhere else')
        case 40%**********************************************************            
            % C     *********************************************
            % C     RE-SCALE Y AFTER A CHANGE OF STEPSIZE   *
            % C     *********************************************
            RH = max([RH,HMIN/abs(H)]);

            MYcase = 50;
            continue

        case 50%***********************************************************
            RH = min([RH,HMAX/abs(H),RMAX]);
            [N,L,RH,Y] = RSCALE(N,L,RH,Y);
            RMAX = 10;
            JCHANG = 1;
            H = H*RH;
            RC = RC*RH;
            if (JSNOLD>IBND)
                %CFAIL = .TRUE.
                CFAIL = 1;
                NEWPAR = 0;
                RC = ZERO;
                % C ***************************************************************
                % C   CFAIL=TRUE AND NEWPAR=0 SHOULD FORCE A NEW J TO BE EVALUATED
                % C   AFTER 7 STEPS WITH AN OLD J, IF WE HAVE HAD A FAILURE OF ANY
                % C   KIND ON THE FIRST, SECOND OR THIRD STAGE OF THE CURRENT STEP
                % C ***************************************************************
            end

            IDOUB = L + 1;
            [Y,YHOLD] = CPYARY(N*L,Y,YHOLD);

            MYcase = 60;
            continue

        case 60 %*********************************************************
              
            if (abs(RC-ONE) > UPBND)
                IWEVAL = MITER;
            end

            HUSED = H;
            % C  ------------------------------------------------------------------
            % C  THIS SECTION COMPUTES THE PREDICTED VALUES OF Y
            % C  AND THE RHS, ARH, FOR USE IN THE NEWTON ITERATION SCHEME.
            % C  RC IS THE RATIO OF THE NEW TO OLD VALUES OF THE COEFFICIENT
            % C  H*EL(1). WHEN RC DIFFERS FROM 1 BY MORE THAN 20 PERCENT, IWEVAL IS
            % C  SET TO MITER TO FORCE THE PARTIALS TO BE UPDATED.
            % C  ------------------------------------------------------------------
            QI = H*EL(1);
            QQ = ONE/QI;

            ARH(1:N)=EL(2)*Y(1:N,1);
            for J1 = 2:NQ
                JP1 =J1+1;
                for I = 1:N
                    ARH(I) = ARH(I) + EL(JP1)*Y(I,J1);
                end
            end

            if (JCHANG==1)
                % C    IF WE HAVE CHANGED STEPSIZE THEN PREDICT A VALUE FOR Y(T+H)
                % C    AND EVALUATE THE DERIVATIVE THERE (STORED IN SAVE2())
                [T,H,Y,L,N,SAVE2,NFE,IPAR,RPAR,IERR] = PRDICT(T, H, Y, ...
                    L, N, SAVE2, NFE,IPAR,RPAR,IERR);
                if (IERR~=0)
                    %GOTO 8000
                    MYcase = 8000;
                    continue
                end

            else
                % C  ELSE USE THE VALUES COMPUTED FOR THE SECOND BDF FROM THE LAST
                % C  STEP. Y( ,LMAX+2) HOLDS THE VALUE FOR THE DERIVATIVE AT (T+H)
                % C  AND Y( ,LMAX+3) HOLDS THE APPROXIMATION TO Y AT THIS POINT.
                LMP2=LMAX+2;
                LMP3=LMAX+3;

                for I = 1:N
                    SAVE2(I) = Y(I,LMP2);
                    Y(I,1) = Y(I,LMP3);
                end

                T = T + H;
            end

            MYcase = 110;
            continue

        case 110 %*********************************************************   
     
            if (IWEVAL <= 0)
                %GO TO 120
            else
                % C -------------------------------------------------------------------
                % C  IF INDICATED, THE MATRIX P = I/(H*EL(2)) - J IS RE-EVALUATED BEFORE
                % C  STARTING THE CORRECTOR ITERATION.  IWEVAL IS SET = 0 TO INDICATE
                % C  THAT THIS HAS BEEN DONE. P IS COMPUTED AND PROCESSED IN PSET.
                % C  THE PROCESSED MATRIX IS STORED IN PW
                % C -------------------------------------------------------------------
                IWEVAL = 0;
                RC = ONE;
                IITER = MEQC1 - MQ1TMP;
                IITER2 = MEQC2 - MQ2TMP;

                if (JNEWIM)
                    if (JSNOLD >= 3)
                        AVNEWJ = TCRAT1/IITER;
                        AVNEW2 = TCRAT2/IITER2;

                    else
                        AVNEWJ = ONE;
                        AVNEW2 = ONE;
                    end

                else
                    % C
                    % C          MATRIX P WAS FORMED WITH A COPY OF J
                    % C
                    if (JSNOLD >= 3)
                        AVOLDJ = TCRAT1/IITER;
                        AVOLD2 = TCRAT2/IITER2;
                        if (AVOLDJ < AVNEWJ)
                            AVNEWJ = AVOLDJ;

                        elseif ( ((abs(AVOLDJ-AVNEWJ)) > 0.3) || ...
                                ((AVOLDJ > 0.85) && (AVOLDJ~=ONE)) )
                            % C
                            % C     SINCE IN CERTAIN INSTANCES AVOLDJ WILL
                            % C     BE 1.0 AND THERE WILL BE NO NEED TO
                            % C     UPDATE J.
                            % C
                            %CFAIL = .TRUE.
                            CFAIL = 1;
                            CRATE1 = 0.1;
                            CRATE2 = 0.1;
                        end

                    else
                        %CFAIL = .TRUE.
                        CFAIL = 1;
                        CRATE1 = 0.1;
                        CRATE2 = 0.1;
                        % C
                        % C  **********************************************
                        % C  IF WE HAVE REACHED HERE THINGS MUST HAVE GONE WRONG
                        % C  **********************************************
                        % C
                    end

                end

                TCRAT1 = ZERO;
                TCRAT2 = ZERO;

                if (CFAIL)
                    NRENEW = 1;
                    NEWPAR = 1;
                    JSINUP = -1;
                    %JNEWIM = .TRUE.
                    JNEWIM = 1;
                end

                %CFAIL = .FALSE.
                CFAIL = 0;
                JSNOLD = 0;
                MQ1TMP = MEQC1;
                MQ2TMP = MEQC2;

                [Y, N, H, T, UROUND, EPSJAC, QI, MITER, MBND, MASBND, ...
                    NIND1, NIND2, NIND3, IER, NRENEW, YMAX, SAVE1, ...
                    SAVE2, PW, PWCOPY, AM, ERROR, IPIV, ITOL, RTOL, ...
                    ATOL, NPSET, NJE, NFE, NDEC, IPAR, RPAR, IERR] = ...
                    PSET(Y, N, H, T, UROUND, EPSJAC, QI, MITER, MBND, ...
                    MASBND, NIND1, NIND2, NIND3, IER, NRENEW, YMAX, ...
                    SAVE1, SAVE2, PW, PWCOPY, AM, ERROR, IPIV, ITOL, ...
                    RTOL, ATOL, NPSET, NJE, NFE, NDEC, IPAR,RPAR,IERR);

                if (IERR~=0)
                    %GOTO 8000
                    MYcase = 8000;
                    continue
                end

                QQQ=QI;

                % C     NOTE THAT ERROR() IS JUST BEING USED AS A WORKSPACE BY PSET
                if (IER~=0)
                    % C     IF IER>0 THEN WE HAVE HAD A SINGULARITY IN THE ITERATION MATRIX
                    IJUS=1;
                    RED=0.5;
                    NFAIL = NFAIL + 1;
                    %GO TO 450
                    MYcase = 450;
                    continue

                end


            end        %120

            for I = 1:N
                SAVE1(I) = Y(I,1);
                ERROR(I) = ZERO;
            end
            M1 = 0;
            % C **********************************************************************
            % C     UP TO 4 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS MADE
            % C     ON THE R.M.S. NORM OF EACH CORRECTION ,USING BND, WHICH DEPENDS
            % C     ON ATOL AND RTOL.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE
            % C     VECTOR  ERROR(I).  THE Y ARRAY IS NOT ALTERED IN THE CORRECTOR
            % C     LOOP. THE UPDATED Y VECTOR IS STORED TEMPORARILY IN SAVE1.
            % C **********************************************************************
            %       IF (.NOT.SAMPLE) THEN
            if (SAMPLE == 0)                                            
               
                [QQQ, Y, N, T, QI, BND, ARH, CRATE1, TCRAT1, M1, WORKED,...
                    YMAX, ERROR, SAVE1, SAVE2, SCALE, PW, MF, MBND, AM, ...
                    MASBND, NIND1, NIND2, NIND3, IPIV, ITOL, RTOL, ATOL,...
                    IPAR, RPAR, HUSED, NBSOL, NFE, NQUSED, IERR] ...
                    = ITRAT2(QQQ, Y, N, T, QI, BND, ARH, CRATE1, TCRAT1,...
                    M1, WORKED, YMAX, ERROR, SAVE1, SAVE2, SCALE, PW, ...
                    MF, MBND, AM, MASBND, NIND1, NIND2, NIND3, IPIV, 1,...
                    ITOL, RTOL, ATOL, IPAR, RPAR, HUSED, NBSOL, NFE, ...
                    NQUSED, IERR);
                
                    
                if (IERR ~= 0)
                    %GOTO 8000
                    MYcase = 8000;
                    continue
                end

                ITST = 2;

            else   

                [QQQ, Y, N, T, QI, BND, ARH, CRATE1, TCRAT1, M1, WORKED,...
                    YMAX, ERROR, SAVE1, SAVE2, SCALE, PW, MF, MBND, AM, ...
                    MASBND, NIND1, NIND2, NIND3, IPIV, ITOL, RTOL, ATOL,...
                    IPAR, RPAR, HUSED, NBSOL, NFE, NQUSED, IERR] ...
                    = ITRAT2(QQQ, Y, N, T, QI, BND, ARH, CRATE1, TCRAT1,...
                    M1, WORKED, YMAX, ERROR, SAVE1, SAVE2, SCALE, PW, ...
                    MF, MBND, AM, MASBND, NIND1, NIND2, NIND3, IPIV, 0,...
                    ITOL, RTOL, ATOL, IPAR, RPAR, HUSED, NBSOL, NFE,...
                    NQUSED, IERR);

                if (IERR ~= 0)
                    %GOTO 8000
                    MYcase = 8000;
                    continue
                end

                ITST = 3;
            end

            MEQC1 = MEQC1 + M1 + 1;

            % C
            % C       NOW TEST TO SEE IF IT WAS SUCCESSFUL OR NOT
            % C
            % C

            %IF (.NOT.WORKED) THEN            
            if (WORKED==0)

                NFAIL = NFAIL + 1;
                % C **********************************************************************
                % C        THE CORRECTOR ITERATION FAILED TO CONVERGE IN 4 TRIES. IF
                % C        PARTIALS ARE NOT UP TO DATE, THEY ARE RE-EVALUATED FOR THE
                % C        NEXT TRY. OTHERWISE THE Y ARRAY IS REPLACED BY ITS VALUES
                % C        BEFORE PREDICTION AND H IS REDUCED IF POSSIBLE. IF NOT A
                % C        NON-CONVERGENCE EXIT IS TAKEN
                % C **********************************************************************
                if (IWEVAL==-1)
                    % C           HAVE BEEN USING OLD PARTIALS, UPDATE THEM AND TRY AGAIN
                    IWEVAL = MITER;
                    %CFAIL = .TRUE.
                    CFAIL = 1;
                    [N,T,Y,SAVE2,IPAR,RPAR,IERR] = F(N,T,Y,SAVE2,IPAR,...
                        RPAR,IERR);
                    if (IERR ~= 0)
                        %GOTO 8000
                        MYcase = 8000;
                        continue
                    end

                    NFE = NFE + 1;
                    %GO TO 110
                    MYcase = 110;
                    continue

                end

                IJUS=0;
                RED=0.5;
                % C    ***    failed at step 1 because of Newton

                %GO TO 450
                MYcase = 450;
                continue

            end

            IWEVAL = -1;
            HUSED = H;
            NQUSED = NQ;           
             
            Y(1:N,1)=(SAVE1(1:N)-ARH(1:N));

            if(IMAS==0)
               
                for I=1:N
                    SAVE2(I) = Y(I,1)*QQ;
                    Y(I,1) = SAVE1(I);
                end
               
            else
                if((MF<= 22)&&(MLMAS==N))
                    for I=1:N
                        SAVE2(I) = 0;
                        for J=1:N
                            SAVE2(I) = SAVE2(I) + AM((J-1)*N+I)*Y(J,1);
                        end
                    end
                else
                    C1 = MLMAS + MUMAS;
                    C2 = C1 + 1;

                    for I=1:N
                        SAVE2(I) = 0;
                        for J=-MLMAS:MUMAS
                            ISUM = I+J;

                            if ((ISUM<1)||(ISUM>N))
                            else

                                KK = C2*I + C1*J-MLMAS;
                                SAVE2(I) = SAVE2(I) + AM(KK)*Y(I+J,1);
                                
                            end
                            
                        end
                    end
                end
                
                for I=1:N
                    Y(I,1)=SAVE1(I);
                    SAVE2(I) = SAVE2(I)*QQ;
                end
            end
            % C
            % C     UPDATE THE DIFFERENCES AT N+1
            % C
            for J = 2:L
                JM1 = J-1;
                for I = 1:N
                    Y(I,J) = Y(I,JM1) - YHOLD(I,JM1);
                end
            end
            
            % C
            % C     COMPUTE ERROR IN THE SOLUTION
            % C
            for I=1:N
                AYI = abs(Y(I,1));
                if(ITOL==1)
                    SCALE(I) = YMAX(I);
                elseif (ITOL==2)
                    SCALE(I) = RTOL(1)*AYI + ATOL(1);
                elseif (ITOL==3)
                    SCALE(I) = RTOL(1)*AYI + ATOL(I);
                elseif (ITOL==4)
                    SCALE(I) = RTOL(I)*AYI + ATOL(1);
                elseif (ITOL==5)
                    SCALE(I) = RTOL(I)*AYI + ATOL(I);
                end
            end

            if (NIND2~=0)
                for I = (NIND1+1):(NIND2+NIND1)
                    SCALE(I)=SCALE(I)/HUSED;
                end
            end

            if (NIND3~=0)
                for I = (NIND1 +NIND2 + 1):(NIND1+NIND2+NIND3)
                    SCALE(I)=SCALE(I)/(HUSED^2);
                end
            end
            % C ****
            % C ****  AMMEND
            % C ****  CHANGE 1,N BELOW TO 1,NVARS
            % C ****
            D = ZERO;
            for I = 1:N
                D = D + ((Y(I,L)-YHOLD(I,L))/SCALE(I))^2;
            end
            % c
            % c    STORING Y FROM FIRST STEP FOR USE IN THIRD STEP.
            % C
            if (ITOL==1)
                D = D/(RTOL(1)^2);
            end

            if (D>E)
                %GOTO 330
                MYcase = 330;
                continue
            end
	    YNHOLD(1:N,1)=Y(1:N);
	    YNHOLD(1:N,2)=SAVE2(1:N);

            KFAIL = 0;
            IREDO = 0;
            
	    ARH(1:N)=EL(2)*Y(1:N,1);

            for J1 = 2:NQ
                JP1 = J1+1;
                for I = 1:N
                    ARH(I) = ARH(I) + EL(JP1)*Y(I,J1);
                end
            end

            [T,H,Y,L,N,SAVE2,NFE,IPAR,RPAR,IERR] = PRDICT(T, H, Y, L, N,...
                SAVE2, NFE, IPAR, RPAR, IERR);

            if(IERR~=0)
                %GOTO 8000
                MYcase = 8000;
                continue
            end

            for I = 1:N
                SAVE1(I) = Y(I,1);
                ERROR(I) = ZERO;
            end

            M2 = 0;
            % C
            % C     FOR NOW WILL ASSUME THAT WE DO NOT WISH TO SAMPLE
            % C     AT THE N+2 STEP POINT
            % C

 
            [QQQ, Y, N, T, QI, BND, ARH, CRATE2, TCRAT2, M2, WORKED, ...
                YMAX, ERROR, SAVE1, SAVE2, SCALE, PW, MF, MBND, AM, ...
                MASBND, NIND1, NIND2, NIND3, IPIV, ITOL, RTOL, ATOL, ...
                IPAR, RPAR, HUSED, NBSOL, NFE, NQUSED, IERR]=ITRAT2(QQQ,...
                Y, N, T, QI, BND, ARH, CRATE2, TCRAT2, M2, WORKED, YMAX,...
                ERROR, SAVE1, SAVE2, SCALE, PW, MF, MBND, AM, MASBND, ...
                NIND1, NIND2, NIND3, IPIV, 1, ITOL, RTOL, ATOL, IPAR, ...
                RPAR, HUSED, NBSOL, NFE, NQUSED, IERR);

            if (IERR~=0)
                %GOTO 8000
                MYcase = 8000;
                continue
            end

            MEQC2 = MEQC2 + M2 + 1;
            % C
            % C       NOW CHECK TO SEE IF IT WAS SUCCESSFUL OR NOT
            % C
            %IF (.NOT.WORKED) THEN
            if (WORKED==0)
                NFAIL = NFAIL + 1;
                IJUS=0;
                RED=0.5;
                % C    ***have failed on step 2
                %GOTO 450
                MYcase = 450;
                continue

            end
            % C
            % C    IF WE ARE DOWN TO HERE THEN THINGS MUST HAVE CONVERGED
            % C
            LMP2=LMAX+2;
            LMP3=LMAX+3;

            Y(1:N,LMP3)=(SAVE1(1:N)-ARH(1:N));

            if (IMAS==0)
                for I=1:N
                    Y(I,LMP2) = Y(I,LMP3)*QQ;
                    Y(I,LMP3) = SAVE1(I);
                end
            else
                if ((MF<=22)&&(MLMAS==N))
                    for I=1:N
                        Y(I,LMP2) = 0;
                        for J=1:N
                            Y(I,LMP2) = Y(I,LMP2) + ...
                                AM((J-1)*N+I)*Y(J,LMP3);
                        end
                    end

                else
                    C1 = MLMAS + MUMAS;
                    C2 = C1 + 1;
                    for I=1:N
                        Y(I,LMP2)=0;
                        for J=-MLMAS:MUMAS
                            ISUM = I+J;
                            if ((ISUM<1)||(ISUM>N))

                            else
                                KK = C2*I + J*C1 -MLMAS;
                                Y(I,LMP2)=Y(I,LMP2)+AM(KK)*Y(I+J,LMP3);
                            end
                        end
                    end

                end

                for I=1:N
                    Y(I,LMP2) = Y(I,LMP2)*QQ;
                    Y(I,LMP3) = SAVE1(I);
                end
            end

            % C
            % C     WE ARE NOW COMPUTING THE THIRD STAGE
            % C
            LL = L + 1;
            T = TOLD + H;
            DELST = ELST(1)-EL(1);
            NQP2 = NQ+2;
            LMP2=LMAX+2;
            for I=1:N
                ARH(I) = H*(ELST(NQP2)*Y(I,LMP2)+DELST*YNHOLD(I,2));
                for J1 = 1:NQ
                    ARH(I) = ARH(I) + ELST(J1+1)*YHOLD(I,J1);
                end
            end

            for I = 1:N
                SAVE2(I) = YNHOLD(I,2);
                Y(I,1) = YNHOLD(I,1);
            end

            M3STEP = 0;

            MYcase = 300;
            continue

        case 300 
            %********************************************************

            if(IMAS==0)
                for I=1:N
                    SAVE1(I) = (-Y(I,1)  + QI*SAVE2(I) + ARH(I)) /QQQ;
                end
            else
                for I=1:N
                    SAVE1(I) = -Y(I,1)+ARH(I)-H*(ELST(NQP2)*Y(I,LMP2)+...
                        DELST*YNHOLD(I,2));
                end

                if ((MF<=22)&&(MLMAS==N))
                    for I=1:N
                        SAVE2(I) = (SAVE2(I)*QI + ...
                            H*(ELST(NQP2)*Y(I,LMP2) + ...
                            DELST*YNHOLD(I,2)));
                        for J=1:N
                            SAVE2(I) = SAVE2(I) + AM((J-1)*N+I)*SAVE1(J);
                        end
                    end
                else
                    C1 = MLMAS + MUMAS;
                    C2 = C1 + 1;
                    for I=1:N
                        SAVE2(I) = (SAVE2(I)*QI + H*(ELST(NQP2)*Y(I,LMP2)+...
                            DELST*YNHOLD(I,2)));
                        for J=-MLMAS:MUMAS
                            ISUM = I+J;
                            if ((ISUM<1)|| (ISUM>N))

                            else
                                KK=C2*I+C1*J-MLMAS;
                                SAVE2(I) = SAVE2(I) + AM(KK)*SAVE1(I+J);
                            end
                        end
                    end
                end

                SAVE1(1:N)=SAVE2(1:N)/QQQ;
            end

            if (MF>=23)
                PW = reshape(PW,MBND(4),N);
                [PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE1] = DGBSL(PW, ...
                    MBND(4), N, MBND(1), MBND(2), IPIV, SAVE1, 0);
                PW = reshape(PW,1,MBND(4)*N);
                NBSOL=NBSOL+1;
            else
            
                PW = reshape(PW,N,N);

                [N,N,PW,SAVE1,IPIV] = SOL(N,N,PW,SAVE1,IPIV);

                PW = reshape(PW,1,N*N);
    
                NBSOL = NBSOL + 1;

            end

            for I=1:N
                AYI = abs(Y(I,1));
                if(ITOL==1)
                    SCALE(I) = YMAX(I);
                elseif (ITOL==2)
                    SCALE(I) = RTOL(1)*AYI + ATOL(1);
                elseif (ITOL==3)
                    SCALE(I) = RTOL(1)*AYI + ATOL(I);
                elseif (ITOL==4)
                    SCALE(I) = RTOL(I)*AYI + ATOL(1);
                elseif (ITOL==5)
                    SCALE(I) = RTOL(I)*AYI + ATOL(I);
                end
            end

            if(NIND2~=0)

                for I = (NIND1+1):(NIND2+NIND1)
                    SCALE(I)=SCALE(I)/HUSED;
                end

            end

            if (NIND3~=0)

                for I = (NIND1+NIND2 + 1):(NIND1+NIND2+NIND3)
                    SCALE(I) = SCALE(I)/(HUSED^2);
                end

            end

            D = ZERO;

            for I = 1:N
                D = D + (SAVE1(I)/SCALE(I))^2;
                Y(I,1) = Y(I,1) + SAVE1(I);
            end

            if (ITOL==1)
                D = D/(RTOL(1)^2);
            end

            if ((D*min([ONE,2*CRATE1]))<=BND)
                %GO TO 360
                MYcase = 360;
                continue
            end

            if (M3STEP==5)
                %c         WRITE (LOUT,9000)
                IJUS=1;
                RED=0.5;
                %C    ****  step 3 fails
                NFAIL = NFAIL + 1;
                %GO TO 450
                MYcase = 450;
                continue
            end

            M3STEP = M3STEP + 1;
            [N,T,Y,SAVE2,IPAR,RPAR,IERR] = F(N,T,Y,SAVE2,IPAR,RPAR,IERR);
            if (IERR~=0)
                %GOTO 8000
                MYcase = 8000;
                continue
            end
            NFE = NFE + 1;
            %GO TO 300
            MYcase = 300;
            continue

        case 330 %********************************************************

            KFAIL = KFAIL - 1;
            % C **********************************************************************
            % C     THE ERROR TEST FAILED. KFAIL KEEPS TRACK OF MULTIPLE FAILURES.
            % C     RESTORE T AND THE Y ARRAY TO THEIR PREVIOUS VALUES AND PREPARE TO
            % C     TRY THE STEP AGAIN. COMPUTE THE OPTIMAL STEP SIZE FOR THIS ORDER
            % C     AND ONE ORDER LOWER.
            % C **********************************************************************
            % C     ***  failed on step 1 because of accuracy
            % C     COMPUTE ERROR IN THE SOLUTION
            % C
            NFAIL = NFAIL + 1;

            for I=1:N
                AYI = abs(Y(I,1));
                if(ITOL==1)
                    SCALE(I) = YMAX(I);
                elseif (ITOL==2)
                    SCALE(I) = RTOL(1)*AYI + ATOL(1);
                elseif (ITOL==3)
                    SCALE(I) = RTOL(1)*AYI + ATOL(I);
                elseif (ITOL==4);
                    SCALE(I) = RTOL(I)*AYI + ATOL(1);
                elseif (ITOL==5)
                    SCALE(I) = RTOL(I)*AYI + ATOL(I);
                end
            end

            if (NIND2~=0)
                for I = (NIND1+1):(NIND2+NIND1)
                    SCALE(I)=SCALE(I)/HUSED;
                end
            end

            if (NIND3~=0)
                for I = (NIND1+NIND2 + 1):(NIND1+NIND2+NIND3)
                    SCALE(I)=SCALE(I)/(HUSED^2);
                end
            end

            DDOWN = ZERO;
            TWODWN = ZERO;

            for I=1:N
                DDOWN = DDOWN + ((Y(I,L))/SCALE(I))^2;
                TWODWN = TWODWN + ((Y(I,L-1))/SCALE(I))^2;
            end

            if (ITOL==1)
                D = D/(RTOL(1)^2);
            end

            T = TOLD;
            HOLD = H;


            if (NQ > 1)
                FFAIL = 0.5/NQ;
            end

            if (NQ>2)
                FRFAIL = 0.5/(NQ-1);
            end

            EFAIL = 0.5/L;
            [YHOLD,Y] = CPYARY(N*L,YHOLD,Y);
            RMAX = 2;

            if (abs(H)<=HMIN*1.00001)
                % C
                % C        REQUESTED ERROR NOT POSSIBLE WITH GIVEN HMIN
                % C
                KFLAG = -1;
                HOLD = H;
                return

            end

            if (KFAIL<=-3)
                %GO TO 340
                MYcase = 340;
                continue
            end

            IREDO = 2;
            % C
            % C     PREDICTING A NEW H AFTER INSUFFICIENT ACCURACY
            % C
            PRFAIL = ((D/(0.2*E))^EFAIL)*1.5 + 1.6E-6;
            PLFAIL = ((DDOWN/(0.2*EDN))^FFAIL)*1.5+1.7E-6;

            if (NQ > 2)
                PLLFAL =((TWODWN/(0.2*EDDN))^FRFAIL)*1.5+1.7E-6;
            end

            if (PLLFAL>PLFAIL)
                PLFAIL=PLLFAL;
            end

            if ((PLFAIL < PRFAIL)&&(NQ~=1))

                NEWQ=NQ-1;
                NQ=NEWQ;
                RH=ONE/(PLFAIL*(-KFAIL));
                L=NQ+1;
                [NQ,EL,ELST,TQ,NCOSET,MAXORD] = COSET(NQ,EL,ELST,TQ,...
                    NCOSET,MAXORD);
                EL(1);
                RC=RC*EL(1)/OLDLO;
                OLDLO=EL(1);
                [N,TQ,EDN,E,EUP,BND,EDDN] = FERRORS(N,TQ,EDN,E,EUP,...
                    BND,EDDN);
            else
                NEWQ = NQ;
                RH = ONE/ (PRFAIL*(-KFAIL));
            end
            %GO TO 40
            MYcase = 40;
            continue


        case 340 %*********************************************************

            % C **********************************************************************
            % C     CONTROL REACHES THIS STAGE IF 3 OR MORE FAILURES HAVE OCCURED.
            % C     IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE Y
            % C     ARRAY HAVE ERRORS OF THE WRONG ORDER. HENCE THE FIRST DERIVATIVE
            % C     IS RE-COMPUTED, AND THE ORDER IS SET TO 1. THEN H IS REDUCED BY A
            % C     FACTOR OF 10, AND THE STEP IS RETRIED. AFTER A TOTAL OF 7
            % C     FAILURES AN EXIT IS TAKEN WITH KFLAG=-2.
            % C **********************************************************************


            if (KFAIL==-7)
                % C        ERROR SMALLER THAN CAN BE HANDLED FOR PROBLEM
                KFLAG = -2;
                HOLD = H;
                return

            end
            % C     *********************************
            % C     START FROM ORDER 1 AGAIN    *
            % C     *********************************

            JCHANG = 1;
            RH = max([HMIN/abs(H),0.1]);
            [RH,H,OVRIDE] = HCHOSE(RH,H,OVRIDE);
            H = H*RH;
            [N,T,YHOLD,SAVE1,IPAR,RPAR,IERR] = F(N,T,YHOLD,SAVE1,IPAR,RPAR,IERR);
            if (IERR~=0)
                %GOTO 8000
                MYcase  = 8000;
                continue
            end

            NFE = NFE + 1;

            for I = 1:N
                Y(I,1) = YHOLD(I,1);
                Y(I,2) = H*SAVE1(I);
                YHOLD(I,2) = Y(I,2);
            end

            IWEVAL = MITER;
            %CFAIL = .TRUE.
            CFAIL = 1;

            % C     SINCE WE HAVE HAD PROBLEMS PROCEED WITH THIS ORDER
            % C     FOR 10 STEPS (IF WE CAN)

            IDOUB = 10;
            if (NQ==1)
                %GO TO 60
                MYcase = 60;
                continue
            end

            NQ = 1;
            L = 2;

            % C     RESET ORDER, RECALCULATE ERROR BOUNDS

            [NQ,EL,ELST,TQ,NCOSET,MAXORD] = COSET(NQ,EL,ELST,TQ,...
                NCOSET,MAXORD);
            LMAX = MAXDER + 1;
            RC = RC*EL(1)/OLDLO;
            OLDLO = EL(1);
            [N,TQ,EDN,E,EUP,BND,EDDN] = FERRORS(N,TQ,EDN,E,EUP,BND,EDDN);

            % C     NOW JUMP TO NORMAL CONTINUATION POINT

            %GO TO 60
            MYcase = 60;
            continue



            % C **********************************************************************
            % C     THE ITERATION FOR THE CORRECTED SOLUTION HAS CONVERGED.
            % C     UPDATE THE Y ARRAY.
            % C **********************************************************************


        case 360 %*********************************************************

            % C   ****
            % C   **** AMMEND ****
            % C   **** CHANGE 1,N BELOW TO 1,NVARS
            % C   ****
            DEMB=0;

            for I=1:N
                DEMB=DEMB+((Y(I,1)-YNHOLD(I,1))/SCALE(I))^2;
            end


            if (DEMB>4*N)
                IEMB=1;
                IJUS=1;
                RED=0.5;
                % C     ***  failed because of embedded error estimate
                NFAIL = NFAIL + 1;
                %GOTO 450
                MYcase = 450;
                continue
            end

            for J2 = 2:LL
                J2M1=J2-1;
                for I = 1:N
                    Y(I,J2) = Y(I,J2M1) - YHOLD(I,J2M1);
                end
            end
            % C ---------------------------------------------------------------------
            % C     IF THE COUNTER IDOUB EQUALS 2 AND WE ARE NOT ALREADY USING THE
            % C     MAXIMUM ALLOWABLE ORDER , STORE Y(I,LMAX+4) WHICH IS USED IN
            % C     ASSESSING THE POSSIBILITY OF INCREASING THE ORDER. IF IDOUB = 0
            % C     CONTROL PASSES TO 480 WHERE AN ATTEMPT TO CHANGE THE STEPSIZE AND
            % C     ORDER IS MADE.
            % C ----------------------------------------------------------------------

            if ((IDOUB==2)&&(L~=LMAX))
                LMP4=LMAX+4;
                Y(1:N,LMP4)=Y(1:N,LL);
            end

            IDOUB = IDOUB - 1;

            TRANGE=(TEND-TOLD-H)*H;

            if (TRANGE< 0)
                IDOUB = IDOUB + 2;
                %GOTO 440
                MYcase = 440;
                continue
            end

            JCHANG = 0;
            if (IDOUB==0)
                %SAMPLE = .FALSE.
                SAMPLE = 0;

                ISAMP = ISAMP + 1;

                if (ISAMP==4)
                    %SAMPLE = .TRUE.
                    SAMPLE = 1;
                    ISAMP = 0;
                end


                % C **********************************************************************
                % C        NOW COMPUTE THE FACTORS PR1, PR2 AND PR3, BY WHICH
                % C        H COULD BE DIVIDED AT ORDER NQ-1, ORDER NQ AND ORDER NQ+1
                % C        RESPECTIVELY. THE SMALLEST OF THESE IS DETERMINED AND THE NEW
                % C        ORDER CHOSEN ACCORDINGLY. IF THE ORDER IS TO BE INCREASED WE
                % C        MUST COMPUTE ONE MORE BACKWARD DIFFERENCE.
                % C **********************************************************************

                PR3 = 1.E+20;
                FAC = 1.5;

                if (IEMB==1)
                    FAC = 1.8;
                end

                for I = 1:N
                    AYI = abs(Y(I,1));
                    if(ITOL==1)
                        VHOLD = YMAX(I);
                    elseif (ITOL==2)
                        VHOLD = RTOL(1)*AYI + ATOL(1);
                    elseif (ITOL==3)
                        VHOLD = RTOL(1)*AYI + ATOL(I);
                    elseif (ITOL==4)
                        VHOLD = RTOL(I)*AYI + ATOL(1);
                    elseif (ITOL==5)
                        VHOLD = RTOL(I)*AYI + ATOL(I);
                    end
                    SCALE(I)=VHOLD;
                end

                if (NIND2~=0)
                    for I=(NIND1+1):(NIND1+NIND2)
                        SCALE(I) = SCALE(I)/HUSED;
                    end
                end

                if (NIND3~=0)
                    for I=(NIND1+NIND2+1):N
                        SCALE(I) = SCALE(I)/(HUSED^2);
                    end
                end


                if (L~=LMAX)
                    LMP4 = LMAX + 4;
                    DUP = ZERO;
                    for I=1:N
                        DUP = DUP + ((Y(I,LL)-Y(I,LMP4))/SCALE(I))^2;
                    end

                    if (ITOL==1)
                        DUP = DUP/(RTOL(1)^2);
                    end

                    ENQ3 = 0.5/(L+1);
                    PR3 = ((DUP/EUP)^ENQ3)*(FAC+0.2) + 1.8E-6;

                end


                ENQ2 = 0.5/L;
                D = ZERO;
                DDOWN=ZERO;

                for I = 1:N
                    D = D + (Y(I,LL)/SCALE(I))^2;
                    % c           DDOWN = DDOWN + (Y(I,LMP4)/SCALE(I))**2
                end

                if (ITOL==1)
                    D = D/(RTOL(1)^2);
                end

                PR2 = ((D/E)^ENQ2)*FAC + 1.6E-6;

                PR1 = 1E20;

                if (NQ>1)
                    DDDOWN=ZERO;
                    DDOWN = ZERO;
                    for I = 1:N
                        DDOWN = DDOWN + (Y(I,L)/SCALE(I))^2;
                        DDDOWN = DDDOWN + (Y(I,L-1)/SCALE(I))^2;
                    end

                    if (ITOL==1)
                        DDOWN = DDOWN/(RTOL(1)^2);
                    end

                    ENQ1 = 0.5/NQ;

                    PR1 = ((DDOWN/EDN)^ENQ1)*(FAC+0.1) + 1.7E-6;

                    if (NQ>2)
                        ENQ0 = 0.5/(NQ-1);

                        PR0 = ((DDDOWN/EDDN)^ENQ0)*(FAC+0.1) + 1.7E-6;

                        if (PR0>PR1)
                            PR1 = PR0;
                        end

                        if (DDDOWN<DDOWN)
                            DDOWN = DDDOWN;
                        end

                    end
                end

                if (L==LMAX)
                    DUP = 0;
                end

                if (NQ<=1)
                else

                    if ((DUP>D)&&(D>DDOWN))
                        PR2=1E+30;
                        PR3=1E+30;
                    end

                end

                if (PR2<=PR3)
                    if (PR2>PR1)
                        NEWQ = NQ - 1;
                        RH = 1/PR1;

                    else
                        NEWQ = NQ;
                        RH = 1/PR2;
                    end

                elseif (PR3<PR1)
                    NEWQ = L;
                    RH = 1/PR3;

                else
                    NEWQ = NQ - 1;
                    RH = 1/PR1;
                end

                IEMB=0;

                if ((RH>1)&&(RH<1.1))
                    IDOUB=10;
                    NQ=NQUSED;
                    L=NQ+1;
                    %GOTO 440
                    MYcase = 440;
                    continue
                end

                RH = min([RH,RMAX]);

                [RH,H,OVRIDE] = HCHOSE(RH,H,OVRIDE);

                if ((JSINUP<=20)&&(KFLAG==0)&&(RH<1.1))
                    % C           WE HAVE RUN INTO PROBLEMS
                    IDOUB = 10;
                    NQ = NQUSED;
                    L = NQ + 1;
                    %GO TO 440
                    MYcase = 440;
                    continue

                end
                % C **********************************************************************
                % C        IF THERE IS A CHANGE IN ORDER, RESET NQ, L AND THE
                % C        COEFFICIENTS. IN ANY CASE H IS RESET  AND THE
                % C        Y ARRAY IS RE-SCALED
                % C **********************************************************************

                if((IMAS~=0)&&(NIND3~=0))
                    if((NQ<=2)&&(PR3<1))
                        NEWQ=NQ+1;
                    end
                end

                if (NEWQ~=NQ)
                    if (NEWQ>NQ)
                        %C              ADD AN EXTRA TERM TO THE HISTORY ARRAY
                        Y(1:N,LL)=Y(1:N,L)-YHOLD(1:N,L);
                    end


                    NQ = NEWQ;
                    L = NQ + 1;

                    % C           RESET ORDER,RECALCULATE ERROR BOUNDS

                    [NQ,EL,ELST,TQ,NCOSET,MAXORD] = COSET(NQ,EL,ELST,...
                        TQ,NCOSET,MAXORD);
                    LMAX = MAXDER + 1;
                    RC = RC*EL(1)/OLDLO;
                    OLDLO = EL(1);
                    [N,TQ,EDN,E,EUP,BND,EDDN] = FERRORS(N,TQ,EDN,E,EUP,...
                        BND,EDDN);
                end


                RH = max([RH,HMIN/abs(H)]);
                RH = min([RH,HMAX/abs(H),RMAX]);
                [N,L,RH,Y] = RSCALE(N,L,RH,Y);
                RMAX = 10;
                JCHANG = 1;
                H = H*RH;
                RC = RC*RH;
                if (JSNOLD>IBND)
                    RC=ZERO;
                end


                IDOUB = L + 1;

            end

            MYcase = 440;
            continue

        case 440 
            %*******************************************************
            % C ----------------------------------------------------------------------
            % C     STORE THE Y ARRAY IN THE MATRIX YHOLD.  STORE IN THE Y ARRAY THE
            % C     INFORMATION NECESSARY TO PERFORM AN INTERPOLATION TO FIND THE
            % C     SOLUTION AT THE SPECIFIED OUTPUT POINT IF APPROPRIATE.
            % C ----------------------------------------------------------------------
            [Y,YHOLD] = CPYARY(N*L,Y,YHOLD);

            NSTEP = NSTEP + 1;
            JSINUP = JSINUP + 1;
            JSNOLD = JSNOLD + 1;
            JSTART = NQUSED;
            T = TOLD + HUSED;
            HOLD = H;
            KFAIL = 0;
            NEWPAR = 0;
            %       CFAIL = .FALSE.
            CFAIL =0;
            return

        case 450
            %*********************************************************
            %FINISH = .FALSE.
            FINISH = 0;

            T=TOLD;
            RMAX=2;
            for J1=1:L
                for I=1:N
                    Y(I,J1)=YHOLD(I,J1);
                end
            end

            if (abs(H)<=HMIN*1.00001)
                % C
                % C     CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED
                % C
                if (NSTEP==0)
                    KFLAG=-1;
                else
                    KFLAG=-3;
                end
                % C
                % C    TO SUPPRESS ERROR MESSAGES AT START AS H MAY
                % C    HAVE BEEN TOO LARGE ON THE FIRST STEP.
                % C
                HOLD=H;
                %FINISH = .TRUE.
                FINISH =1;

            end
            RH = RED;
            IREDO=1;
            % C
            % C     TRY AGAIN WITH UPDATED PARTIALS
            % C

            MYcase = 8000;
            continue

        case 8000

            if (IERR~=0)
                disp(['IERR IS NON-ZERO BECAUSE OF AN ILLEGAL' ...
                    ' FUNCTION CALL']);
                %write(LOUT,1975)             
                h= h/2;
                if (H<EPSJAC/100)
                    fdisp() %WRITE(6,9161)

                    %9161       FORMAT (/,/,'STEPSIZE IS TOO SMALL')
                    IDID = -7;
                    return
                end

                T= TOLD;

                if ((T-TOUT)*H >= 0)
                    % C           HAVE OVERSHOT TOUT
                    disp('%f \n %f \n %f \n',T,TOUT,H);
                    %WRITE (LOUT,*) T,TOUT,H
                    [N,JSTART,H,T,Y,TOUT,Y0] = INTERP(N,JSTART,H,T,...
                        Y,TOUT,Y0);
                    HO = H;
                    T0 = TOUT;
                    IDID = -5;
                    return
                end

                IERR = 0;
                JSTART = -1;
                %goto 30
                MYcase = 30;
                continue
            end

            if (IJUS==0)
                [RH,H,OVRIDE] = HCHOSE(RH,H,OVRIDE);
            end

            %IF(.NOT.FINISH) THEN
            if (FINISH==0)
                %GO TO 40
                MYcase = 40;
                continue
            else

                return
            end
    end
end
%========================================================================
%========================================================================

function [N,L,RH,Y] = RSCALE(N,L,RH,Y)
% 
% C     SUBROUTINE IS FOR RESCALING THE HISTORY ARRAY AFTER A CHANGE IN
% C     STEPSIZE
% C
% C     N      ORDER OF THE PROBLEM
% C     L      NUMBER OF TERMS IN THE HISTORY ARRAY TO BE RESCALED
% C     RH     RATIO OF THE STEPSIZE CHANGE (I.E. RH = HNEW/HOLD)
% C     Y()    THE HISTORY ARRAY
% C
% 
ZERO = 0;
% C     ..

DI(2,2) = RH;

if (L > 2)
    TA = RH*RH;
    DI(2,3) = RH* (1-RH)/2;
    DI(3,3) = TA;
    if (L > 3)
        TB = TA*RH;
        DI(2,4) = RH* ((RH-3)*RH+2)/6;
        DI(3,4) = TA* (1-RH);
        DI(4,4) = TB;
        if (L > 4)
            TC = TB*RH;
            DI(2,5) = - (((RH-6)*RH+11)*RH-6)*RH/24;
            DI(3,5) = TA* ((7*RH-18)*RH+11)/12;
            DI(4,5) = 1.5*TB*(1-RH);
            DI(5,5) = TC;
            if (L>5)
                TD = TC*RH;
                DI(2,6) = ((((RH-10)*RH+35)*RH-50)*RH+24.0D+0)*RH/120;
                DI(3,6) = - (((3*RH-14)*RH+21)*RH-10)*TA/12;
                DI(4,6) = ((5*RH-12)*RH+7)*TB/4;
                DI(5,6) = 2*TC* (1-RH);
                DI(6,6) = TD;
                if (L > 6)
                    TE = TD*RH;
                    DI(2,7) = -RH* (RH-1)* (RH-2)*(RH-3)*(RH-4)*(RH-5)/720;
                    DI(3,7) = TA* ((((62*RH-450)*RH+1190)*RH-1350)...
                        *RH+548)/720;
                    DI(4,7) = TB* (((-18*RH+75)*RH-102)*RH+45)/24;
                    DI(5,7) = TC* ((13*RH-30)*RH+17)/6;
                    DI(6,7) = 2.5*TD* (1-RH);
                    DI(7,7) = TE;
                    if (L>7)
                        TF = TE*RH;
                        DI(2,8) = RH*(RH-1)*(RH-2)*(RH-3)*(RH-4)*(RH-5)...
                            *(RH-6)/5040;
                        DI(3,8) = TA * ((((((-126*RH)+1302)*RH-5250)*...
                            RH+10290)*RH-9744)* RH+3528)/5040;
                        DI(4,8) = TB* ((((43*RH-270)*RH+625)*RH-630)...
                            *RH+232)/120;
                        DI(5,8) = TC* (((-10*RH+39)*RH-50)*RH+21)/6;
                        DI(6,8) = TD* ((20*RH-45)*RH+25)/6;
                        DI(7,8) = 3*TE*(1-RH);
                        DI(8,8) = TF;
                    end
                end
            end
        end
    end
end

for I = 1:N
    for J = 2:L
        ZZ = ZERO;
        for J1 = J:L
            ZZ = ZZ + DI(J,J1)*Y(I,J1);
        end
        Y(I,J) = ZZ;
    end
end

%========================================================================
%========================================================================

function [SOURCE, TARGET] = CPYARY(NELEM, SOURCE, TARGET)

TARGET(1:NELEM)=SOURCE(1:NELEM);

%========================================================================
%========================================================================


function [RH, H, OVRIDE] = HCHOSE(RH, H, OVRIDE)

global STPSZE HSTPSZ

% C     FIRST MOVE ALL ELEMENTS DOWN ONE PLACE
% C

if (H~=HSTPSZ(2,1))
    I = 12;
    while (I ~= 1)
        %for I=12:-1:2
        I2=I-1;
        HSTPSZ(1,I)=HSTPSZ(1,I2);
        HSTPSZ(2,I)=HSTPSZ(2,I2);
        I = I - 1;
    end
    % C
    % C     NOW INSERT VALUE OF H USED BEFORE THIS CALL
    % C
    HSTPSZ(1,2)=H/HSTPSZ(2,1);
    HSTPSZ(2,1)=H;
end
% C
% C NOW DECIDE ON THE NEW CHANGE
% C
if (RH>1)
    %OVRIDE=.FALSE.
    OVRIDE = 0;

elseif (HSTPSZ(1,2)<=1)
    %OVRIDE=.FALSE.
    OVRIDE = 0;
elseif ((RH*H)<=HSTPSZ(2,2))
    %OVRIDE=.FALSE.
    OVRIDE = 0;
else
    RH=HSTPSZ(2,2)/H;
    %OVRIDE=.TRUE.
    OVRIDE = 1;
end
HSTPSZ(1,1)=RH;

