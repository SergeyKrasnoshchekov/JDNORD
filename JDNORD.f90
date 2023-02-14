!-----------------------------------------------------------------------------------------------
!  This project consisting of two files can be compiled by:
!     Intel Fortran (IFORT) with an essential option: Project Properties (Alt-Enter) > 
!     Configuration Properties > Fortran > Data > Default Real KIND := 16.
!  Alternatively, use Gfortran:
!     gfortran -c -O3 -fdefault-real-16 *.f90
!     gfortran *.o -o JDNORD
!
!  CONTENTS:
!
!  MODULE: PROBE_JREDUCE                                                       EDIT: 14 Feb 2023
!  MODULE: JREDUCE                                                             EDIT: 14 Feb 2023
!  MODULE: ROTOPKET                                                            EDIT: 14 Feb 2023
!  MODULE: PROBE_DREDUCE                                                       EDIT: 14 Feb 2023
!  MODULE: DREDUCE                                                             EDIT: 14 Feb 2023
!  MODULE: WIGOPKET                                                            EDIT: 14 Feb 2023
!  MODULE: BINOM_TABLE                                                         EDIT: 14 Feb 2023
!  MODULE: STIRL_TABLE                                                         EDIT: 14 Feb 2023
!-----------------------------------------------------------------------------------------------

      PROGRAM JDNORD  !  Pronounced: ['dzhee-'dai nord]
!-----------------------------------------------------------------------------------------------
!  INTRODUCTION:
!     This demonstration program pre-calculates binomial coefficients and Stirling numbers
!     of the first kind and then invokes two routines for testing the correctness of the
!     the reduction formula for six-term products of angular momentum ladder operators to the
!     normal form:
!
!               a     b     c     d     e     f    N           k(j)     l(j)     m(j)
! (1) J~6 = J(z)  J(+)  J(-)  J(z)  J(+)  J(-)  = SUM A(j) J(z)     J(+)     J(-)    ;
!                                                 j=1
!
!     and the reduction formula for the normal-ordered product of rotational ladder operators
!     and the Wigner function D1(0,eps), eps = -1,0,+1 to the normal form:
!
!         a     b     c              N           (abcd)             r(ij)     s(ij)     t(ij)
! (2) J(z)  J(+)  J(-)  D1(0,eps) = SUM   SUM   C       D1(0,j) J(z)      J(+)      J(-)     .
!                                   i=1 j=-1..1  ij
!
!     Both tests are based on comparing the results of action of operators before and after 
!     normal ordering on rigid rotor rotational wave functions |J,K>, J=0..J(max), K=-J..+J.
!
!     During these tests very small relative errors can occur due to the limited accuracy of
!     floating point numbers (currently Real(16)). An individual test is marked as `dubious'
!     if the difference between the supposedly exact results is not zero, but a small number.
!
!     The user can specify the max powers a,b,c,d,e,f (MAXPOW) and J(max) (JMAX).
!     As the derived coefficients can become very large, the recommended values for case (1) 
!     are 6/9 and for case (2) are 8/14.
!
!  CONTENTS:
! (A) BINOM_TABLE -- Pre-calculate binomial coefficients.
! (B) STIRL_TABLE -- Pre-calculate Stirling numbers of the first kind.
! (C) JREDUCE     -- Reduction of six-term products of angular momentum ladder operators to the
!                    normal form;
! (D) ROTOPKET    -- Evaluate the result of action of a rotational normal-ordered operator 
!                    product on rigid rotor rotational wave functions |J,K>;
! (E) DREDUCE     -- Reduction of the normal-ordered product of rotational ladder operators
!                    and the Wigner function D1(0,eps), eps = -1,0,+1 to the normal form:
! (F) WIGOPKET    -- Evaluate the result of action of a composite Wigner D-function and the 
!                    rotational normal-ordered operator product on wave functions |J,K>.
!-----------------------------------------------------------------------------------------------
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)
!----    Logical unit numbers for global use:
!        NINP -- File input (currently not used);
!        NOUT -- File output (general listing);
!        NAUX -- Auxiliary binary file(s);
!        NTMP -- Auxiliary temporary binary file(s);
!        NSTI -- Standard input (keyboard);
!        NSTO -- Standard output (console).
      INTEGER (4) :: NINP,NOUT,NAUX,NTMP,NSTI,NSTO
      COMMON /IOLUN/ NINP,NOUT,NAUX,NTMP,NSTI,NSTO

      INTEGER (4) :: BINOM(0:12,0:12)
      COMMON /BINOMTAB/ BINOM
      INTEGER (4) :: STIRL
      COMMON /STIRLTAB/ STIRL(0:12,0:12)

      CHARACTER (64) BUFFER
!----    Initialization of logical unit numbers:
      DATA NINP /1/, NOUT /2/, NAUX /3/, NTMP /4/, NSTI /5/, NSTO /6/
 
      RENUM = REAL (2) / REAL (3)
      WRITE (UNIT = BUFFER, FMT = "(F40.36)") RENUM
      N = 40 - 36
      IF (BUFFER(N+31:N+34) /= '6666') THEN
         WRITE (NSTO,"(A,A)") 'Actual Real model is not REAL(16): ',BUFFER
         WRITE (NSTO,"(A)") 'Fatal error, set IFORT/Gfortran option Default Real KIND = 16`'
         PAUSE
         STOP
      ENDIF

      WRITE (NSTO,1000)
 1000 FORMAT (/'Demonstration and verification of normal ordering of angular momentum'/        &
     &   'ladder operators and their products with Wigner D-functions.'/                       &
     &   'Partial results are shown on the screen, full output can be found in Listing.log')

      OPEN (UNIT = NOUT, FILE = 'Listing.log')

      CALL BINOM_TABLE
      CALL STIRL_TABLE

!----    Larger values are not recommended due to a potential overflow of Int(8) coefficients
      MAXPOW = 6
      JMAX   = 9
      CALL PROBE_JREDUCE (MAXPOW, JMAX, IERR)
      WRITE (NSTO,1100) IERR
 1100 FORMAT (/'#### PROBE_JREDUCE ####'/'Exit code = ',I0)

!----    Larger values are not recommended due to a potential overflow of Int(8) coefficients
      MAXPOW = 8
      JMAX   = 14
!----    Diangostic information is printed with this step
      NSTEP  = 5000
      CALL PROBE_DREDUCE (MAXPOW, JMAX, NSTEP, IERR)
      WRITE (NSTO,1200) IERR
 1200 FORMAT (/'#### PROBE_DREDUCE ####'/'Exit code = ',I0)

      CLOSE (1)
      WRITE (NSTO,'(A)') 'Press [Enter] for exit ... '
      PAUSE 
      STOP
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: PROBE_JREDUCE                                                       EDIT: 14 Feb 2023
!
!  PURPOSE:
!     Reduce a product of six powers of ladder rotational operators to the normal form:
!
!               a     b     c     d     e     f    N           k(j)     l(j)     m(j)
!     J~6 = J(z)  J(+)  J(-)  J(z)  J(+)  J(-)  = SUM A(j) J(z)     J(+)     J(-)    .
!                                                 j=1
!
!  THEORY:                        
!  1) MAIN REDUCTION FORMULA
!            d                j  min(c,e)           c!      k
!     J~6 = SUM C(d,j) (c - b)     SUM    C(e,k) --------  SUM S(k,l) x
!           j=0                    k=0           (c - k)!  l=0
!
!            l      m                l-m     (a+d-j+m)     (b+c-k)     (c+f-k)
!         x SUM (-2)  C(l,m) (2b-c+e)    J(z)          J(+)        J(-)       ,
!           m=0               
!
!     where:
!              / n \       n!     
!     C(n,k) = |   | = ----------- , and S(k,l) - Stirling numbers of the first kind.
!              \ k /   k! (n - k)!
!
!  2) MATRIX ELEMENTS
!     Evaluation of transformed rotational ket function | K, J > after the action of normally
!     ordered product of powers of rotational ladder operators:
!
!     J(z)^a * J(+)^b * J(-)^c | J, K > = Coeff | J', K' >,
!
!     where rotational ladder operators are defined in the following way from Cartesian form:
!     J(z) (ladder) = J(z), J(+) = J(x) - i J(y), J(-) = J(x) + i J(y).
!
!     The total angular momentum quantum number J can have values 0, +1, +2, +3, etc.
!     The additional quantum number K can have values -J, -J+1, -J+2, ... 0, 1, 2, ... +J.
!
!  a. Rising operator has the following definition: J(+) = J(x) - i J(y),
!     with the property of rising action on a rotational ket-function:
!                       _________________
!     J(+) | J, K > = \/ J(J+1) - K(K+1)  | J, K+1 >.
!
!  b. Lowering operator has the following definition: J(-) = J(x) + i J(y),
!     with the property of lowering action on a rotational ket-function:
!                       _________________
!     J(-) | J, K > = \/ J(J+1) - K(K-1)  | J, K-1 >.
!
!  c. Z-axis number operator coincides with its Cartesian form: J(z) (ladder) = J(z),
!     with the property of extracting angular quantum number K from a rotational ket-function:
!
!     J(z) | J,K > = K | J, K >.
!
!  INPUT:
!     MAXPOW -- Maximum power of operators J(+), J(-), J(z);
!     JMAX   -- Maximum value of rotational quantum number J.
!
!  OUTPUT:
!     IERR   -- Output error code, = 0 for success.
!-----------------------------------------------------------------------------------------------
      SUBROUTINE PROBE_JREDUCE (MAXPOW,JMAX,IERR)
      PARAMETER (MAXDIM = 1000)
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)

      INTEGER (4) :: NINP,NOUT,NAUX,NTMP,NSTI,NSTO
      COMMON /IOLUN/ NINP,NOUT,NAUX,NTMP,NSTI,NSTO

      INTEGER (4) :: NOPMAX,JJJJJJ(6)
      COMMON /STATS/ NOPMAX,JJJJJJ

      INTEGER (1) JOPROT(3),JPOW(6)
      INTEGER (4) KETROT(2),KETJJJ(2),KETJJJ1(2),KETJJJ2(2)
      INTEGER (1) JOPRED(3,MAXDIM)
      INTEGER (8) CONRED(MAXDIM)

      INTEGER (4) A,B,C,D,E,F

      DATA TOL /1.0q-32/      !  Tolerance for vanishing coefficients in DO-loop 100
      DATA TOLREL /1.0q-16/   !  Tolerance for `dubious' cases -- relative error

      NOPMAX = 0
      JJJJJJ = 0

      WRITE (NOUT,2000) MAXPOW
      WRITE (NOUT,2010) JMAX

      N = 0
      NERR = 0
      MAXTER = 0
      RELMAX = 0.0q0

      DO 400 J = 0, JMAX
      DO 400 K = -J, J
         KETROT(1) = J
         KETROT(2) = K

         N = N + 1
         WRITE (NSTO,1010) N,J,K
         WRITE (NOUT,1010) N,J,K

         JERR = 0
         M = 0
         DO 300 A = 0, MAXPOW
         DO 300 B = 0, MAXPOW
         DO 300 C = 0, MAXPOW

         DO 300 D = 0, MAXPOW
         DO 300 E = 0, MAXPOW
         DO 300 F = 0, MAXPOW
            M = M + 1

         DO 200 II = 1, 2
            IF (II == 2) WRITE (NOUT,1000) M,A,B,C,D,E,F
            IF (II == 2) WRITE (NOUT,1010) N,J,K

!----          Apply right part of J~6: J(z)^d J(+)^e J(-)^f to |J,K> and get C1 |J1,K1>
            JOPROT(1) = D
            JOPROT(2) = E
            JOPROT(3) = F
            CALL ROTOPKET (JOPROT,KETROT,COEF1,KETJJJ1)
            IF (II == 2) WRITE (NOUT,1020) JOPROT,KETROT,KETJJJ1,COEF1

!----          Apply left part of J~6: J(z)^a J(+)^b J(-)^c to |J1,K1> and get C2 |J2,K2>
            JOPROT(1) = A
            JOPROT(2) = B
            JOPROT(3) = C
            CALL ROTOPKET (JOPROT,KETJJJ1,COEF2,KETJJJ2)
            IF (II == 2) WRITE (NOUT,1030) JOPROT,KETJJJ1,COEF1,KETJJJ2,COEF1*COEF2

!----          Print the result: J~6 |J,K> = C1 C2 |J2,K2>
            COEF12 = COEF1 * COEF2
            IF (II == 2) WRITE (NOUT,1040) A,B,C,D,E,F,KETROT,KETJJJ2,COEF12

!----          Reduce J~6 to a linear combination of terms: SUM(j=1,N) c(j) J~3(j)
            JPOW(1) = A;  JPOW(4) = D  
            JPOW(2) = B;  JPOW(5) = E  
            JPOW(3) = C;  JPOW(6) = F  
            CALL JREDUCE (JPOW,MAXDIM,NTERM,JOPRED,CONRED,IRET)
            IF (IRET /= 0) GO TO 500
            MAXTER = MAX (MAXTER, NTERM)

            IF (II == 2) THEN
               WRITE (NOUT,1100) NTERM
            ENDIF

!----          SUM(j=1,N) c(j) J~3(j) |J,K> = SUM(j=1,N) c(j) C(j) |J',K'>
!              Verify: SUM(j=1,N) c(j) C(j) = C1 * C2  and  |J',K'> = |J2,K2>
            COEFSM = 0.0q0
            IEQUAL = 1
            DO 100 I = 1, NTERM
               JOPROT(:) = JOPRED(:,I)
               CALL ROTOPKET (JOPROT,KETROT,COEF,KETJJJ)
               COEFSM = COEFSM + COEF * REAL (CONRED(I), KIND = 16)

!----             If the trial ket |K,J> vanishes after JJJJJJ/JJJ(i), do not compare kets
!                 Example: J(-) J(+)^3  | 1, -1> = 0;
!                 Reduced: C(1) J(+)^3 J(-) + C(2) J(+)^2 + C(3) J(z) J(+)^2
               IF (ABS (COEF12) < TOL .OR. ABS (COEF) < TOL) CYCLE
               IF (KETJJJ(1) /= KETJJJ2(1) .OR. KETJJJ(2) /= KETJJJ2(2)) IEQUAL = 0

               IF (II == 2 .AND. IEQUAL == 0) THEN
                  WRITE (NOUT,1110) I,CONRED(I),(JOPRED(L,I),L=1,3)
                  WRITE (NOUT,1120) KETROT,CONRED(I),KETJJJ,CONRED(I)*COEF
                  WRITE (NOUT,1130) KETROT,KETJJJ,KETROT,KETJJJ2
               ENDIF
  100       CONTINUE

            DELTA = ABS (COEF12 - COEFSM)
            IF (COEF12 /= 0.0q0) THEN
               RELERR = DELTA / COEF12
               RELMAX = MAX (RELERR, RELMAX)
            ENDIF

            IF (II == 2) WRITE (NOUT,1200) DELTA,RELERR,IEQUAL,COEF12,COEFSM

!----          Skip second pass if no errors encountered
            IF (IEQUAL == 1 .AND. RELERR < TOLREL) THEN
               EXIT
            ELSE
               JERR = JERR + 1
               NERR = NERR + 1
            ENDIF
  200    CONTINUE  !  II = 1, 2

  300    CONTINUE  !  Loop over powers of JJJJJJ

         WRITE (NSTO,1210) JERR
         WRITE (NOUT,1210) JERR
  400 CONTINUE

      WRITE (NSTO,3000) M,N,NERR,100.0q0*REAL(NERR)/(REAL(M)*REAL(N)),MAXTER,RELMAX
      WRITE (NSTO,3100) MAXDIM,NOPMAX,MAXTER,JJJJJJ
      WRITE (NOUT,3000) M,N,NERR,100.0q0*REAL(NERR)/(REAL(M)*REAL(N)),MAXTER,RELMAX
      WRITE (NOUT,3100) MAXDIM,NOPMAX,MAXTER,JJJJJJ

      IF (RELMAX < TOLREL) THEN
         IERR = 0
         WRITE (NOUT,3200) 
      ELSE
         IERR = 1
         WRITE (NOUT,3210) 
      ENDIF

      RETURN
  500 IERR = -1
      RETURN

 2000 FORMAT (/96('-')//'#### PROBE_JREDUCE ####'/                                             &
     &   'Verification of reducing a product of ladder operators to normal form:'/             &
     &   'Operators: J(z)^a J(+)^b J(-)^c J(z)^d J(+)^e J(-)^f, where max(a,b,c,d,e,f) =',I3)
 2010 FORMAT (/'Trial wave functions:  | J = 0...J(max), K = -J...J >,  where J(max) = ',I0)

 1000 FORMAT (/96('-')//'Reduce a product number ',I8,' of ladder operators to normal form:'/  &
     &   'Operator J(z)^a J(+)^b J(-)^c J(z)^d J(+)^e J(-)^f,  where a,b,c,d,e,f = (',6I3,')')
 1010 FORMAT (/'Trial rotational wave function #',I0,':  | J = ',I0,', K = ',I0,' >')
 1020 FORMAT (/'Part 1: Operator J(z)^d J(+)^e J(-)^f,  where d,e,f = (',3I3,')'/              &
     &         'Original ket:       | J   = ',I0,', K   =',I3,' >'/                            &
     &         'Modified ket: Coeff | J`  = ',I0,', K`  =',I3,' >',',  Coeff =',E42.32)
 1030 FORMAT (/'Part 2: Operator J(z)^a J(+)^b J(-)^c,  where a,b,c = (',3I3,')'/              &
     &         'Original ket: Coeff | J`  = ',I0,', K`  =',I3,' >',',  Coeff =',E42.32/        &
     &         'Modified ket: Coeff | J`` = ',I0,', K`` =',I3,' >',',  Coeff =',E42.32)
 1040 FORMAT (/'Summary: J(z)^a J(+)^b J(-)^c J(z)^d J(+)^e J(-)^f, where a,b,c,d,e,f =',6I3/  &
     &         'Original ket:       | J   = ',I0,', K   =',I3,' >'/                            &
     &         'Modified ket: Coeff | J`` = ',I0,', K`` =',I3,' >',',  Coeff =',E42.32)

 1100 FORMAT (/'The number of reduced terms: ',I2)
 1110 FORMAT (/'Term #',I2,': ',F18.2,' x J(z)^k J(+)^l J(-)^m,  where k,l,m =',3I3)
 1120 FORMAT ( 'Original ket: Coeff | J   = ',I0,', K   =',I3,' >',',  Coeff =',E42.32/        &
     &         'Modified ket: Coeff | J`` = ',I0,', K`` =',I3,' >',',  Coeff =',E42.32)
 1130 FORMAT ( 'Verify  kets:'/  &
     &         'JJJJJJ  | J = ',I0,', K =',I3,' > ~ | J`` = ',I0,', K`` =',I3,' >'/            &
     &         'JJJ(i)  | J = ',I0,', K =',I3,' > ~ | J`` = ',I0,', K`` =',I3,' >')

 1200 FORMAT ('Delta = ',ES8.1,',  Relative error = ',ES8.1,',  Same Kets (0/1) = ',I0/        &
     &   'Coeff1 = ',ES40.32/'Coeff2 = ',ES40.32)
 1210 FORMAT ('Testing is accomplished, the number of dubious cases = ',I0)

 3000 FORMAT (/96('=')//'The number of studied operators =',I8/                                &
     &         'The number of studied kets |JK> =',I8/                                         &
     &         'The number of non-exact cases   =',I8,F12.6,' %'/                              &
     &         'Max number of reduced terms     =',I8/                                         &
     &         'Maximum relative deviation      =',ES12.1)
 3100 FORMAT (/'Limit of dimension of reduced terms before collection = ',I4/                  &
     &         'Max actual number of reduced terms before collection  = ',I4/                  &
     &         'Max actual number of reduced terms after  collection  = ',I4/                  &
     &         'The operator that caused largest number of terms      = ',6I4)
 3200 FORMAT (/'#### PROBE_JREDUCE ####'/'The test was successful.')
 3210 FORMAT (/'#### PROBE_JREDUCE ####'/'The test FAILED.')
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: JREDUCE                                                             EDIT: 14 Feb 2023
!
!  PURPOSE:
!     Reduce a product of six powers of ladder rotational operators to the normal form:
!
!               a     b     c     d     e     f    N           k(j)     l(j)     m(j)
!     J~6 = J(z)  J(+)  J(-)  J(z)  J(+)  J(-)  = SUM A(j) J(z)     J(+)     J(-)    .
!                                                 j=1
!
!  THEORY:                        
!     The general reduction formula:
!
!            d                j  min(c,e)           c!      k
!     J~6 = SUM C(d,j) (c - b)     SUM    C(e,k) --------  SUM S(k,l) x
!           j=0                    k=0           (c - k)!  l=0
!
!            l      m                l-m     (a+d-j+m)     (b+c-k)     (c+f-k)
!         x SUM (-2)  C(l,m) (2b-c+e)    J(z)          J(+)        J(-)       ,
!           m=0               
!
!     where:
!              / n \       n!     
!     C(n,k) = |   | = -----------  are binomial coefficients, and
!              \ k /   k! (n - k)!
!
!     S(k,l) are signed Stirling numbers of the first kind.
!
!  TECHNICAL:
!     Binomial coefficients and Stirling numbers of the first kind must be pre-calculated
!     and made available through /COMMON/ blocks, see the code.
!
!  INPUT:
!     JINP   -- Powers of six rotational operators: a,b,c,d,e,f;
!     MAXDIM -- Maximum dimension of output arrays COEF, JOUT;
!               equal to the maximum number of terms before collecting like terms.
!
!  OUTPUT:
!     NTERM  -- The number of generated terms;
!     JOUT   -- Array of triples (k,l,m) of powers of normilized rotational operators;
!     COEF   -- Integer(8): Coefficients of normalized terms, cannot be reduced to Integer(4);
!     IERR   -- Return code, = 0 for success; and 
!                            = 1 if the the number of generated terms exceeds MAXDIM.
!-----------------------------------------------------------------------------------------------
      SUBROUTINE JREDUCE (JINP,MAXDIM,NTERM,JOUT,COEF,IERR)
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)

!----    Initialized in STIRL_TABLE ()
      INTEGER (4) :: STIRL(0:12,0:12)
      COMMON /STIRLTAB/ STIRL
!----    Initialized in BINOM_TABLE ()
      INTEGER (4) :: BINOM(0:12,0:12)
      COMMON /BINOMTAB/ BINOM
!----    For statistical purpose only
      INTEGER (4) :: NOPMAX,JJJJJJ(6)
      COMMON /STATS/ NOPMAX,JJJJJJ

!----    Input/output variables
      INTEGER (1) JINP(6)
      INTEGER (1) JOUT(3,MAXDIM)
      INTEGER (8) COEF(MAXDIM)

!----    Internal variables
      INTEGER (4) A,B,C,D,E,F
!----    Coefficients can be large, use extended Integer (8):
      INTEGER (8) IPROD,IPROD1,IPROD2,IPROD3,IPROD4

!----    Clean output results from previous calls
      JOUT = 0
      COEF = 0

      A = JINP(1)
      B = JINP(2)
      C = JINP(3)
      D = JINP(4)
      E = JINP(5)
      F = JINP(6)

!----    Trivial cases
      NTERM = 1
      COEF(1) = 1
      IF (D + E + F == 0) THEN
         JOUT(1,1) = A
         JOUT(2,1) = B
         JOUT(3,1) = C
         GO TO 400
      ELSE IF (A == 0 .AND. B == 0 .AND. C == 0) THEN
         JOUT(1,1) = D
         JOUT(2,1) = E
         JOUT(3,1) = F
         GO TO 400
      ELSE IF (B == 0 .AND. C == 0) THEN
         JOUT(1,1) = A + D
         JOUT(2,1) = E
         JOUT(3,1) = F
         GO TO 400
      ELSE IF (C == 0 .AND. D == 0) THEN
         JOUT(1,1) = A
         JOUT(2,1) = B + E
         JOUT(3,1) = F
         GO TO 400
      ELSE IF (D == 0 .AND. E == 0) THEN
         JOUT(1,1) = A
         JOUT(2,1) = B
         JOUT(3,1) = C + F
         GO TO 400
      ENDIF

!----    General case
      N = 0
      COEF = 0

      DO 100 J = 0, D
         IPROD1 = BINOM(D,J)
         IF (J > 0) IPROD1 = IPROD1 * (C - B) ** J

      DO 100 K = 0, E  !  min(C,E)
         IF (K > C) CYCLE

         IPROD2 = 1
         DO I = 1, K
            IPROD2 = IPROD2 * (C - K + I)
         ENDDO
         IPROD2 = BINOM(E,K) * IPROD2

      DO 100 L = 0, K
         IPROD3 = STIRL(K,L)
         IF (IPROD3 == 0) CYCLE

      DO 100 M = 0, L
         IPROD4 = BINOM(L,M)
         IPROD4 = IPROD4 * (-2) ** M
         IF (M < L) IPROD4 = IPROD4 * (2 * B - C + E) ** (L - M)

!----       Assemble all multipliers together
         IPROD = IPROD1 * IPROD2 * IPROD3 * IPROD4
         IF (IPROD == 0) CYCLE

         IF (N == MAXDIM) GO TO 500
         N = N + 1
         COEF(N) = IPROD
         JOUT(1,N) = A + D - J + M
         JOUT(2,N) = B + E - K
         JOUT(3,N) = C + F - K
  100 CONTINUE

!----    Gathering statistics: memorize max number of terms before reduction and 
!        the operator responsible for this case
      IF (NOPMAX < N) THEN
         NOPMAX = N
         JJJJJJ = JINP
      ENDIF

!----    Collect like terms
      DO 300 I = 1, N - 1
         DO 200 J = I + 1, N
            IF (COEF(I) == 0 .OR. COEF(J) == 0) CYCLE
            DO L = 1, 3
               IF (JOUT(L,I) /= JOUT(L,J)) GO TO 200
            ENDDO
            COEF(I) = COEF(I) + COEF(J)
            COEF(J) = 0
  200    CONTINUE
  300 CONTINUE

      K = 0
      DO I = 1, N
         IF (COEF(I) == 0) CYCLE
         K = K + 1
         COEF(K) = COEF(I)
         IF (K == I) CYCLE
         DO J = 1, 3
            JOUT(J,K) = JOUT(J,I)
         ENDDO
      ENDDO
      NTERM = K

  400 IERR = 0
      RETURN
  500 IERR = 1
      RETURN
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: ROTOPKET                                                            EDIT: 14 Feb 2023
!
!  PURPOSE:
!     Evaluate transformed rotational ket function | K, J > after the action of normally
!     ordered product of powers of rotational ladder operators:
!
!     J(z)^a * J(+)^b * J(-)^c | J, K > = Coeff | J', K' >,
!
!     where rotational ladder operators are defined in the following way from Cartesian form:
!
!     J(z) (ladder) = J(z),
!     J(+) = J(x) - i J(y),
!     J(-) = J(x) + i J(y).
!
!     After step-by-step actions of rotational operators on the original ket function the
!     overall numerical coefficient is returned, along with the modified ket function.
!
!  THEORY:
!     The total angular momentum quantum number J can have values 0, +1, +2, +3, etc.
!     The additional quantum number K can have values -J, -J+1, -J+2, ... 0, 1, 2, ... +J.
!
!  1. Rising operator has the following definition:
!
!     J(+) = J(x) - i J(y),
!
!     with the property of rising action on a rotational ket-function:
!                       _________________
!     J(+) | J, K > = \/ J(J+1) - K(K+1)  | J, K+1 >.
!
!     The arithmetic expression J(J+1)-K(K+1) with 3 addition and 2 multiplication operators
!     can be written in an equivalent economic form (J-K)(J+K+1) with only one multiplication:
!                       _______________
!     J(+) | J, K > = \/ (J-K) (J+K+1)  | J, K+1 >.
!
!     It may be useful to have another form of this rule with the modified ket operator:
!                         _______________
!     J(+) | J, K-1 > = \/ (J+K) (J-K+1)  | J, K >.
!
!  2. Lowering operator has the following definition:
!
!     J(-) = J(x) + i J(y),
!
!     with the property of lowering action on a rotational ket-function:
!                       _________________
!     J(-) | J, K > = \/ J(J+1) - K(K-1)  | J, K-1 >.
!
!     The arithmetic expression J(J+1)-K(K-1) with 3 addition and 2 multiplication operators
!     can be written in an equivalent economic form (J-K)(J+K+1) with only one multiplication:
!                       _______________
!     J(-) | J, K > = \/ (J+K) (J-K+1)  | J, K-1 >.
!
!     It may be useful to have another form of this rule with the modified ket operator:
!                         _______________
!     J(-) | J, K+1 > = \/ (J-K) (J+K+1)  | J, K >.
!
!  3. Z-axis number operator coincides with its Cartesian form:
!
!     J(z) (ladder) = J(z),
!
!     with the property of extracting angular quantum number K from a rotational ket-function:
!
!     J(z) | J,K > = K | J, K >.
!
!  ALGORITHM:
!     At first, the numerical coefficient R is set equal to 1.0 and then
!     (1) the lowering operator J(-) is applied c times to ket function | J, K >
!         with the result R' | J', K' >;
!     (2) the rising operator J(+) is applied b times to ket function | J', K' >
!         with the result R'' | J'', K'' >;
!     (3) the number operator J(z) is applied a times to ket function | J'', K'' >
!         with the result R''' | J''', K''' >.
!
!     The calculation can be optimized using pre-calculated square roots of integer numbers.
!
!  INPUT:
!     JOPROT -- Powers (a,b,c) of ladder operators in J(z)^a J(+)^b J(-)^c;
!     KETROT -- Rotational quantum numbers J,K for the ket wave function |J,K>;
!
!  OUTPUT:
!     COEF   -- The real coefficient after applying the operator;
!     KETJJJ -- The quantum numbers of the new ket |J',K'> after applying the operator.
!-----------------------------------------------------------------------------------------------
      SUBROUTINE ROTOPKET (JOPROT,KETROT,COEF,KETJJJ)
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)

      REAL (16) :: COEF
      INTEGER (1) :: JOPROT(3)
      INTEGER (4) :: KETROT(2),KETJJJ(2)
      REAL (16) :: ARG

      COEF = 1.0q0
      J = KETROT(1)
      K = KETROT(2)

!----                      ______________
!        J(-) | J, K > = \/ (J+K)(J-K+1)  | J, K-1 >
      DO I = 1, JOPROT(3)
!----    Convert argument of SQRT() to Real(16) first to avoid the loss of accuracy
         ARG = (J + K) * (J - K + 1)
         IF (K - 1 < -J) GO TO 400
         K = K - 1
         COEF = COEF * SQRT (ARG)
      ENDDO

!----                      ______________
!        J(+) | J, K > = \/ (J-K)(J+K+1)  | J, K+1 >
      DO I = 1, JOPROT(2)
         ARG = (J + K + 1) * (J - K)
         IF (K + 1 > +J) GO TO 400
         K = K + 1
         COEF = COEF * SQRT (ARG)
      ENDDO

!----    J(z) | J, K > = K | J, K >
      DO I = 1, JOPROT(1)
         COEF = COEF * K
      ENDDO

      KETJJJ(1) = J
      KETJJJ(2) = K
      RETURN

!----   If K is out of range, return nil eigenvalue and set | J = 0, K = 0 >
  400 CONTINUE
      COEF = 0.0q0
      KETJJJ(1) = 0
      KETJJJ(2) = 0

      RETURN
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: PROBE_DREDUCE                                                       EDIT: 14 Feb 2023
!
!  PURPOSE:
!     Verify normal ordering of JJJ D products.
!
!  THEORY:
!         a     b     c               N           (abcd)             r(ij)     s(ij)     t(ij)
!     J(z)  J(+)  J(-)  D1(0,eps) = SUM   SUM   C       D1(0,j) J(z)      J(+)      J(-)      ,
!                                    i=1 j=-1..1  ij
!
!  ALGORITHM:
!     Compare results of action of JJJ D products on |JK> before and after normal ordering.
!
!  CALLS:
!     WIGOPKET (); DREDUCE (); ROTOPKET ();
!
!  INPUT:
!     MAXPOW -- Maximum individual powers of J(z)^a J(+)^b J(-)^c;
!     JMAX   -- Maximum values of a rotational quantyum number J. K = -J ... +J;
!     NSTEP  -- Output a case results with a certain step;
!
!  OUTPUT:
!     IERR   -- = 0, the test was fully successful; = 1, some cases failed.
!-----------------------------------------------------------------------------------------------
      SUBROUTINE PROBE_DREDUCE (MAXPOW,JMAX,NSTEP,IERR)
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)

      INTEGER (4) :: NINP,NOUT,NAUX,NTMP,NSTI,NSTO
      COMMON /IOLUN/ NINP,NOUT,NAUX,NTMP,NSTI,NSTO

      INTEGER (4) :: WOPMAX(-1:+1),KETMAX(2)
      COMMON /STATD/ RELE3J,WOPMAX,KETMAX

      INTEGER (1) :: WIGOPER(-1:+1),JOPER(3)
      INTEGER (4) :: KETROT(2),KETJJJ(2)
      INTEGER (4) :: KETWIG(2,-1:+1)  !  |J'(j),K'(j)>, j = -1, 0, +1
      REAL (16) :: WCOEF(-1:+1)

      INTEGER (1), ALLOCATABLE :: JOPERARR(:,:),WOPERARR(:,:)
      INTEGER (4), ALLOCATABLE :: KETOP1(:,:),KETOP2(:,:),KETOP3(:,:)
      REAL   (16), ALLOCATABLE :: DCOEF(:),DCOEF1(:),DCOEF2(:),DCOEF3(:)

      CHARACTER (64) BUFFER
      LOGICAL PASSED
!-----------------------------------------------------------------------------------------------
!     SETTINGS:
!     (a) Tolerance for comparing kets: TOL = 10^(-8);
!     (b) The maximum number of terms after normal ordering, MaxTerm = 256.
!-----------------------------------------------------------------------------------------------
      DATA TOL /1.0q-6/
      DATA MTERM /256/   !  Same setting as in DREDUCE

      RELE3J = 0.0q0

      WRITE (NOUT,1000)
      WRITE (NSTO,1010)

!----    Verify D-Wigner function
      DO IWIG = -1, +1
         WIGOPER = (/ 0,0,0 /)
         WIGOPER(IWIG) = 1

         WRITE (NOUT,1020) WIGOPER(-1), WIGOPER(0), WIGOPER(1)

         NCASE = 0
         DO JROT = 0, JMAX
         DO KROT = -JROT, JROT
            KETROT(1) = JROT
            KETROT(2) = KROT
            
            NCASE = NCASE + 1
            
            CALL WIGOPKET (WIGOPER,KETROT,KETWIG,WCOEF)
            
            WRITE (NOUT,1100) NCASE,WIGOPER(-1),WIGOPER(0),WIGOPER(1),JROT,KROT
            DO I = -1, 1
               IF (ABS (WCOEF(I)) == 0.0q0) CYCLE
               WRITE (UNIT = BUFFER, FMT = "(SP,F16.8)") WCOEF(I)
               WRITE (NOUT,1200) TRIM(ADJUSTL(BUFFER)),KETWIG(1,I),KETWIG(2,I)
            ENDDO
            WRITE (NOUT,1300)
         ENDDO; ENDDO

      ENDDO

!-----------------------------------------------------------------------------------------------
!     METHOD 1: Reduce the product JJJ D to a sum of operators SUM D' JJJ'
!
!     D^1(0,d) | J K > = SUM[eps=-1,0,+1] C_eps | J+eps K+d >.
!-----------------------------------------------------------------------------------------------
      OPEN (UNIT = NTMP, FORM = 'UNFORMATTED')

      ALLOCATE (JOPERARR(3,MTERM),WOPERARR(3,MTERM))
      ALLOCATE (DCOEF(MTERM))

      NCASE = 0
      NPASS = 0
      DEVMAX = 0.0q0
      RELMAX = 0.0q0

      DO 700 IWIG = -1, 1
         WIGOPER = (/ 0,0,0 /)
         WIGOPER(IWIG) = 1

      WRITE (NSTO,1400) IWIG,IWIG
      WRITE (NOUT,1400) IWIG,IWIG

      NOPER = 0

      DO 600 I1 = 0, MAXPOW
      DO 600 I2 = 0, MAXPOW
      DO 600 I3 = 0, MAXPOW
         JOPER(1) = I1
         JOPER(2) = I2
         JOPER(3) = I3
         NOPER = NOPER + 1

!----    Reduce non-normalized product into a sum of normalized products --
!        Jz^a J+^b J-^c D1(0,eps) = 
!        = SUM[i(a,b,c,eps)=1...n] SUM[j=-1...1] C(ij) D1(0,j) J(z)^r(ij) J(+)^s(ij) J(-)^t(ij)
!        =======================================================================================
         CALL DREDUCE (JOPER,WIGOPER, MTERM,NTERM, JOPERARR,WOPERARR,DCOEF)
!        =======================================================================================

      NPROB = 0
      DO 500 JROT = 0, JMAX
      DO 500 KROT = -JROT, JROT
         KETROT(1) = JROT
         KETROT(2) = KROT

         NPROB = NPROB + 1
         NCASE = NCASE + 1

         IF (MOD (NCASE,NSTEP) == 0) WRITE (NOUT,1410) NCASE,NOPER,IWIG,I1,I2,I3,NPROB,JROT,KROT
         IF (MOD (NCASE,NSTEP) == 0) WRITE (NSTO,1410) NCASE,NOPER,IWIG,I1,I2,I3,NPROB,JROT,KROT

         NRED = 0
         REWIND (NTMP)
         DO 200 I = 1, NTERM
!----       Apply i-th operator: C(ij) D1(0,j) J(z)^r(ij) J(+)^s(ij) J(-)^t(ij) to | J,K >

!----          | J K > transforms to: C_rot | J' K' >
            CALL ROTOPKET (JOPERARR(:,I),KETROT,RCOEF,KETJJJ)
!----          D^1(0,del) | J' K' > = SUM[eps=-1,0,+1] C_eps | J'+eps K'+del >
            CALL WIGOPKET (WOPERARR(:,I),KETJJJ,KETWIG,WCOEF)

!----          Save three transformed kets: SUM[eps=-1,0,+1] C_red C_rot C_eps | J'+eps K'+d >
            DO L = -1, 1
               COEF = DCOEF(I) * RCOEF * WCOEF(L)
               IF (ABS (COEF) < TOL) CYCLE
               WRITE (NTMP) COEF,KETWIG(1,L),KETWIG(2,L)
               NRED = NRED + 1
            ENDDO
  200    CONTINUE

!-----------------------------------------------------------------------------------------------
!        Collect like terms
!-----------------------------------------------------------------------------------------------
         REWIND (NTMP)
         ALLOCATE (DCOEF1(NRED),KETOP1(2,NRED))
         DO I = 1, NRED
            READ (NTMP) DCOEF1(I),KETOP1(1,I),KETOP1(2,I)
         ENDDO
         IF (NRED == 1) GO TO 230

         DO 210 M = 1, NRED - 1
         DO 210 N = M + 1, NRED
            IF (DCOEF1(M) == 0.0q0 .OR. DCOEF1(N) == 0.0q0) CYCLE
            DO L = 1, 2
               IF (KETOP1(L,M) /= KETOP1(L,N)) GO TO 210
            ENDDO
            DCOEF1(M) = DCOEF1(M) + DCOEF1(N)
            DCOEF1(N) = 0.0q0
  210    CONTINUE

         N = 0
         DO 220 I = 1, NRED
            IF (ABS (DCOEF1(I)) < TOL) CYCLE
            N = N + 1
            IF (N == I) CYCLE
            DCOEF1(N) = DCOEF1(I)
            DO L = 1, 2
               KETOP1(L,N) = KETOP1(L,I)
            ENDDO
  220    CONTINUE
         NRED = N

!----    Normal ordering: the result
  230    CONTINUE  
         IF (MOD (NCASE, NSTEP) == 0) THEN
            WRITE (NOUT,1500)
            IF (NRED == 0) WRITE (NOUT,1510)
            DO I = 1, NRED
               WRITE (NOUT,1520) DCOEF1(I), KETOP1(1,I), KETOP1(2,I)
            ENDDO
            WRITE (NSTO,1500)
            IF (NRED == 0) WRITE (NSTO,1510)
            DO I = 1, NRED
               WRITE (NSTO,1520) DCOEF1(I), KETOP1(1,I), KETOP1(2,I)
            ENDDO
         ENDIF

!===============================================================================================
!        Direct Calculation -- Jz^a J+^b J-^c D1(0,del)
!        Step 1: D1(0,del) | J,K > = SUM[eps=-1,0,+1] C(eps) | J+eps, K+del >
!        Step 2: Jz^a J+^b J-^c | J+eps, K+del >
!===============================================================================================
         CALL WIGOPKET (WIGOPER,KETROT,KETWIG,WCOEF)

         REWIND (NTMP)
         DO I = -1, 1
            CALL ROTOPKET (JOPER,KETWIG(:,I),RCOEF,KETJJJ)
            WRITE (NTMP) WCOEF(I) * RCOEF, KETJJJ
         ENDDO

         NDIR = 3
         ALLOCATE (DCOEF2(NDIR),KETOP2(2,NDIR))
         DCOEF2 = 0
         KETOP2 = 0

         REWIND (NTMP)
         DO I = 1, NDIR
            READ (NTMP) DCOEF2(I),KETOP2(1,I),KETOP2(2,I)
         ENDDO

!-----------------------------------------------------------------------------------------------
!        Collect like terms
!-----------------------------------------------------------------------------------------------
         DO 310 M = 1, NDIR - 1
         DO 310 N = M + 1, NDIR
            IF (DCOEF2(M) == 0.0q0 .OR. DCOEF2(N) == 0.0q0) CYCLE
            DO L = 1, 2
               IF (KETOP2(L,M) /= KETOP2(L,N)) GO TO 310
            ENDDO
            DCOEF2(M) = DCOEF2(M) + DCOEF2(N)
            DCOEF2(N) = 0.0q0
  310    CONTINUE

         N = 0
         DO 320 I = 1, NDIR
            IF (ABS (DCOEF2(I)) < TOL) CYCLE
            N = N + 1
            IF (N == I) CYCLE
            DCOEF2(N) = DCOEF2(I)
            DO L = 1, 2
               KETOP2(L,N) = KETOP2(L,I)
            ENDDO
  320    CONTINUE
         NDIR = N

         IF (MOD (NCASE, NSTEP) == 0) THEN
            WRITE (NOUT,1600)
            IF (NDIR == 0) WRITE (NOUT,1610)
            DO I = 1, NDIR
               WRITE (NOUT,1620) DCOEF2(I), KETOP2(1,I), KETOP2(2,I)
            ENDDO
            WRITE (NSTO,1600)
            IF (NDIR == 0) WRITE (NSTO,1610)
            DO I = 1, NDIR
               WRITE (NSTO,1620) DCOEF2(I), KETOP2(1,I), KETOP2(2,I)
            ENDDO
         ENDIF

!-----------------------------------------------------------------------------------------------
!        Compare by subtracting the first (ordered) operator from the unordered
!-----------------------------------------------------------------------------------------------
         ALLOCATE (DCOEF3(NRED+NDIR),KETOP3(2,NRED+NDIR))

         VALMAX = 0.0q0
         N = 0
         DO I = 1, NRED
            N = N + 1
            VALMAX = MAX (VALMAX, DCOEF1(I))
            DCOEF3(N) = DCOEF1(I)
            KETOP3(1,N) = KETOP1(1,I)
            KETOP3(2,N) = KETOP1(2,I)
         ENDDO
         DO I = 1, NDIR
            N = N + 1
            VALMAX = MAX (VALMAX, DCOEF2(I))
            DCOEF3(N) = -DCOEF2(I)
            KETOP3(1,N) = KETOP2(1,I)
            KETOP3(2,N) = KETOP2(2,I)
         ENDDO
         NTOT = N

         DO 410 M = 1, NTOT - 1
         DO 410 N = M + 1, NTOT
            IF (DCOEF3(M) == 0.0q0 .OR. DCOEF3(N) == 0.0q0) CYCLE
            DO L = 1, 2
               IF (KETOP3(L,M) /= KETOP3(L,N)) GO TO 410
            ENDDO
            DCOEF3(M) = DCOEF3(M) + DCOEF3(N)
            DCOEF3(N) = 0.0q0
  410    CONTINUE

         N = 0
         DO 420 I = 1, NTOT
            DEVMAX = MAX (DEVMAX, ABS (DCOEF3(I)))
            IF (VALMAX > 0.0q0) RELMAX = MAX (RELMAX, ABS (DCOEF3(I)/VALMAX))
            IF (ABS (DCOEF3(I)) < TOL) CYCLE
            N = N + 1
            IF (N == I) CYCLE
            DCOEF3(N) = DCOEF3(I)
            DO L = 1, 2
               KETOP3(L,N) = KETOP3(L,I)
            ENDDO
  420    CONTINUE
         NTOT = N

         IF (NTOT == 0) THEN
            PASSED = .TRUE.
         ELSE
            PASSED = .FALSE.
            WRITE (NOUT,1700) TOL
            DO 430 I = 1, NTOT
               WRITE (NOUT,1710) I,DCOEF3(I),KETOP3(1,I),KETOP3(2,I)
  430       CONTINUE
         ENDIF
         DEALLOCATE (DCOEF3,KETOP3)

         IF (PASSED) THEN
            NPASS = NPASS + 1
            IF (MOD (NCASE, NSTEP) == 0) WRITE (NOUT,1810)
            IF (MOD (NCASE, NSTEP) == 0) WRITE (NSTO,1810)
         ELSE
            PAUSE
            WRITE (NOUT,1820)
         ENDIF

         DEALLOCATE (KETOP1,DCOEF1)
         DEALLOCATE (KETOP2,DCOEF2)
  500 CONTINUE
  600 CONTINUE
  700 CONTINUE

      CLOSE (NTMP, STATUS = 'DELETE')

      DEALLOCATE (JOPERARR,WOPERARR,DCOEF)

      WRITE (NOUT,1900) MAXPOW,JMAX,NPASS,NCASE,100.0q0*REAL(NPASS)/REAL(NCASE),DEVMAX,RELMAX
      WRITE (NSTO,1900) MAXPOW,JMAX,NPASS,NCASE,100.0q0*REAL(NPASS)/REAL(NCASE),DEVMAX,RELMAX

      WRITE (NOUT,2000) RELE3J,WOPMAX,KETMAX
      WRITE (NSTO,2000) RELE3J,WOPMAX,KETMAX
 2000 FORMAT (/'Relative diff. between the 3-j symbol by library/local routines: ',ES12.2/  &
     &         'The particular Wigner D-function caused this result (-1, 0, +1): ',3I3/     &
     &         'The particular rigid rotor ket |JK> function caused this result: ',2I3)

      WRITE (NOUT,9000)
      WRITE (NSTO,9000)
 9000 FORMAT (/'PROBE_DREDUCE: Normal Termination.')

      IERR = 0
      IF (NPASS < NCASE) IERR = 1

      RETURN

 1000 FORMAT (/96('=')//                                                                       &
     &   'Normal ordering of the product of rotational ladder and Wigner D-operators:--'//     &
     &   'J(z)^a J(+)^b J(-)^c D^1_(0,d) ='/                                                   &
     &   'SUM [i=1..N] SUM[j=-1..1] C(abcd,ij) D(j) J(z)^r(ij) J(+)^s(ij) J(-)^t(ij)')
 1010 FORMAT (/'Verify normal ordering for rotational ladder and Wigner D-operators ...'/      &
     &   'Full report is provided with the step of ',I0,' test cases.')
 1020 FORMAT (/96('=')//'Verify action of a pure Wigner function D1(0,eps) on | J,K > '/       &
     &   'D1(0,-1)^',I1,' D1(0,0)^',I1,' D1(0,+1)^',I1//96('=')/)

 1100 FORMAT ('Case Number : ',I0/  &
     &   'Action of Wigner function D1(0,-1)^',I1,' D1(0,0)^',I1,' D1(0,+1)^',I1,' on ',       &
     &   '| J = ',I0,', K = ',I0,' >')
 1200 FORMAT (A,' | J = ',I0,', K = ',I0,' >')
 1300 FORMAT (96('-'))

 1400 FORMAT (/96('=')//                                                                       &
     &   'Normal ordering of the product of rotational ladder and Wigner D-operators:--'//     &
     &   'Case ',I0,' of -1, 0, +1'//'J(z)^a J(+)^b J(-)^c D1(0,',SP,I2,') ='//                &
     &   '= SUM[i=1..n] SUM[j=-1,0,+1] C(ij) D1(0,j) J(z)^r(ij) J(+)^s(ij) J(-)^t(ij)'/)
 1410 FORMAT (96('-')/'Global case: ',I0/                                                      &
     &   'Normal ordering of the product of rotational ladder and Wigner D-operators:--'/      &
     &   'Operator #',I0,':  J(z)^a J(+)^b J(-)^c D1(0,',I0,'), a,b,c =',3I3/                  &
     &   'Rot. ket #',I0,':  | J = ',I0,', K = ',I0,' >')

 1500 FORMAT (/'Transformed kets AFTER the normal ordering, SUM D1(0,eps) Jz^r J+^s J-^t:')
 1510 FORMAT ('   The operator is nil.')
 1520 FORMAT (ES40.32,' | J = ',I0,', K = ',I0,' > ':' +')

 1600 FORMAT ('Transformed kets BEFORE normal ordering of Jz^a J+^b J-^c D1(0,eps):')
 1610 FORMAT ('   The operator is nil.')
 1620 FORMAT (ES40.32,' | J = ',I0,', K = ',I0,' > ':' +')

 1700 FORMAT (/'Operator difference did not vanish!  Tolerance = ',ES8.1)
 1710 FORMAT (I4,ES40.32,' | J = ',I0,', K = ',I0' > ':' +')

 1810 FORMAT (/'SUCCESS: Test passed.')
 1820 FORMAT (/'ERROR: Test was not passed.')

 1900 FORMAT (/96('=')/  &
     &   'Testing of normal ordering of J(z)^a J(+)^b J(-)^c D^1_(0,d) is complete:'/          &
     &   'J_z^a J_+^b J_-^c, max(a,b,c) =',I2,', Rotational quantum number J_max = ',I0//      &
     &   'The results: ',I8,' of ',I8,' cases have passed test, i.e. ',F12.8,' %'/             &
     &   'Maximum non-vanishing difference coefficient = ',ES8.1/                              &
     &   'Max relative non-vanishing diff. coefficient = ',ES8.1)
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: DREDUCE                                                             EDIT: 14 Feb 2023
!
!  PURPOSE:
!     Reduce the operator product of ladder J-operators and Wigner D-function to normal form:
!
!         a     b     c              N           (abcd)             r(ij)     s(ij)     t(ij)
!     J(z)  J(+)  J(-)  D1(0,eps) = SUM   SUM   C       D1(0,j) J(z)      J(+)      J(-)
!                                   i=1 j=-1..1  ij
!     
!  INPUT:
!     JPOW   -- Array of powers (a,b,c) of angular momentum operator: Jz^a J+^b J+^c;
!     LPOW   -- Array of powers (d,e,f, d+e+f = 0 or 1) of the Wigner D^1_{-1,0,+1};
!     MAXDIM -- The max length of arrays for storing results, equal to the number of terms;
!
!  OUTPUT:
!     NTERM  -- The number of terms produced;
!     JOPER  -- (3,1:NTERM) Array of powers (a_j, b_j, c_j) of normally ordered operators:
!               Jz^a_j J+^b_j J+^c_j;
!     WOPER  -- (3,1:NTERM) Array of powers (d_j, e_j, f_j) of normally ordered D-operators;
!     DCOEF  -- (1:NTERM) Array of real coefficients of normally ordered operators;
!-----------------------------------------------------------------------------------------------
      SUBROUTINE DREDUCE (JPOW,LPOW,MAXDIM,NTERM,JOPER,WOPER,DCOEF)
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)

      INTEGER (4) :: NINP,NOUT,NAUX,NTMP,NSTI,NSTO
      COMMON /IOLUN/ NINP,NOUT,NAUX,NTMP,NSTI,NSTO

      INTEGER (4) :: BINOM(0:12,0:12)
      COMMON /BINOMTAB/ BINOM

      INTEGER (1) :: JPOW(3),LPOW(3)
      INTEGER (1) :: JOPER(3,MAXDIM),WOPER(3,MAXDIM)
      REAL (16) :: DCOEF(MAXDIM)
      INTEGER (4) :: WMIN,WZER,WPLU

      DATA WMIN,WZER,WPLU /1, 2, 3/

      IPOWA = JPOW(1)
      IPOWB = JPOW(2)
      IPOWC = JPOW(3)

      NTERM = 0
      JOPER = 0
      WOPER = 0
      COEF  = 0.0q0
      DCOEF = 0.0q0

      LSM = LPOW(1) + LPOW(2) + LPOW(3)
      IF (LSM == 0) THEN
         NTERM = 1
         JOPER(1,1) = JPOW(1)
         JOPER(2,1) = JPOW(2)
         JOPER(3,1) = JPOW(3)
         RETURN
      ENDIF
      IF (LSM > 1) THEN
         WRITE (NSTO,'(A)') '#### DREDUCE ####  Total power of D-function > 1'
         PAUSE
         STOP
      ENDIF
      IF (LPOW(1) == 1) GO TO 100
      IF (LPOW(2) == 1) GO TO 200
      IF (LPOW(3) == 1) GO TO 300

!-----------------------------------------------------------------------------------------------
!     Case 1:  Jz^a J+^b J-^c D1(0,-1)
!-----------------------------------------------------------------------------------------------
  100 CONTINUE
!----    Reduce (Jz - 1)^a to normal form
      DO 110 K = 0, IPOWA
         COEF = REAL (BINOM (IPOWA, K), KIND = 16)
         IF (MOD (IPOWA - K, 2) > 0) COEF = -COEF

         NTERM = NTERM + 1
         DCOEF(NTERM) = COEF
         WOPER(WMIN,NTERM) = 1  !  D^1_{0,-1}
         JOPER(1,NTERM) = K
         JOPER(2,NTERM) = IPOWB
         JOPER(3,NTERM) = IPOWC
  110 CONTINUE

      IF (IPOWB > 0) THEN
         NTERM = NTERM + 1
         DCOEF(NTERM) = - SQRT (2.0q0) * IPOWB
         WOPER(WZER,NTERM) = 1  !  D^1_{0,0}
         JOPER(1,NTERM) = IPOWA
         JOPER(2,NTERM) = IPOWB - 1
         JOPER(3,NTERM) = IPOWC
      ENDIF

!----    Reduce (Jz + 1)^a to normal form
      IF (IPOWB > 1) THEN
         IPROD = IPOWB * (IPOWB - 1)
         DO 120 K = 0, IPOWA
            COEF = IPROD * BINOM (IPOWA, K)
            NTERM = NTERM + 1
            DCOEF(NTERM) = COEF
            WOPER(WPLU,NTERM) = 1  !  D^1_{0,1}
            JOPER(1,NTERM) = K
            JOPER(2,NTERM) = IPOWB - 2
            JOPER(3,NTERM) = IPOWC
  120    CONTINUE
      ENDIF

      GO TO 800

!-----------------------------------------------------------------------------------------------
!     Case 2: Jz^a J+^b J-^c D1(0,0)
!-----------------------------------------------------------------------------------------------
  200 CONTINUE
!----    Reduce (Jz - 1)^a to normal form
      IF (IPOWC > 0) THEN
      DO 210 K = 0, IPOWA
         COEF = BINOM (IPOWA, K)
         IF (MOD (IPOWA - K, 2) > 0) COEF = -COEF

         NTERM = NTERM + 1
         DCOEF(NTERM) = - COEF * SQRT (2.0q0) * IPOWC
         WOPER(WMIN,NTERM) = 1
         JOPER(1,NTERM) = K
         JOPER(2,NTERM) = IPOWB
         JOPER(3,NTERM) = IPOWC - 1
  210 CONTINUE
      ENDIF

!----    Second term of two parts: D^1_0 (J_z^a J_+^b J_-^c + 2 b c J_z^a J_+^{b-1} J_-^{c-1})
      NTERM = NTERM + 1
      DCOEF(NTERM) = 1.0q0
      WOPER(WZER,NTERM) = 1
      JOPER(1,NTERM) = IPOWA
      JOPER(2,NTERM) = IPOWB
      JOPER(3,NTERM) = IPOWC

      IPROD = 2 * IPOWB * IPOWC
      IF (IPROD > 0) THEN
         NTERM = NTERM + 1
         DCOEF(NTERM) = IPROD
         WOPER(WZER,NTERM) = 1
         JOPER(1,NTERM) = IPOWA
         JOPER(2,NTERM) = IPOWB - 1
         JOPER(3,NTERM) = IPOWC - 1
      ENDIF

!----    Third term of two parts
      IF (IPOWB > 0) THEN
      DO 220 K = 0, IPOWA
         COEF = BINOM (IPOWA, K)

         NTERM = NTERM + 1
         DCOEF(NTERM) = - COEF * SQRT (2.0q0) * IPOWB
         WOPER(WPLU,NTERM) = 1
         JOPER(1,NTERM) = K
         JOPER(2,NTERM) = IPOWB - 1
         JOPER(3,NTERM) = IPOWC

         IF (IPOWB > 1) THEN
            NTERM = NTERM + 1
            IPROD = IPOWB * (IPOWB - 1) * IPOWC
            DCOEF(NTERM) = - COEF * SQRT (2.0q0) * IPROD
            WOPER(WPLU,NTERM) = 1
            JOPER(1,NTERM) = K
            JOPER(2,NTERM) = IPOWB - 2
            JOPER(3,NTERM) = IPOWC - 1
         ENDIF
  220 CONTINUE
      ENDIF

      GO TO 800

!-----------------------------------------------------------------------------------------------
!     Case 3:  Jz^a J+^b J-^c D1(0,+1)
!-----------------------------------------------------------------------------------------------
  300 CONTINUE
!----    Reduce (Jz - 1)^a to normal form
      IF (IPOWC > 0) THEN
      DO 310 K = 0, IPOWA
         COEF = BINOM (IPOWA, K)
         IF (MOD (IPOWA - K, 2) > 0) COEF = -COEF

         NTERM = NTERM + 1
         DCOEF(NTERM) = COEF * IPOWC * (IPOWC - 1)
         WOPER(WMIN,NTERM) = 1
         JOPER(1,NTERM) = K
         JOPER(2,NTERM) = IPOWB
         JOPER(3,NTERM) = IPOWC - 2
  310 CONTINUE
      ENDIF

!----    Second term of two parts:
!        D^1_0 2^1/2 c (J_z^a J_+^b J_-^{c-1} + b (c-1) J_z^a J_+^{b-1} J_-^{c-2})
      IF (IPOWC > 0) THEN
         NTERM = NTERM + 1
         DCOEF(NTERM) = - SQRT (2.0q0) * IPOWC
         WOPER(WZER,NTERM) = 1
         JOPER(1,NTERM) = IPOWA
         JOPER(2,NTERM) = IPOWB
         JOPER(3,NTERM) = IPOWC - 1
      ENDIF

      IF (IPOWB > 0 .AND. IPOWC > 1) THEN
         NTERM = NTERM + 1
         IPROD = IPOWB * IPOWC * (IPOWC - 1)
         DCOEF(NTERM) = - SQRT (2.0q0) * IPROD
         WOPER(WZER,NTERM) = 1
         JOPER(1,NTERM) = IPOWA
         JOPER(2,NTERM) = IPOWB - 1
         JOPER(3,NTERM) = IPOWC - 2
      ENDIF

!----    Third term of three parts:
!        D^1_1 2^1/2 c (J_z^a J_+^b J_-^{c-1} + b (c-1) J_z^a J_+^{b-1} J_-^{c-2})
      DO 320 K = 0, IPOWA
         COEF = BINOM (IPOWA, K)

         NTERM = NTERM + 1
         DCOEF(NTERM) = COEF
         WOPER(WPLU,NTERM) = 1
         JOPER(1,NTERM) = K
         JOPER(2,NTERM) = IPOWB
         JOPER(3,NTERM) = IPOWC

         IF (IPOWB > 0 .AND. IPOWC > 0) THEN
            NTERM = NTERM + 1
            IPROD = 2 * IPOWB * IPOWC
            DCOEF(NTERM) = COEF * IPROD
            WOPER(WPLU,NTERM) = 1
            JOPER(1,NTERM) = K
            JOPER(2,NTERM) = IPOWB - 1
            JOPER(3,NTERM) = IPOWC - 1
         ENDIF

         IF (IPOWB > 1 .AND. IPOWC > 1) THEN
            NTERM = NTERM + 1
            IPROD = IPOWB * (IPOWB - 1) * IPOWC * (IPOWC - 1)
            DCOEF(NTERM) = COEF * IPROD
            WOPER(WPLU,NTERM) = 1
            JOPER(1,NTERM) = K
            JOPER(2,NTERM) = IPOWB - 2
            JOPER(3,NTERM) = IPOWC - 2
         ENDIF
  320 CONTINUE

  800 RETURN
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: WIGOPKET                                                            EDIT: 14 Feb 2023
!
!  PURPOSE:
!     One of three possible Wigner D-functions, D^1(0,-1), D^1(0,0) or D^1(0,+1), acts on a
!     rotational ket | J,K > and produces three new rotational kets multiplied by real
!     coefficients.
!
!  INPUT:
!     WIGOPER - Array of powers of three operators: D^1(0,-1), D^1(0,0) and D^1(0,+1),
!               and operators are parametrized by single indices: [-1,0,+1].
!               Allowed combinations: 1,0,0; 0,1,0 or 0,0,1; or error message is generated.
!     KETROT -- Integer(4) array [1,2] of quantum numbers for | J,K > operator;
!
!  OUTPUT:
!     WCOEF  -- Real(16) Array [3]: coefficients for modified kets: |J'K'> _, |J'K'> 0, |J'K'> +
!     KETROTMOD Array [2,(-1,0,+1)] of modified quantum numbers for three kets.
!-----------------------------------------------------------------------------------------------
      SUBROUTINE WIGOPKET (WIGOPER,KETROT,KETROTMOD,WCOEF)
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)

      INTEGER (4) :: NINP,NOUT,NAUX,NTMP,NSTI,NSTO
      COMMON /IOLUN/ NINP,NOUT,NAUX,NTMP,NSTI,NSTO

      INTEGER (4) :: WOPMAX(-1:+1),KETMAX(2)
      COMMON /STATD/ RELE3J,WOPMAX,KETMAX

!----    Input parameters
      INTEGER (1) :: WIGOPER(-1:+1)
      INTEGER (4) :: KETROT(2)
!----    Output parameters
      INTEGER (4) :: KETROTMOD(2,-1:+1)  !  |J'(j),K'(j)>, j = -1, 0, +1
      REAL (16) :: WCOEF(-1:+1)
!----    Local parameters
      REAL (16) :: FACT2 = 1.0q0 / SQRT (2.0q0)
      REAL (16) :: W3JA(200)

      INTEGER (4) BRA_J,BRA_K,KET_J,KET_K
      LOGICAL CONDIT1, CONDIT2, CONDIT3

      DATA TOLERR /1.0q-16/

      J = KETROT(1)
      K = KETROT(2)

      WCOEF = 0.0q0
      KETROTMOD = 0

      N = WIGOPER(-1) + WIGOPER(0) + WIGOPER(+1)
      IF (N /= 0 .AND. N /= 1) THEN
         WRITE (NSTO,'(A)') '#### WIGOPKET #### Exactly one Wigner D(-1,0,+1) operator expected'
         PAUSE
         RETURN
      ENDIF

      IF (WIGOPER(-1) == 1) GO TO 100
      IF (WIGOPER( 0) == 1) GO TO 200
      IF (WIGOPER(+1) == 1) GO TO 300
      GO TO 900

!-----------------------------------------------------------------------------------------------
!     D^1_{0,-1}
!-----------------------------------------------------------------------------------------------
  100 CONTINUE
      DO 110 M = -1, 1
!----       Negative values of J are impossible
         IF (J + M < 0) THEN
            WCOEF(M) = 0.0q0
            KETROTMOD(1,M) = 0
            KETROTMOD(2,M) = 0
            CYCLE
         ENDIF

!----       Modified kets | J', K' > are produced regardless of the 3-j coefficients
         IF (ABS (K - 1) > J + M) THEN
            WCOEF(M) = 0.0q0
            KETROTMOD(1,M) = 0
            KETROTMOD(2,M) = 0
            CYCLE
         ELSE
            KETROTMOD(1,M) = J + M
            KETROTMOD(2,M) = K - 1
         ENDIF

!-----------------------------------------------------------------------------------------------
!        Prepare components of the 3-j symbol
!        / j1 j2 j3 \  /  1  J  J+M    \
!        \ m1 m2 m3 /  \ -1  K  -(K-1) /  
!-----------------------------------------------------------------------------------------------
         j1 =  1;  j2 = J;  j3 = J + M
         m1 = -1;  m2 = K;  m3 = -(K - 1)

!----       Common condition for all Wigner 3-j symbols for =/= 0:
         CONDIT1 = ABS (j1 + j2) >= j3
         CONDIT2 = j3 >= ABS (j1 - j2)
         CONDIT3 = ABS(m1) <= j1 .AND. ABS(m2) <= j2 .AND. ABS(m3) <= j3
         IF (.NOT. (CONDIT1 .AND. CONDIT2 .AND. CONDIT3)) CYCLE
!----       Evaluate Wigner 3-j
         IZERO = 0
         W3J = FACT2
         IF (M == -1 .AND. J > 0) THEN
            IZERO = 1
            IF (MOD (ABS (J-K), 2) > 0) W3J = -W3J
            ANOM = (J + K - 1) * (J + K)
            DENO = J * (2 * J - 1) * (2 * J + 1)
         ELSE IF (M == 0 .AND. J > 0) THEN
            IZERO = 1
            IF (MOD (ABS (J-K+1), 2) > 0) W3J = -W3J
            ANOM = (J - K + 1) * (J + K)
            DENO = J * (J + 1) * (2 * J + 1)
         ELSE IF (M == 1 .AND. J >= 0) THEN
            IZERO = 1
            IF (MOD (ABS (-J-K), 2) > 0) W3J = -W3J
            ANOM = (J - K + 1) * (J - K + 2)
            DENO = (J + 1) * (2 * J + 1) * (2 * J + 3)
         ENDIF
         IF (IZERO == 0) THEN
            W3J = 0.0q0
         ELSE
            W3J = W3J * SQRT (ANOM / DENO)
         ENDIF

!-----------------------------------------------------------------------------------------------
!        GitHub function: Wigner3j
!        Input: j2, j3, m1, m2, m3
!        Output: w3j[1 : jmax-jmin+1];  jmin = max(|j2-j3|, |m1|);  jmax = j2 + j3
!-----------------------------------------------------------------------------------------------
         CALL Wigner3j (W3JA, jmin, jmax, j2, j3, m1, m2, m3)
         l = 2 - jmin
         DIFF = ABS (W3JA(l) - W3J)
         IF (ABS (W3J) > 1.0q-16) THEN
            RELE3J = MAX (RELE3J, ABS (DIFF/W3J))
            WOPMAX = WIGOPER
            KETMAX = KETROT
         ENDIF
         IF (DIFF > TOLERR) THEN
            WRITE (NSTO,"('Case 1: 3-j diff (lib,our): ',2ES20.12,9I3)")  &
     &         W3JA(l),W3J,j1,j2,j3,m1,m2,m3,jmin,jmax,l
            PAUSE
         ENDIF
!----       This above portion of code is necessary as an independent check of 3-j symbols.
!           Can be deleted for actual applications after testing

         IF (W3J == 0.0q0) CYCLE
         ARG = 2 * J + 1
         PREFIX = SQRT (ARG)
         IF (MOD (ABS (K), 2) > 0) PREFIX = -PREFIX
         ARG = 2 * (J + M) + 1
         WCOEF(M) = PREFIX * W3J * SQRT (ARG)
  110 CONTINUE
      GO TO 900

!-----------------------------------------------------------------------------------------------
!     D^1_{0,0}
!-----------------------------------------------------------------------------------------------
  200 CONTINUE
      DO 210 M = -1, 1
!----       Negative values of J are impossible
         IF (J + M < 0) THEN
            WCOEF(M) = 0.0q0
            KETROTMOD(1,M) = 0
            KETROTMOD(2,M) = 0
            CYCLE
         ENDIF

!----       Modified kets | J', K' > are produced regardless of the 3-j coefficients
         KETROTMOD(1,M) = J + M
         KETROTMOD(2,M) = K

!-----------------------------------------------------------------------------------------------
!        Prepare components of the 3-j symbol
!        / j1 j2 j3 \  / 1  J  J+M \
!        \ m1 m2 m3 /  \ 0  K  -K  /  
!-----------------------------------------------------------------------------------------------
         j1 =  1;  j2 = J;  j3 = J + M
         m1 =  0;  m2 = K;  m3 = -K

!----       Common condition for all Wigner 3-j symbols for =/= 0:
         CONDIT1 = ABS (j1 + j2) >= j3
         CONDIT2 = j3 >= ABS (j1 - j2)
         CONDIT3 = ABS(m1) <= j1 .AND. ABS(m2) <= j2 .AND. ABS(m3) <= j3
         IF (.NOT. (CONDIT1 .AND. CONDIT2 .AND. CONDIT3)) CYCLE
         IF (K == 0 .AND. MOD (ABS (j1+j2+j3),2) > 0) CYCLE

!----       Evaluate Wigner 3-j
         IZERO = 0
         W3J = 1.0q0
         IF (M == -1 .AND. J > 0) THEN
            IZERO = 1
            IF (MOD (ABS (J-K), 2) > 0) W3J = -W3J
            ANOM = (J - K) * (J + K)
            DENO = J * (-1 + 2 * J) * (1 + 2 * J)
         ELSE IF (M == 0 .AND. J > 0) THEN
            IZERO = 1
            IF (MOD (ABS (J-K), 2) > 0) W3J = -W3J
            ANOM = K * K
!----          K is outside the square root, check its sign and modify W3J accordingly
            IF (K .LT. 0) W3J = -W3J
            DENO = J * (1 + J) * (1 + 2 * J)
         ELSE IF (M == 1 .AND. J >= 0) THEN
            IZERO = 1
            IF (MOD (ABS (1-J+K), 2) > 0) W3J = -W3J
            ANOM = (1 + J - K) * (1 + J + K)
            DENO = (1 + J) * (1 + 2 * J) * (3 + 2 * J)
         ENDIF
         IF (IZERO == 0) THEN
            W3J = 0.0q0
         ELSE
            W3J = W3J * SQRT (ANOM / DENO)
         ENDIF

!-----------------------------------------------------------------------------------------------
!        GitHub function: Wigner3j
!        Input: j2, j3, m1, m2, m3
!        Output: w3j[1 : jmax-jmin+1];  jmin = max(|j2-j3|, |m1|);  jmax = j2 + j3
!-----------------------------------------------------------------------------------------------
         CALL Wigner3j (W3JA, jmin, jmax, j2, j3, m1, m2, m3)
         l = 2 - jmin
         DIFF = ABS (W3JA(l) - W3J)
         IF (ABS (W3J) > 1.0q-16) THEN
            RELE3J = MAX (RELE3J, ABS (DIFF/W3J))
            WOPMAX = WIGOPER
            KETMAX = KETROT
         ENDIF
         IF (DIFF > TOLERR) THEN
            WRITE (NSTO,"('Case 2: 3-j diff (lib,our): ',2ES20.12,9I3)")  &
     &         W3JA(l),W3J,j1,j2,j3,m1,m2,m3,jmin,jmax,l
            PAUSE
         ENDIF
!----       This above portion of code is necessary as an independent check of 3-j symbols.
!           Can be deleted for actual applications after testing

         ARG = 2 * J + 1
         PREFIX = SQRT (ARG)
         IF (MOD (ABS(K), 2) > 0) PREFIX = -PREFIX
         ARG = 2 * (J + M) + 1
         WCOEF(M) = PREFIX * W3J * SQRT (ARG)
  210 CONTINUE
      GO TO 900

!-----------------------------------------------------------------------------------------------
!     D^1_{0,+1}
!-----------------------------------------------------------------------------------------------
  300 CONTINUE
      DO 310 M = -1, 1
!----       Negative values of J are impossible
         IF (J + M < 0) THEN
            WCOEF(M) = 0.0q0
            KETROTMOD(1,M) = 0
            KETROTMOD(2,M) = 0
            CYCLE
         ENDIF

!----       Modified kets | J', K' > are produced regardless of the 3-j coefficients
         IF (ABS (K + 1) > J + M) THEN
            WCOEF(M) = 0.0q0
            KETROTMOD(1,M) = 0
            KETROTMOD(2,M) = 0
            CYCLE
         ELSE
            KETROTMOD(1,M) = J + M
            KETROTMOD(2,M) = K + 1
         ENDIF

!-----------------------------------------------------------------------------------------------
!        Prepare components of the 3-j symbol
!        / j1 j2 j3 \  / 1  J  J+M    \
!        \ m1 m2 m3 /  \ 1  K  -(K+1) /  
!-----------------------------------------------------------------------------------------------
         j1 =  1;  j2 = J;  j3 = J + M
         m1 =  1;  m2 = K;  m3 = -(K + 1)

!----       Common condition for all Wigner 3-j symbols for =/= 0:
         CONDIT1 = ABS (j1 + j2) >= j3
         CONDIT2 = j3 >= ABS (j1 - j2)
         CONDIT3 = ABS(m1) <= j1 .AND. ABS(m2) <= j2 .AND. ABS(m3) <= j3
         IF (.NOT. (CONDIT1 .AND. CONDIT2 .AND. CONDIT3)) CYCLE

!----       Evaluate Wigner 3-j
         IZERO = 0
         W3J = FACT2
         IF (M == -1 .AND. J > 0) THEN
            IZERO = 1
            IF (MOD (ABS (J-K), 2) > 0) W3J = -W3J
            ANOM = (-1 + J - K) * (J - K)
            DENO = J * (-1 + 2 * J) * (1 + 2 * J)
         ELSE IF (M == 0 .AND. J > 0) THEN
            IZERO = 1
            IF (MOD (ABS (-J+K), 2) > 0) W3J = -W3J
            ANOM = (J - K) * (1 + J + K)
            DENO = J * (1 + J) * (1 + 2 * J)
         ELSE IF (M == 1 .AND. J >= 0) THEN
            IZERO = 1
            IF (MOD (ABS (-J+K), 2) > 0) W3J = -W3J
            ANOM = (1 + J + K) * (2 + J + K)
            DENO = (1 + J) * (1 + 2 * J) * (3 + 2 * J)
         ENDIF
         IF (IZERO == 0) THEN
            W3J = 0.0q0
         ELSE
            W3J = W3J * SQRT (ANOM / DENO)
         ENDIF

!-----------------------------------------------------------------------------------------------
!        GitHub function: Wigner3j
!        Input: j2, j3, m1, m2, m3
!        Output: w3j[1 : jmax-jmin+1];  jmin = max(|j2-j3|, |m1|);  jmax = j2 + j3
!-----------------------------------------------------------------------------------------------
         CALL Wigner3j (W3JA, jmin, jmax, j2, j3, m1, m2, m3)
         l = 2 - jmin
         DIFF = ABS (W3JA(l) - W3J)
         IF (ABS (W3J) > 1.0q-16) THEN
            RELE3J = MAX (RELE3J, ABS (DIFF/W3J))
            WOPMAX = WIGOPER
            KETMAX = KETROT
         ENDIF
         IF (DIFF > TOLERR) THEN
            WRITE (NSTO,"('Case 3: 3-j diff (lib,our): ',2ES20.12,9I3)")  &
     &         W3JA(l),W3J,j1,j2,j3,m1,m2,m3,jmin,jmax,l
            PAUSE
         ENDIF
!----       This above portion of code is necessary as an independent check of 3-j symbols.
!           Can be deleted for actual applications after testing

         ARG = 2 * J + 1
         PREFIX = SQRT (ARG)
         IF (MOD (ABS (K),2) > 0) PREFIX = -PREFIX
         ARG = 2 * (J + M) + 1
         WCOEF(M) = PREFIX * W3J * SQRT (ARG)
  310 CONTINUE

  900 CONTINUE

      RETURN
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: BINOM_TABLE                                                         EDIT: 14 Feb 2023
!
!  PURPOSE:
!     Calculate the table of binomial coefficients C(n,k), where n >= k >= 0.
!
!  THEORY:
!     The definition of the binomial coefficient C(n,k): 
!
!              / n \       n!
!     C(n,k) = |   | = -----------
!              \ k /   k! (n - k)!
!
!  EXAMPLE:
!            0     1     2     3     4     5     6     7     8     9    10    11    12
!     --------------------------------------------------------------------------------
!      0     1
!      1     1     1
!      2     1     2     1
!      3     1     3     3     1
!      4     1     4     6     4     1
!      5     1     5    10    10     5     1
!      6     1     6    15    20    15     6     1
!      7     1     7    21    35    35    21     7     1
!      8     1     8    28    56    70    56    28     8     1
!      9     1     9    36    84   126   126    84    36     9     1
!     10     1    10    45   120   210   252   210   120    45    10     1
!     11     1    11    55   165   330   462   462   330   165    55    11     1
!     12     1    12    66   220   495   792   924   792   495   220    66    12     1
!
!  ALGORITHM:
!     Binomial coefficients can be calculated recursively,
!
!     C(n,k) = C(n-1,k-1) + C(n-1,k)
!
!     with the initializations:
!
!     C(0,j) = 0 for j = 1..n and C(n,0) = 1 for i = 0...n
!
!  OUTPUT:
!     BINOM -- Table of the binomial coefficients, available through COMMON /BINOMTAB/.
!-----------------------------------------------------------------------------------------------
      SUBROUTINE BINOM_TABLE
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)
!----    Output table of binomial coefficients C(n,k)
      INTEGER (4) :: BINOM(0:12,0:12)
      COMMON /BINOMTAB/ BINOM

!----    Logical unit numbers for global use:
!        NINP -- File input (currently not used);
!        NOUT -- File output (general listing);
!        NAUX -- Auxiliary binary file(s);
!        NTMP -- Auxiliary temporary binary file(s);
!        NSTI -- Standard input (keyboard);
!        NSTO -- Standard output (console).
      INTEGER (4) :: NINP,NOUT,NAUX,NTMP,NSTI,NSTO
      COMMON /IOLUN/ NINP,NOUT,NAUX,NTMP,NSTI,NSTO

      DATA NTABLE /12/

      BINOM = 0
      BINOM(1,1) = 1
      BINOM(NTABLE,NTABLE) = 1

      DO J = 1, NTABLE
         BINOM(0,J) = 0
      ENDDO
      DO I = 0, NTABLE
         BINOM(I,0) = 1
      ENDDO

      DO I = 1, NTABLE
         DO J = 1, NTABLE
            BINOM(I,J) = BINOM(I-1,J-1) + BINOM(I-1,J)
         ENDDO
      ENDDO

      WRITE (NOUT,"(/'Pre-calculated binomial coefficients:')")
      DO I = 0, NTABLE
         WRITE (NOUT,'(I4,I6,12I6)') I,(BINOM(I,J),J=0,I)
      ENDDO

      RETURN
      END


!-----------------------------------------------------------------------------------------------
!  MODULE: STIRL_TABLE                                                         EDIT: 14 Feb 2023
!
!  PURPOSE:
!     Pre-calculate the Stirling number of the first kind:
!                                                           n
!     (x)_n  = x * (x - 1) * (x - 2) * ... * (x - n + 1) = SUM s(n, k) * x^k
!                                                          k=0
!     The coefficients s(n,k) can be calculated using the following recursive relation:
!
!     s(n, k) = s(n - 1, k - 1) - (n - 1) * s(n - 1, k);  s(0, 0) = 1.
!
!  EXAMPLE:
!            0           1           2           3           4           5           6         7
!      -----------------------------------------------------------------------------------------
!      0     1
!      1     0           1
!      2     0          -1           1
!      3     0           2          -3           1
!      4     0          -6          11          -6           1
!      5     0          24         -50          35         -10           1
!      6     0        -120         274        -225          85         -15           1
!      7     0         720       -1764        1624        -735         175         -21         1
!      8     0       -5040       13068      -13132        6769       -1960         322       -28
!      9     0       40320     -109584      118124      -67284       22449       -4536       546
!     10     0     -362880     1026576    -1172700      723680     -269325       63273     -9450
!     11     0     3628800   -10628640    12753576    -8409500     3416930     -902055    157773
!     12     0   -39916800   120543840  -150917976   105258076   -45995730    13339535  -2637558
!
!                        8           9          10          11          12
!     --------------------------------------------------------------------
!      8                 1                                           
!      9               -36           1                               
!     10               870         -45           1                   
!     11            -18150        1320         -55           1       
!     12            357423      -32670        1925         -66           1
!
!  TECHNICAL:
!     Intel® Fortran Compiler Classic and Intel® Fortran Compiler ... 
!     INTEGER(4) values range from -2,147,483,648 to 2,147,483,647
!
!  OUTPUT:
!     STIRL -- Table of the Stirling numbers, available through COMMON /STIRLTAB/.
!-----------------------------------------------------------------------------------------------
      SUBROUTINE STIRL_TABLE
      IMPLICIT REAL (16) (A-H,O-Z), INTEGER (4) (I-N)

!----    Output table of Stirling numbers, Integer (8) can be necessary for larger table!
      INTEGER (4) :: STIRL(0:12,0:12)
      COMMON /STIRLTAB/ STIRL

!----    Logical unit numbers for global use:
!        NINP -- File input (currently not used);
!        NOUT -- File output (general listing);
!        NAUX -- Auxiliary binary file(s);
!        NTMP -- Auxiliary temporary binary file(s);
!        NSTI -- Standard input (keyboard);
!        NSTO -- Standard output (console).
      INTEGER (4) :: NINP,NOUT,NAUX,NTMP,NSTI,NSTO
      COMMON /IOLUN/ NINP,NOUT,NAUX,NTMP,NSTI,NSTO

      DATA NTABLE /12/

      STIRL = 0
      STIRL(0,0) = 1

      DO I = 1, NTABLE
         DO J = 1, NTABLE
            IF (J > I) CYCLE
            STIRL(I,J) = STIRL(I-1,J-1) - (I-1) * STIRL(I-1,J)
         ENDDO
      ENDDO

      WRITE (NOUT,"(/'Pre-calculated Stirling numbers of the first kind:')")
      DO I = 0, NTABLE
         WRITE (NOUT,'(I4,I6,6I12,6I10)') I,(STIRL(I,J),J=0,I)
      ENDDO

      RETURN
      END

!-----------------------------------------------------------------------------------------------
!     End of file.
!-----------------------------------------------------------------------------------------------