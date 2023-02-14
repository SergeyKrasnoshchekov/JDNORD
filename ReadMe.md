This demonstration program pre-calculates binomial coefficients and Stirling numbers
of the first kind and then invokes two routines for testing the correctness of the
the reduction formula for six-term products of angular momentum ladder operators to the
normal form:

     a     b     c     d     e     f    N           k(j)     l(j)     m(j)
J~6 = J(z)  J(+)  J(-)  J(z)  J(+)  J(-)  = SUM A(j) J(z)     J(+)     J(-)    ;          (1)
                                            j=1

and the reduction formula for the normal-ordered product of rotational ladder operators
and the Wigner function D1(0,eps), eps = -1,0,+1 to the normal form:

    a     b     c              N           (abcd)             r(ij)     s(ij)     t(ij)
J(z)  J(+)  J(-)  D1(0,eps) = SUM   SUM   C       D1(0,j) J(z)      J(+)      J(-)      . (2)
                              i=1 j=-1..1  ij

Both tests are based on comparing the results of action of operators before and after 
normal ordering on rigid rotor rotational wave functions |J,K>, J=0..J(max), K=-J..+J.

During these tests very small relative errors can occur due to the limited accuracy of
floating point numbers (currently Real(16)). An individual test is marked as `dubious'
if the difference between the supposedly exact results is not zero, but a small number.

The user can specify the max powers a,b,c,d,e,f (MAXPOW) and J(max) (JMAX).
As the derived coefficients can become very large, the recommended values for case (1) 
are 6/9 and for case (2) are 8/14.

LIBRARY ROUTINES:
BINOM_TABLE -- Pre-calculate binomial coefficients.
STIRL_TABLE -- Pre-calculate Stirling numbers of the first kind.
JREDUCE     -- Reduction of six-term products of angular momentum ladder operators to the
               normal form;
ROTOPKET    -- Evaluate the result of action of a rotational normal-ordered operator 
               product on rigid rotor rotational wave functions |J,K>;
DREDUCE     -- Reduction of the normal-ordered product of rotational ladder operators
               and the Wigner function D1(0,eps), eps = -1,0,+1 to the normal form:
WIGOPKET    -- Evaluate the result of action of a composite Wigner D-function and the 
               rotational normal-ordered operator product on wave functions |J,K>.

EXTERNAL ROUTINE:
Wigner3j    -- For comparison/verification of the correctness of Wigner 3-j symbols.

