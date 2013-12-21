with C_Integer_Arrays,C_Double_Arrays;  use C_Integer_Arrays,C_Double_Arrays;

procedure phc_solver ( n,m : in integer; mc : in C_intarrs.Pointer;
                       ns : in integer; s : in C_intarrs.Pointer;
                       nc : in integer; c : in C_dblarrs.Pointer );

-- DESCRIPTION :
--   Given a coefficient-support representation of a polynomial system,
--   the blackbox solver of PHCpack is called.

-- ON ENTRY :
--   n         number of variables in every equation;
--   m         number of polynomials in the system;
--   mc        number of monomials in every polynomial, range 0..m-1;
--   ns        total number of exponents in the system;
--   s         support vector for the system, range 0..ns-1;
--   nc        number of coefficients;
--   c         coefficient vector for the system, range 0..nc-1.

