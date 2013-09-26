with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;

package Reduction_of_Polynomials is

-- DESCRIPTION :
--   This package implements S-polynomials and R-polynomials.

  function Spoly ( p,q : poly ) return Poly;

  -- DESCRIPTION :
  --   Returns the S-polynomial of p and q :
  --              lcm(in(p),in(q))              lcm(in(p),in(q))
  --   S =  c_q * ----------------  p  -  c_p * ---------------- q
  --                   in(p)                          in(q)
  --   where lcm stands for the least common multiple,
  --         in(p) is the leading term of the polynomial p
  --     and the coefficients c_q and c_p are chosen such that
  --         their moduli are smaller than or equal to 1.

  function Rpoly ( p,q : Poly ) return Poly;

  -- DESCRIPTION :
  --   Returns the R-polynomial of the polynomials p and q :
  --            c_p   lcm(in(p),term(q))
  --   R = p -  --- * ------------------ * q
  --            c_q        term(q)
  --   such that the leading term of p vanishes.

end Reduction_of_Polynomials;
