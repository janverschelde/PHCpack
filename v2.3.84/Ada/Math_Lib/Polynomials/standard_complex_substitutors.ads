with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Standard_Complex_Substitutors is

-- DESCRIPTION :
--   This package contains routines for substituting 
--   equations into polynomials and polynomial systems.

  function  Substitute ( k : integer32; c : Complex_Number; t : Term )
                       return Term;
  procedure Substitute ( k : in integer32; c : in Complex_Number; 
                         t : in out Term );

  function  Substitute ( k : integer32; c : Complex_Number; p : Poly )
                       return Poly;
  procedure Substitute ( k : in integer32; c : in Complex_Number;
                         p : in out Poly );

  function  Substitute ( k : integer32; c : Complex_Number; p : Poly_Sys )
                       return Poly_Sys;
  procedure Substitute ( k : in integer32; c : in Complex_Number;
                         p : in out Poly_Sys );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the term t or
  --   polynomial (system) p by a complex constant c.

  -- ON ENTRY :
  --   k         an unknown in the polynomial p;
  --   c         a complex constant;
  --   t         a term;
  --   p         a polynomial (system).

  -- ON RETURN :
  --   t         a term where the kth unknown has been replaced by the
  --             complex constant c;
  --   p         a polynomial (system) where the kth unknown has been
  --             replaced by the complex constant c.

  function  Substitute ( k : integer32; h : Vector; p : Poly ) return Poly;
  procedure Substitute ( k : in integer32; h : in Vector; p : in out Poly );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomial p
  --   by a linear equation defined by h.

  -- ON ENTRY :
  --   k          an unknown in the polynomial p;
  --   h          a vector h(0..n), n = Number_of_Unknowns(p),
  --              defines h(0) + h(1)*x1 + ... + h(n)*xn;
  --   p          a polynomial.

  -- REQUIRED : h(k) /= 0.

  -- ON RETURN :
  --   p          a polynomial where the kth unknown has been replaced
  --              by a linear equation defined by h.

  function  Substitute ( k : integer32; s,p : Poly ) return Poly;
  procedure Substitute ( k : in integer32; s : in Poly; p : in out Poly );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomial p
  --   by a polynomial s.

  -- ON ENTRY :
  --   k          an unknown in the polynomial p;
  --   s          a polynomial so that xk = s;
  --   p          a polynomial.

  -- ON RETURN :
  --   p          a polynomial where the kth unknown has been replaced
  --              by the polynomial s.

  function  Substitute ( k : integer32; h : Vector; p : Poly_Sys )
                       return Poly_Sys;
  procedure Substitute ( k : in integer32; h : in Vector;
                         p : in out Poly_Sys );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomials in the
  --   system p by a linear equation defined by h.

  -- ON ENTRY :
  --   k          an unknown in the polynomials in the system p;
  --   h          a vector h(0..n), n = Number_of_Unknowns(p(i)),
  --              defines h(0) + h(1)*x1 + ... + h(n)*xn;
  --   p          a polynomial system.

  -- REQUIRED : h(k) /= 0

  -- ON RETURN :
  --   p          a polynomial system where the kth unknown has been replaced
  --              by a linear equation defined by h.

  function  Substitute ( k : integer32; s : Poly;
                         p : Poly_Sys ) return Poly_Sys;
  procedure Substitute ( k : in integer32; s : in Poly; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomials in p
  --   by a polynomial s.

  -- ON ENTRY :
  --   k          an unknown in the polynomial p;
  --   s          a polynomial so that xk = s;
  --   p          a polynomial.

  -- ON RETURN :
  --   p          a polynomial system where the kth unknown has been replaced
  --              by the polynomial s.

end Standard_Complex_Substitutors;
