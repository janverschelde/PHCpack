with Generic_Lists;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Binomial_Factors;          use Standard_Binomial_Factors;

package Standard_Puiseux_Certificates is

-- DESCRIPTION :
--   A Puiseux certificate for a regular common factor of two multivariate
--   polynomials consists in the exponent and coefficient of a second term
--   in the Puiseux series expansion of the factor.

-- DATA STRUCTURES :

  type Germ is record
    t : Term;
    c : Complex_Number;
    w : integer32;
  end record;

  package Germ_Lists is new Generic_Lists(Germ);
  type List_of_Germs is new Germ_Lists.List;

-- I. look for matching terms

  function Second_Power ( p : Poly; k,e : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the lowest power of the second variable with x(k)^e,
  --   where k equals the pivot of the pretropism, either 1 or 2.

  function Shift ( p : Poly; k : integer32 ) return Poly;

  -- DESCRIPTION :
  --   Returns the shifted polynomial so its initial form is a polynomial
  --   with nonzero constant term.  The k is the pivot of the pretropism,
  --   for k = 1, the standard direction is (1,0),
  --   for k = 2, the standard direction is (0,1).

  function Number_of_Initial_Roots
             ( p : Poly; k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of initial roots for a polynomial p, transformed
  --   using direction with povot equal to k and shifted properly, so the
  --   number on return is the largest degree in y (no x) when k = 1,
  --   or in case k = 2, the largest degree in x (no power in y).

  function Initial_Coefficients
             ( p : Poly; k : integer32 )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..Number_of_Initial_roots(p,k)
  --   with the coefficient of the initial form for p,
  --   for a transformed and properly shifted p.

  function Derivative ( p : Standard_Complex_Vectors.Vector;
                        x : Complex_Number ) return Complex_Number;

  -- DESCRIPTION : 
  --   Returns the value of the derivative at x of the univariate
  --   polynomial with coefficients in p.

  procedure Match ( f,g : in Poly; f1,g1 : in Complex_Number;
                    tol : in double_float; output : in boolean;
                    cff : out Complex_Number; exp : out integer32 );

  -- DESCRIPTION :
  --   Given are two univariate polynomials f, g and two terms with
  --   positive second exponent of the initial form.

  -- ON ENTRY :
  --   f,g      two univariate polynomials, transformed, shifted
  --            and evaluated for the second variable;
  --   f1,g1    derivatives of the initial forms of f and g
  --            evaluated at the initial root;
  --   tol      to decide whether a complex number is zero or not;
  --   output   if true, then intermediate output will appear.

  -- ON RETURN :
  --   cff      coefficient of matching term, zero if no match;
  --   exp      exponent of second term, zero if none exists.

  procedure Second_Term
              ( f,g : in Poly; t : in Term; tol : in double_float;
                output : in boolean;
                c : out Complex_Number; w : out integer32 );

  -- DESCRIPTION :
  --   Computes the second term of the Puiseux series expansion.

  -- ON ENTRY :
  --   f        a Laurent polynomial in two variables;
  --   g        another Laurent polynomial in two variables;
  --   t        potential initial term of a regular common factor;
  --   tol      tolerance to decide whether a number equals zero;
  --   output   if true, then intermediate output will print.

  -- ON RETURN :
  --   c        coefficient of the maching term, is zero if no match;
  --   w        exponent of the second term, if zero then nothing found.

  procedure Second_Terms
              ( f,g : in Poly; t : in List_of_Terms;
                tol : in double_float; output : in boolean;
                s : out List_of_Germs; fail : out boolean );

  -- DESCRIPTION :
  --   Computes the second term for a regular common factor of f and g,
  --   starting at the initial terms in t.

  -- ON ENTRY :
  --   f        a Laurent polynomial in two variables;
  --   g        another Laurent polynomial in two variables;
  --   t        potential initial terms of a regular common factor,
  --            the terms is t should not represent binomial factors;
  --   tol      tolerance to decide whether a number equals zero;
  --   output   set to true for intermediate output.

  -- ON RETURN :
  --   s        list of second terms;
  --   fail     if true, then all given terms in t are isolated roots,
  --            if false, then g contains certificates for a factor.

-- II. evaluate 

  function Evaluate ( t : Term; g : Germ ) return Poly;
  function Evaluate ( p : Poly; g : Germ ) return Poly;

  -- DESCRIPTION :
  --   Returns the value of the germ g at the term t or polynomial p,
  --   but only the leading and the second-order terms, not evaluating
  --   terms of order(t^(g.w^2)).

  function Order ( p : Poly; tol : double_float ) return integer32;

  -- DESCRIPTION :
  --   Returns the smallest degree over all terms t in p for which
  --   the coefficient of t is larger than tol in absolute value.

end Standard_Puiseux_Certificates;
