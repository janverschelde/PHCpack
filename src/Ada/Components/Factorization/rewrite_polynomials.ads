with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Rewrite_Polynomials is

-- DESCRIPTION :
--   Instead of polynomials of high degree, we better work with systems
--   of polynomials with modest degrees.  This package offers routines
--   to symbolically rewrite polynomials into lower degree systems.

  procedure Binary ( k,d : in natural32; deco : out Link_to_Vector );

  -- DESCRIPTION :
  --   Writes the binary decomposition of d on standard output.
  --   The step counter is k.  The binary decomposition is in deco,
  --   starting with the least significant bit.
  --   This routine was useful for developing, but is not as user
  --   friendly as the function Binary.

  -- REQUIRED : d > 0.

  function Binary_Length ( d : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the length of the binary representation of d,
  --   this is the number of bits plus one.

  function Recursive_Binary ( k,d : natural32 ) return Vector;

  -- DESCRIPTION :
  --   This is a recursive auxiliary procedure to compute the
  --   binary decomposition of the number d.  The k counts the
  --   number of divisions by 2.  The routine should be called
  --   with k = 0.

  -- REQUIRED : d > 0.

  function Binary ( d : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the binary decomposition of the natural number d.
  --   The result contains the coefficients of the binary respresentation
  --   of d, starting with the least significant bit.

  function Multi_Degrees
             ( p : Poly ) return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vectors with n-homogeneous degrees of the polynomial.

  procedure Number_of_Variables
               ( deg : in Standard_Natural_Vectors.Vector;
                 nvr : out Standard_Natural_Vectors.Vector;
                 tnv : out natural32 );

  -- DESCRIPTION :
  --   Given the multi-degrees of the polynomial p, this procedure returns
  --   the number of variables for each original variable in nvr, and the
  --   total number of variables in tnv.

  function Rewrite_Univariate_Term ( n : natural32; t : Term ) return Term;

  -- DESCRIPTION :
  --   Writes the term t as a term in n variables.

  -- REQUIRED : n >= length of binary decomposition of t.dg. 
  --            t is term in one variable.

  function Rewrite_Multivariate_Term
             ( n : natural32; t : Term;
               nvr : Standard_Natural_Vectors.Vector ) return Term;

  -- DESCRIPTION :
  --   Writes the multivariate term with auxiliary variables.
  --   For every original variable i there are nvr(i) auxiliary variables.
  --   The total number of variables equals n.

  function Rewrite_Univariate_Poly ( n : natural32; p : Poly ) return Poly;

  -- DESCRIPTION :
  --   Rewrites all terms of the polynomial, n is the number of
  --   variables of the new polynomial on return.

  function Rewrite_Multivariate_Poly
             ( n : natural32; p : Poly;
               nvr : Standard_Natural_Vectors.Vector ) return Poly;

  -- DESCRIPTION :
  --   Rewrites every term in the polynomial using auxiliary variables.

  procedure Telescope ( sys : in out Poly_Sys; n : in natural32 );

  -- DESCRIPTION :
  --   Fills in the equations x(i+1) - x(i)^2 into the system sys,
  --   where n is the number of variables used.  Thus i < n.

  procedure Telescope ( sys : in out Poly_Sys; n : in natural32;
                        nvr : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Defines the quadratic equations in the telescope, where n is
  --   the total number of variables and nvr(i) is the number of 
  --   variables in the i-th block.

  function Rewrite_Univariate_Polynomial ( p : Poly ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system that is equivalent to p.

  -- REQUIRED : p must be a univariate polynomial.

  function Rewrite_Multivariate_Polynomial ( p : Poly ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system that is equivalent to the multivariate
  --   polynomial p, using auxiliary variables.  Random hyperplanes are
  --   padded to the system to have as many equations as unknowns.

  procedure Enlarge_Symbol_Table ( n : in natural32; sb1 : in Symbol );

  -- DESCRIPTION :
  --   Enlarges the symbol table with n symbols of the form
  --   z1,..,zn, where z is contained in the symbol sb1.

  procedure Define_Symbol_Table
              ( n : in natural32; nvr : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Defines a double index symbol table with a total of n symbols
  --   grouped in blocks of nvr(i) variables.

end Rewrite_Polynomials;
