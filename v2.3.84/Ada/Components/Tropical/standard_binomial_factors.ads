with Generic_Lists;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Matrices;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;

package Standard_Binomial_Factors is

-- DESCRIPTION :
--   Given two Laurent polynomials in two variables, the computation of the
--   initial roots of the candidate tropisms suffice to identify the sparsest
--   nontrivial common factors, i.e.: those that are binomial. 

-- DATA STRUCTURE :

  package Term_Lists is new Generic_Lists(Term);
  type List_of_Terms is new Term_Lists.List;

-- I. compute tropisms

  function Common_Normals
               ( N,M : Standard_Integer64_Matrices.Matrix ) return List;

  -- DESCRIPTION :
  --   Returns the list of normals common to the columns of N and M.

  function Common_Inner_Normals
               ( A,B : Standard_Integer64_Matrices.Matrix ) return List;

  -- DESCRIPTION :
  --   Returns the common inner normals for the points with coordinates
  --   in the columns of A and B, sorted in lexicographic descreasing order.

  procedure Common_Inner_Normals
               ( A,B : in Standard_Integer64_Matrices.Matrix;
                 output : in boolean;
                 V,W,N,M,T : out List; fail : out boolean );

  -- DESCRIPTION :
  --   Returns the common inner normals to the points with coordinates in
  --   the columns of A and B with extra output.

  -- ON ENTRY :
  --   A         columns contain coordinates of points,
  --             lexicographically ordered in decreasing order;
  --   B         columns contain coordinates of points,
  --             lexicographically ordered in decreasing order;
  --   output    true if diagnostics have to be written on screen.

  -- ON RETURN :
  --   V        vertex set of A, two consecutive vertices span an edge;
  --   W        vertex set of B, two consecutive vertices span an edge;
  --   N        inner normals to V to the corresponding edges;
  --   M        inner normals to W to the corresponding edges;
  --   T        common inner normals of V and W;
  --   fail     true if checks with H-representations failed, false otherwise.

  function Common_Inner_Normals ( f,g : Poly ) return List;

  -- DESCRIPTION :
  --   Returns the inner normals common to two edges of respective
  --   Newton polygons of f and g.  The list on return contains the
  --   candidate tropisms for a common factor.

  procedure Common_Inner_Normals
              ( f,g : in Poly; output : in boolean;
                A,B,V,W,N,M,T : out List; fail : out boolean );

  -- DESCRIPTION :
  --   Returns in T the candidate tropisms for a common factor of f and g,
  --   with additional output for diagnostics.

  -- ON ENTRY :
  --   f        a Laurent polynomial in two variables;
  --   g        a Laurent polynomial in two variables;
  --   output   true if the routine should write to screen.

  -- ON RETURN :
  --   A        support of f in lexicographical decreasing order;
  --   B        support of g in lexicographical decreasing order;
  --   V        vertex set of A, two consecutive vertices span an edge;
  --   W        vertex set of B, two consecutive vertices span an edge;
  --   N        inner normals to V to the corresponding edges;
  --   M        inner normals to W to the corresponding edges;
  --   T        common inner normals of V and W;
  --   fail     true if checks with H-representations failed, false otherwise.

-- II. compute initial roots

  function Coefficients ( p : Poly ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficient vector of the univariate Laurent polynomial p
  --   as a vector of range Minimal_Degree(p,1)..Maximal_Degree(p,1).

  function Normalize ( c : Standard_Complex_Vectors.Vector )
                     return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   A normalized coefficient vector of a polynomial starts at 0.
  --   The vector on return is the coefficient vector of a polynomial
  --   with the same nonzero roots as the original Laurent polynomial.

  function Roots ( c : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector;
  function Roots ( c : Standard_Complex_Vectors.Vector;
                   max : natural32; tol : double_float )
                 return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the roots of the polynomial with normalized coefficient
  --   vector as computed by the method of Durand-Kerner.
  --   Default values for max and tol will be used if not provided.

  -- ON ENTRY :
  --   c        coefficient vector of a polynomials;
  --   max      maximum number of iterations in the root finder;
  --   tol      tolerance on the residual for the root finder.

  procedure Roots ( c : in Standard_Complex_Vectors.Vector;
                    max : in natural32; tol : in double_float;
                    z : out Standard_Complex_Vectors.Vector;
                    nb : out natural32; res : out double_float;
                    fail : out boolean );

  -- DESCRIPTION :
  --   Applies the Durand_Kerner method to find all roots of a polynomial,
  --   and returns more output than the other Roots functions above.

  -- REQUIRED : z'range = 1..c'last.
 
  -- ON ENTRY :
  --   c        coefficient vector of a polynomials;
  --   max      maximum number of iterations in the root finder;
  --   tol      tolerance on the residual for the root finder.

  -- ON RETURN :
  --   z        approximate values for the roots ;
  --   nb       number of Durand-Kerner iterations;
  --   res      maximal residual value of the roots;
  --   fail     true if tolerance is not achieved.


  function Match ( v : Standard_Integer_Vectors.Vector;
                   z1,z2 : Standard_Complex_Vectors.Vector;
                   tol : double_float ) return List_of_Terms;

  -- DESCRIPTION :
  --   Two roots match if they differ by less than the given tolerance.

  -- ON ENTRY :
  --   v        current tropism;
  --   z1       initial roots for first polynomial;
  --   z2       initial roots for second polynomial;
  --   tol      tolerance on difference between roots.

  -- ON RETURN :
  --   List of terms (may be empty) with matching initial roots
  --   as coefficient and tropism v as exponent vector.

  function Common_Initial_Roots
              ( f,g : Poly; v : Standard_Integer_Vectors.Vector;
                tol : double_float ) return List_of_Terms;

  -- DESCRIPTION :
  --   Computes the roots of the polynomials f and g and returns
  --   the matching roots as initial terms.

  -- REQUIRED : f and g are polynomials in one variable.

  procedure Common_Initial_Roots
              ( f,g : in Poly; v : in Standard_Integer_Vectors.Vector;
                output : in boolean;
                max : in natural32; eps,tol : in double_float;
                t : out List_of_Terms; fail : out boolean );

  -- DESCRIPTION :
  --   Computes roots of two univariate Laurent polynomials.

  -- ON ENTRY :
  --   f        a Laurent polynomial in one variable;
  --   g        another Laurent polynomial in one variable;
  --   v        candidate tropism;
  --   output   true if output should be written to screen;
  --   max      maximal number of iterations for root finder;
  --   eps      tolerance for the residual of the root finder;
  --   tol      tolerance for two roots to be consider equal.

  -- ON RETURN :
  --   t        list of initial terms for matching roots;
  --   fail     true if root finder failed.

  function Initial_Terms
              ( f,g : Poly; v : Standard_Integer_Vectors.Vector )
              return List_of_Terms;

  -- DESCRIPTION :
  --   Returns the initial terms with the candidate tropism v.
  --   For each common root of the initial forms, one term is returned
  --   with exponents equal to v and as coefficient the initial root.
  --   The list on return is empty if there are no common roots.

  procedure Initial_Terms
              ( f,g : in Poly; v : in Standard_Integer_Vectors.Vector;
                output : in boolean;
                max : in natural32; eps,tol : in double_float;
                t : out List_of_Terms; fail : out boolean );

  -- DESCRIPTION :
  --   Computes initial roots of f and g for the given tropism v.

  -- ON ENTRY :
  --   f        a Laurent polynomial in two variables;
  --   g        another Laurent polynomial in two variables;
  --   v        candidate tropism;
  --   output   true if output should be written to screen;
  --   max      maximal number of iterations for root finder;
  --   eps      tolerance for the residual of the root finder;
  --   tol      tolerance for two roots to be consider equal.

  -- ON RETURN :
  --   t        list of initial terms for matching roots;
  --   fail     true if root finder failed.

  function Initial_Terms ( f,g : Poly; t : List ) return List_of_Terms;

  -- DESCRIPTION :
  --   Computes the leading terms of Puiseux series expansions for
  --   common factors of the polynomials f and g.

  procedure Initial_Terms
              ( f,g : in Poly; v : in List; output : in boolean;
                max : in natural32; eps,tol : in double_float;
                t : out List_of_Terms; fail : out boolean );

  -- DESCRIPTION :
  --   Computes the leading terms of Puiseux series expansions for
  --   common factors of the polynomials f and g.

  -- ON ENTRY :
  --   f        a Laurent polynomial in two variables;
  --   g        another Laurent polynomial in two variables;
  --   v        a list of candidate tropisms;
  --   output   true if output should be written to screen;
  --   max      maximal number of iterations for root finder;
  --   eps      tolerance for the residual of the root finder;
  --   tol      tolerance for two roots to be consider equal.

  -- ON RETURN :
  --   t        list of initial terms for matching roots;
  --   fail     true if root finder failed.

-- III. evaluate

  function Evaluate ( t,x : Term ) return Term; 
  function Evaluate ( p : Poly; x : Term ) return Poly; 

  -- DESCRIPTION :
  --   Evaluates the term t or polynomial at the initial term x
  --   of the Puiseux series.
  --   The result is in general a polynomial in the free variable of
  --   the series and equals zero if x represents a binomial factor.

  function Residual ( p : Poly ) return double_float;

  -- DESCRIPTION :
  --   Returns the magnitude of the largest coefficient of p,
  --   as criterion to decide whether an initial term represents
  --   a binomial factor or not.

  procedure Split ( f,g : in Poly; t : in List_of_Terms;
                    tol : in double_float; b,r : out List_of_Terms );

  -- DESCRIPTION :
  --   Splits the list of terms t in two lists: b and r,
  --   in b the binomial factors and in r the rest.

  -- ON ENTRY :
  --   f,g      two Laurent polynomials in two variables;
  --   t        a list of initial terms of Puiseux series;
  --   tol      tolerance on the residual of evaluation.

  -- ON RETURN :
  --   b        all terms in t that evaluate to residual less than tol;
  --   r        remainder of t not representing binomial factors.

end Standard_Binomial_Factors;
