with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Abstract_Ring;
with Generic_Vectors;
with Generic_VecVecs;
with Generic_Matrices;
with Generic_Monomials;
with Generic_Monomial_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Monomials is new Generic_Monomials(Ring,Vectors,VecVecs);
  with package Polynomials is
    new Generic_Monomial_Vectors(Ring,Vectors,VecVecs,Monomials);

package Generic_Polynomial_Vectors is

-- DESCRIPTION :
--   A polynomial vector is a vector of monomial vectors.

  type Polynomial_Vector is 
    array ( integer32 range <> ) of Polynomials.Link_to_Polynomial;

  type System ( dim,nbr : integer32 ) is record
   -- dim is the ambient dimension, the total number of variables
   -- nbr is the number of monomials, excluding the constant term
    pols : Polynomial_Vector(1..nbr); -- polynomials in the system
    deg1 : boolean;                   -- all variables appear with degree one?
    maxexp : Standard_Natural_Vectors.Vector(1..dim);
   -- maxexp(i) stores the largest exponent of the i-th variable
    powtab : VecVecs.Link_to_VecVec;
   -- if deg1, then powtab caches the power table of the variables
   -- which appear with power larger than one
  end record;

  type Link_to_Polynomial_Vector is access Polynomial_Vector;
  type Link_to_System is access System;

-- CONSTRUCTORS :

  procedure Power_Update ( s : in out System );
  procedure Power_Update ( s : in out Link_to_System );

  -- DESCRIPTION :
  --   Based on the monomials stored in p, determines the flag p.deg1,
  --   the largest exponents in p.maxexp.

  function Compute_Deg1 ( p : Polynomial_Vector ) return boolean;
  function Compute_Deg1 ( p : Link_to_Polynomial_Vector ) return boolean;
  function Compute_Deg1 ( s : System ) return boolean;
  function Compute_Deg1 ( s : Link_to_System ) return boolean;

  -- DESCRIPTION :
  --   The deg1 flag of a polynomial vector is true if no exponents
  --   higher than one appear.  Otherwise, deg1 is false.

  -- REQUIRED :
  --   For every polynomial in p or in s.pols, the deg1 flag is properly set.

  procedure Largest_Exponents
              ( p : in Polynomial_Vector;
                e : out Standard_Natural_Vectors.Vector );
  procedure Largest_Exponents
              ( p : in Link_to_Polynomial_Vector;
                e : out Standard_Natural_Vectors.Vector );
  procedure Largest_Exponents
              ( s : in System;
                e : out Standard_Natural_Vectors.Vector );
  procedure Largest_Exponents
              ( s : in Link_to_System;
                e : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in e the largest exponents of the variables in p.

  -- REQUIRED :
  --   For every polynomial in p or s.pols, maxexp is properly set.

  function Allocate_Power_Table
              ( maxexp : Standard_Natural_Vectors.Vector )
              return VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns an allocated power table for maximal exponents in maxexp.
  --   The vector on return has the same range as maxexp.

  procedure Allocate_Power_Table
              ( powtab : in out VecVecs.VecVec;
                maxexp : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Allocates memory for the power table for maximal exponents maxexp.

  -- REQUIRED : powtab'range = maxexp'range.

  procedure Update_Power_Table
              ( powtab : in out VecVecs.VecVec; x : in Vectors.Vector );
  procedure Update_Power_Table
              ( powtab : in VecVecs.Link_to_VecVec; x : in Vectors.Vector );

  -- DESCRIPTION :
  --   Given an allocated power table for the maximal exponents,
  --   updates the powers with new values x for the variables.

  -- REQUIRED : powtab'range = x'range and 
  --   powtab is allocated: powtab(row) runs from 0 to maxexp(row),
  --   where maxexp was used to allocate the power table.

-- EVALUATION and DIFFERENTIATION :

  function Eval ( p : Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector;
  function Eval ( s : System; x : Vectors.Vector ) return Vectors.Vector;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate p or s at x.

  function Eval ( p : Link_to_Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector;
  function Eval ( s : Link_to_System;
                  x : Vectors.Vector ) return Vectors.Vector;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate p or s at x.

  -- REQUIRED : p /= null and s /= null.

  procedure Diff ( p : in Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix );
  procedure Diff ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix );
  procedure Diff ( s : in System; x : in Vectors.Vector;
                   m : out Matrices.Matrix );
  procedure Diff ( s : in Link_to_System; x : in Vectors.Vector;
                   m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate all partial
  --   derivatives of p or s at x.  On return is the a matrix m.

  -- REQUIRED : m'range(1) = p'range = s.pols'range and m'range(2) = x'range.

  procedure Speel ( s : in System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( s : in System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix );
  procedure Speel ( s : in Link_to_System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( s : in Link_to_System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the Speelpenning algorithm to evaluate all partial
  --   derivatives of s at x, returned in m and p(x) is in y.

  -- ON ENTRY :
  --   s        system of polynomials in several variables;
  --   x        a vector of values for the variables of s.

  -- ON RETURN :
  --   y        the values of the polynomials in s at x;
  --   m        all partial derivatives of s at x;
  --   yd       vector with the partial derivatives of the current polynomial;
  --   wrk      work space used in the computation of the yd vectors.

  procedure Speel_on_Product
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector );
  procedure Speel_on_Product
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix );
  procedure Speel_on_Product
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector );
  procedure Speel_on_Product
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the Speelpenning algorithm to evaluate all partial
  --   derivatives of s at x, returned in m and p(x) is in y.

  -- REQUIRED : s.deg1 is true, no variables appear with exponent >= 2.
  --   The ranges of yd and wrk are equal to x'range.

  -- ON ENTRY :
  --   s        system of polynomials in several variables;
  --   x        a vector of values for the variables of s.

  -- ON RETURN :
  --   y        the values of the polynomials in s at x;
  --   m        all partial derivatives of s at x;
  --   yd       vector with the partial derivatives of the current polynomial;
  --   wrk      work space used in the computation of the yd vectors.

  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix );
  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the Speelpenning algorithm to evaluate all partial
  --   derivatives of p at x, returned in m and p(x) is in y.

  -- ON ENTRY :
  --   p        polynomials in several variables;
  --   x        a vector of values for the variables of p;
  --   powtab   contains the powers of the values of the variables,
  --            needed for the higher degree powers in the monomial,
  --            used to compute the common factor b for each monomial in p.

  -- ON RETURN :
  --   y        the values of the polynomials in p at x;
  --   m        all partial derivatives of p at x;
  --   yd       vector with the partial derivatives of the current polynomial;
  --   wrk      work space used in the computation of the yd vectors.

  procedure Speel_without_System_Cache
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector );
  procedure Speel_without_System_Cache
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix );
  procedure Speel_without_System_Cache
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector );
  procedure Speel_without_System_Cache
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the Speelpenning algorithm to evaluate all partial
  --   derivatives of s at x, returned in m and p(x) is in y,
  --   without using the cached powers of variables for the system.
  --   These procedures are for debugging purposes only.

  -- ON ENTRY :
  --   s        system of polynomials in several variables;
  --   x        a vector of values for the variables of s.

  -- ON RETURN :
  --   y        the values of the polynomials in s at x;
  --   m        all partial derivatives of s at x;
  --   yd       vector with the partial derivatives of the current polynomial;
  --   wrk      work space used in the computation of the yd vectors.

  procedure Speel_without_System_Cache
              ( p : in Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector );
  procedure Speel_without_System_Cache
              ( p : in Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix );
  procedure Speel_without_System_Cache
              ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector );
  procedure Speel_without_System_Cache
              ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the Speelpenning algorithm to evaluate all partial
  --   derivatives of p at x, returned in m and p(x) is in y,
  --   without using the cached powers of variables for the system.
  --   These procedures are for debugging purposes only.

  -- ON ENTRY :
  --   p        polynomials in several variables;
  --   x        a vector of values for the variables of p.

  -- ON RETURN :
  --   y        the values of the polynomials in p at x;
  --   m        all partial derivatives of p at x;
  --   yd       vector with the partial derivatives of the current polynomial;
  --   wrk      work space used in the computation of the yd vectors.

-- DESTRUCTORS :

  procedure Clear ( p : in out Polynomial_Vector );
  procedure Clear ( p : in out Link_to_Polynomial_Vector );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the polynomial vector p.

  procedure Clear ( s : in out System );
  procedure Clear ( s : in out Link_to_System );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the system s.

end Generic_Polynomial_Vectors;
