with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Abstract_Ring;
with Generic_Vectors;
with Generic_VecVecs;
with Generic_Monomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Monomials is new Generic_Monomials(Ring,Vectors,VecVecs);

package Generic_Monomial_Vectors is

-- DESCRIPTION :
--   A monomial vector is a vector of monomials in several variables.

  type Monomial_Vector is 
    array ( integer32 range <> ) of Monomials.Link_to_Monomial;

  type Polynomial ( dim,nbr : integer32 ) is record
   -- dim is the ambient dimension, the total number of variables
   -- nbr is the number of monomials, excluding the constant term
    cff0 : Ring.number;             -- constant term of the polynomial
    mons : Monomial_Vector(1..nbr); -- monomials with at least one variable
    deg1 : boolean;                 -- all variables appear with degree one?
    maxexp : Standard_Natural_Vectors.Vector(1..dim);
   -- maxexp(i) stores the largest exponent of the i-th variable
    powtab : VecVecs.Link_to_VecVec;
   -- if deg1, then powtab caches the power table of the variables
   -- which appear with power larger than one
  end record;

  type Link_to_Monomial_Vector is access Monomial_Vector;
  type Link_to_Polynomial is access Polynomial;

-- CONSTRUCTORS :

  procedure Power_Update ( p : in out Polynomial );
  procedure Power_Update ( p : in out Link_to_Polynomial );

  -- DESCRIPTION :
  --   Based on the monomials stored in p, determines the flag p.deg1,
  --   the largest exponents in p.maxexp and, if not p.deg1, allocates
  --   the power table in p.powtab.

  function Compute_Deg1 ( v : Monomial_Vector ) return boolean;
  function Compute_Deg1 ( v : Link_to_Monomial_Vector ) return boolean;
  function Compute_Deg1 ( p : Polynomial ) return boolean;
  function Compute_Deg1 ( p : Link_to_Polynomial ) return boolean;

  -- DESCRIPTION :
  --   The deg1 flag of a monomial vector or polynomial is true
  --   if the n_base attribute for all monomials equals zero,
  --   otherwise, if variables appear with an exponent higher than one,
  --   then the deg1 flag is false.
  --   If deg1, then no power table is needed to evaluate and differentiate
  --   the monomials.

  procedure Largest_Exponents
              ( v : in Monomial_Vector;
                e : out Standard_Natural_Vectors.Vector );
  procedure Largest_Exponents
              ( v : in Link_to_Monomial_Vector;
                e : out Standard_Natural_Vectors.Vector );
  procedure Largest_Exponents
              ( p : in Polynomial;
                e : out Standard_Natural_Vectors.Vector );
  procedure Largest_Exponents
              ( p : in Link_to_Polynomial;
                e : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in e the largest exponents which occur in v or p.

  -- REQUIRED :
  --   The range of e should be 1..v(v'first).dim = p.dim.

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

-- SELECTORS :

  function Degree ( v : Monomial_Vector ) return integer32;
  function Degree ( v : Link_to_Monomial_Vector ) return integer32;
  function Degree ( p : Polynomial ) return integer32;
  function Degree ( p : Link_to_Polynomial ) return integer32;

  -- DESCRIPTION :
  --   Returns -1 for an empty monomial vector or polynomial,
  --   otherwise the degree of the vector or polynomial is returned.

-- EVALUATION and DIFFERENTIATION :

  function Eval ( v : Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number;
  function Eval ( v : Link_to_Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate v at x.

  function Eval ( p : Polynomial;
                  x : Vectors.Vector ) return Ring.number;
  function Eval ( p : Link_to_Polynomial;
                  x : Vectors.Vector ) return Ring.number;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate p at x.

  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector );
  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector );
  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector );
  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to compute all derivatives
  --   of the polynomial defined by the sum of the monomials in v at x.

  -- REQUIRED : the ranges for yd and wrk are equal to x'range.

  -- ON ENTRY :
  --   v        a monomial vector;
  --   x        a vector of range 1..m.dim.

  -- ON RETURN :
  --   yd       all partial derivatives of v at x;
  --   wrk      work space to hold intermediate values of the derivatives
  --            of the monomials (optional).

  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates p at x and computes all partial derivatives
  --   using the Speelpenning algorithm.

  -- REQUIRED :
  --   if not p.deg1, then the power table of p is allocated.

  -- ON ENTRY :
  --   p        a polynomial in several variables;
  --   x        a vector of range 1..p.dim.

  -- ON RETURN :
  --   y        the value of the monomial vector at x;
  --   yd       all partial derivatives of v at x;
  --   wrk      work space to hold intermediate values of the derivatives
  --            of the monomials (optional).

  procedure Speel_on_Product
                  ( v : in Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel_on_Product
                  ( v : in Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel_on_Product
                  ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel_on_Product
                  ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the monomials
  --   in v and all its derivatives at x.

  -- REQUIRED : v(i).nvr > 0, at least one variable has exponent 1 and
  --   v(i).n_base = 0, i.e.: no variables appear with exponent 2 or higher.
  --   The ranges of yd and wrk are equal to x'range.

  -- ON ENTRY :
  --   v        a monomial vector;
  --   x        a vector of range 1..v(i).dim.

  -- ON RETURN :
  --   y        the value of the monomial vector at x;
  --   yd       all partial derivatives of v at x;
  --   wrk      work space to hold intermediate values of the derivatives
  --            of the monomials.

  procedure Speel_on_Product
                  ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel_on_Product
                  ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel_on_Product
                  ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel_on_Product
                  ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the polynomial p
  --   and all its derivatives at x.

  -- REQUIRED : p.mons(i).nvr > 0, at least one variable has exponent 1 and
  --   p.mons(i).n_base = 0, i.e.: no variables appear with exponent >= 2.
  --   The ranges of yd and wrk are equal to x'range.

  -- ON ENTRY :
  --   p       a polynomial;
  --   x       a vector of range 1..p.dim.

  -- ON RETURN :
  --   y       the value of the monomial at x;
  --   yd      all partial derivatives of p at x;
  --   wrk     work space to hold intermediate values of the derivatives
  --           of the monomials.

  procedure Speel ( v : in Monomial_Vector;
                    x : in Vectors.Vector; powtab : in VecVecs.VecVec;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel ( v : in Monomial_Vector;
                    x : in Vectors.Vector; powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel ( v : in Link_to_Monomial_Vector;
                    x : in Vectors.Vector; powtab : in VecVecs.VecVec;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel ( v : in Link_to_Monomial_Vector;
                    x : in Vectors.Vector; powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the monomials m and
  --   all its derivatives at x, for general monomials with a power table.

  -- REQUIRED : v(i).nvr > 0, at least one variable has exponent 1.

  -- ON ENTRY :
  --   v        a monomial vector;
  --   x        a vector of range 1..v(i).dim;
  --   powtab   contains the powers of the values of the variables,
  --            needed for the higher degree powers in the monomial,
  --            used to compute the common factor b for each monomial in v.

  -- ON RETURN :
  --   y        the value of the monomial at x;
  --   yd       all partial derivatives of v at x;
  --   wrk      work space to hold intermediate values of the derivatives
  --            of the monomials.

  procedure Speel_without_Cache
              ( v : in Monomial_Vector; x : in Vectors.Vector;
                y : in out Ring.number;
                yd,wrk : in out Vectors.Vector );
  procedure Speel_without_Cache
              ( p : in Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel_without_Cache
              ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Similar to the general Speel, but without p.powtab.
  --   For debugging purposes only.

-- DESTRUCTORS :

  procedure Clear ( v : in out Monomial_Vector );
  procedure Clear ( v : in out Link_to_Monomial_Vector );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial vector v.

  procedure Clear ( p : in out Polynomial );
  procedure Clear ( p : in out Link_to_Polynomial );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the polynomial p.

end Generic_Monomial_Vectors;
