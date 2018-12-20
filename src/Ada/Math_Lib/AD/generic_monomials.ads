with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Abstract_Ring;
with Generic_Vectors;
with Generic_VecVecs;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);

package Generic_Monomials is

-- DESCRIPTION :
--   A monomial in several variables has a coefficient and exponents.
--   Higher exponents are stored separately for efficient differentiation.

  type Monomial is record
    cff : Ring.number; -- coefficient of the monomial
    nvr : natural32;   -- number of variables with positive exponent
    dim : natural32;   -- ambient dimension, total number of variables
    pos : Standard_Natural_Vectors.Link_to_Vector;
     -- array of range 1..nvr with positions of the variables
    exp : Standard_Natural_Vectors.Link_to_Vector;
     -- array of range 1..nvr with exponents of the variables
    n_base : natural32; -- number of variables with exponent larger than one
    pos_base : Standard_Natural_Vectors.Link_to_Vector;
     -- array of size 1..n_base with positions of variables
    exp_base : Standard_Natural_Vectors.Link_to_Vector;
     -- array of size 1..n_base with exponents of variables
    exp_tbl_base : Standard_Natural_Vectors.Link_to_Vector;
     -- stores the values of the exponents minus 2
  end record;

  type Link_to_Monomial is access Monomial;

-- CONSTRUCTORS :

  function Create ( c : Ring.number; dim : natural32 ) return Monomial;

  -- DESCRIPTION :
  --   Returns the constant monomial with coefficient c
  --   with ambient dimension dim. 

  function Create ( c : Ring.number;
                    e : Standard_Natural_Vectors.Vector ) return Monomial;
  function Create ( c : Ring.number;
                    e : Standard_Natural_Vectors.Vector )
                  return Link_to_Monomial;

  -- DESCRIPTION :
  --   Returns the constant monomial with coefficient c
  --   and exponents in e.

-- SELECTORS :

  function Degree ( m : Monomial ) return integer32;
  function Degree ( m : Link_to_Monomial ) return integer32;

  -- DESCRIPTION :
  --   Returns -1 for a constant or empty monomial,
  --   otherwise, returns the degree of the monomial.

  function Largest_Exponent ( m : Monomial ) return natural32;
  function Largest_Exponent ( m : Link_to_Monomial ) return natural32;

  -- DESCRIPTION :
  --   Returns the largest value of the exponent in m.
  --   If m is an empty monomial, then 0 is returned.

-- EVALUATION and DIFFERENTIATION :

  function Eval ( m : Monomial; x : Vectors.Vector ) return Ring.number;
  function Eval ( m : Link_to_Monomial;
                  x : Vectors.Vector ) return Ring.number;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate m at x.

  procedure Diff ( m : in Monomial; x : in Vectors.Vector;
                   yd : in out Vectors.Vector );
  procedure Diff ( m : in Link_to_Monomial; x : in Vectors.Vector;
                   yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to differentiate m at x.
  --   The evaluated derivatives of m at x are in the vector yd,
  --   stored in the first m.nvr positions.

  procedure Speel ( m : in Monomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( m : in Link_to_Monomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the monomial m and
  --   all its derivatives at x, for monomials that are variable products.
  --
  -- NOTE that this procedure allocates, computes, and then deallocates
  --   the power table for a one time use.  For evaluating many monomials
  --   at the same vector x, the power table should be given on input and
  --   shared between all values, as done by the last Speel procedure.

  -- REQUIRED : m.nvr > 0, at least one variable has exponent 1.

  -- ON ENTRY :
  --   m        a monomial;
  --   x        a vector of range 1..m.dim.

  -- ON RETURN :
  --   y        the value of the monomial at x;
  --   yd       all partial derivatives of m at x,
  --            stored in the first m.nvr positions.

  procedure Speel_on_Product
              ( m : in Monomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel_on_Product
              ( m : in Link_to_Monomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the monomial m and
  --   all its derivatives at x, for monomials that are variable products.

  -- REQUIRED : m.nvr > 0, at least one variable has exponent 1 and
  --   m.n_base = 0, i.e.: no variables appear with exponent 2 or higher.

  -- ON ENTRY :
  --   m        a monomial;
  --   x        a vector of range 1..m.dim.

  -- ON RETURN :
  --   y        the value of the monomial at x;
  --   yd       all partial derivatives of m at x,
  --            stored in the first m.nvr positions.

  procedure Power_Table
              ( m : in monomial; x : Vectors.Vector;
                powtab : out VecVecs.VecVec );
  procedure Power_Table
              ( m : in Link_to_monomial; x : Vectors.Vector;
                powtab : out VecVecs.VecVec );

  -- DESCRIPTION :
  --   Allocates memory for the rows in the power table and computes the
  --   highest powers of the values of x according to the degrees in m.
  --   On entry, powtab has as many rows as m.dim (same range as x'range)
  --   and as powtab has many columns as the corresponding highest power
  --   of the variable at the row position.
  --   The second, column index, starts at 0, with the value of x(row)
  --   at the start position 0.

  procedure Common_Factor
              ( m : in Monomial; powtab : in VecVecs.VecVec;
                y : in out Ring.number );
  procedure Common_Factor
              ( m : in Monomial; powtab : in VecVecs.Link_to_VecVec;
                y : in out Ring.number );
  procedure Common_Factor
              ( m : in Link_to_Monomial; powtab : in VecVecs.VecVec;
                y : in out Ring.number );
  procedure Common_Factor
              ( m : in Link_to_Monomial; powtab : in VecVecs.Link_to_VecVec;
                y : in out Ring.number );

  -- DESCRIPTION :
  --   Given the powers of the values of the variables which occur with
  --   higher power in the monomial m, returns the product of the
  --   coefficient of m with the product of all the evaluated powers.
  --   This value is the common factor b in the general Speel procedure.

  -- ON ENTRY :
  --    m         a monomial in several variables;
  --    powtab    powtab[i][e] contains x(i) to the power e,
  --              where e is defined in m.exp_tbl_base.

  -- ON RETURN :
  --    y         product of the coefficient of m with the powers,
  --              corresponding to m.pos_base and m.exp_tbl_base.

  procedure Speel ( m : in Monomial;
                    x : in Vectors.Vector; b : in Ring.number;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( m : in Link_to_Monomial;
                    x : in Vectors.Vector; b : in Ring.number;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the monomial m and
  --   all its derivatives at x, for general monomials with evaluated base.

  -- REQUIRED : m.nvr > 0, at least one variable has exponent 1.

  -- ON ENTRY :
  --   m        a monomial;
  --   x        a vector of range 1..m.dim;
  --   b        is the common factor, defined by all higher degree powers
  --            of the variables in the monomial, evaluated at x,
  --            and multiplied with the coefficient of the monomial.

  -- ON RETURN :
  --   y        the value of the monomial at x;
  --   yd       all partial derivatives of m at x,
  --            stored in the first m.nvr positions.

  procedure Speel ( m : in Monomial;
                    x : in Vectors.Vector; powtab : in VecVecs.VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( m : in Monomial;
                    x : in Vectors.Vector; powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( m : in Link_to_Monomial;
                    x : in Vectors.Vector; powtab : in VecVecs.VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( m : in Link_to_Monomial;
                    x : in Vectors.Vector; powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the monomial m and
  --   all its derivatives at x, for general monomials with a power table.

  -- REQUIRED : m.nvr > 0, at least one variable has exponent 1.

  -- ON ENTRY :
  --   m        a monomial;
  --   x        a vector of range 1..m.dim;
  --   powtab   contains the powers of the values of the variables,
  --            needed for the higher degree powers in the monomial,
  --            used to compute the common factor b.

  -- ON RETURN :
  --   y        the value of the monomial at x;
  --   yd       all partial derivatives of m at x,
  --            stored in the first m.nvr positions.

-- DESTRUCTORS :

  procedure Clear ( m : in out Monomial );
  procedure Clear ( m : in out Link_to_Monomial );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial m.

end Generic_Monomials;
