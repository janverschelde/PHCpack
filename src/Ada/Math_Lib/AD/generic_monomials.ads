with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Abstract_Ring;
with Generic_Vectors;
with Generic_Matrices;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);

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

  -- DESCRIPTION :
  --   Returns the constant monomial with coefficient c
  --   and exponents in e.

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
  --   all its derivatives at x.

  -- REQUIRED : m.nvr > 0, at least one variable has exponent 1 and
  --   m.n_base = 0, i.e.: no variables appear with exponent 2 or higher.

  -- ON ENTRY :
  --   m       a monomial;
  --   x       a vector of range 1..m.dim.

  -- ON RETURN :
  --   y       the value of the monomial at x;
  --   yd      all partial derivatives of m at x,
  --           stored in the first m.nvr positions.

-- DESTRUCTORS :

  procedure Clear ( m : in out Monomial );
  procedure Clear ( m : in out Link_to_Monomial );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial m.

end Generic_Monomials;
