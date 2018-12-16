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

  function Create ( c : Ring.number; dim : natural32 ) return Monomial;

  -- DESCRIPTION :
  --   Returns the constant monomial with coefficient c
  --   with ambient dimension dim. 

  function Create ( c : Ring.number;
                    e : Standard_Natural_Vectors.Vector ) return Monomial;

  -- DESCRIPTION :
  --   Returns the constant monomial with coefficient c
  --   and exponents in e.

  procedure Clear ( m : in out Monomial );
  procedure Clear ( m : in out Link_to_Monomial );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial m.

end Generic_Monomials;
