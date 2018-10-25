with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Jacobian_Rabinowitsch_Trick is

-- DESSCRIPTION :
--   Formulating a polynomial system along the lines of the trick of 
--   Rabinowitsch allows to move singular solutions to infinity.

  function Identity_Matrix
              ( n : integer32 ) return Standard_Complex_Matrices.Matrix;
  function Identity_Matrix
              ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix;
  function Identity_Matrix
              ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

  procedure Add_Last_Multiplier
              ( p : in out Standard_Complex_Polynomials.Poly;
                d : in integer32 );
  procedure Add_Last_Multiplier
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                d : in integer32 );
  procedure Add_Last_Multiplier
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                d : in integer32 );

  -- DESCRIPTION :
  --   Updates the last equation as p*y - 1,
  --   where y is the variable at index d.

  function Jacobian_Rabinowitsch
              ( p : Standard_Complex_Poly_Systems.Poly_Sys )
              return Standard_Complex_Poly_Systems.Poly_Sys;
  function Jacobian_Rabinowitsch
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
              return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Jacobian_Rabinowitsch
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
              return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the system augmented with the Jacobian matrix which
  --   puts singular solutions at infinity.

  function Jacobian_Rabinowitsch
              ( s : Standard_Complex_Solutions.Solution )
              return Standard_Complex_Solutions.Solution;
  function Jacobian_Rabinowitsch
              ( s : DoblDobl_Complex_Solutions.Solution )
              return DoblDobl_Complex_Solutions.Solution;
  function Jacobian_Rabinowitsch
              ( s : QuadDobl_Complex_Solutions.Solution )
              return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Adds as many zero coordinates to the solution as s.n
  --   and adds one extra coordinate 1.0 for the last variable.

  function Jacobian_Rabinowitsch
              ( s : Standard_Complex_Solutions.Solution_List )
              return Standard_Complex_Solutions.Solution_List;
  function Jacobian_Rabinowitsch
              ( s : DoblDobl_Complex_Solutions.Solution_List )
              return DoblDobl_Complex_Solutions.Solution_List;
  function Jacobian_Rabinowitsch
              ( s : QuadDobl_Complex_Solutions.Solution_List )
              return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Adds as many zero coordinates to each solution in s 
  --   as the number of coordinates in the solution,
  --   and adds one extra coordinate 1.0 for the last variable.

  procedure Add_Trick_Symbols ( nvar : in natural32 );

  -- DESCRIPTION :
  --   Adds nvar+1 symbols to the symbol table for the multiplier variables
  --   and the last variable used in the Rabinowitsch trick.

end Jacobian_Rabinowitsch_Trick;
