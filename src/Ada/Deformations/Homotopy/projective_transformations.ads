with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Projective_Transformations is

-- DESCRIPTION :
--   This package exports methods to perform projective transformations
--   on polynomials, polynomial systems, solutions, and lists of solutions,
--   in standard double, double double, and quad double precision.

  function Projective_Transformation 
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly;
  function Projective_Transformation
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  procedure Projective_Transformation
              ( p : in out Standard_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out DoblDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out QuadDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   An additional unknown is added so that each term has the same degree.

  function Projective_Transformation
             ( s : Standard_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : DoblDobl_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : QuadDobl_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Adds 1 as last component to the solution.

  function Projective_Transformation
             ( sols : Standard_Complex_Solutions.Solution_List )
           return Standard_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : DoblDobl_Complex_Solutions.Solution_List )
           return DoblDobl_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : QuadDobl_Complex_Solutions.Solution_List )
           return QuadDobl_Complex_Solutions.Solution_List;
  procedure Projective_Transformation
              ( sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   An additional component, equal to 1, is added to the solution vector.

  function Affine_Transformation
             ( s : Standard_Complex_Solutions.Solution )
           return Standard_Complex_Solutions.Solution;
  function Affine_Transformation
             ( s : DoblDobl_Complex_Solutions.Solution )
           return DoblDobl_Complex_Solutions.Solution;
  function Affine_Transformation
             ( s : QuadDobl_Complex_Solutions.Solution )
           return QuadDobl_Complex_Solutions.Solution;
  procedure Affine_Transformation
              ( sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Affine_Transformation
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Affine_Transformation
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   All components of the solution vector will be divided by the last
  --   component, which is afterwards cut off.

end Projective_Transformations;
