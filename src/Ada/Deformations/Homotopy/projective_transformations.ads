with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with PentDobl_Complex_Solutions;
with OctoDobl_Complex_Solutions;
with DecaDobl_Complex_Solutions;
with HexaDobl_Complex_Solutions;

package Projective_Transformations is

-- DESCRIPTION :
--   This package exports methods to perform projective transformations
--   on polynomials, polynomial systems, solutions, and lists of solutions,
--   in double, double double, triple double, quad double, penta double,
--   octo double, deca double, and hexa double precision.

  function Projective_Transformation 
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : TripDobl_Complex_Polynomials.Poly )
             return TripDobl_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : PentDobl_Complex_Polynomials.Poly )
             return PentDobl_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : OctoDobl_Complex_Polynomials.Poly )
             return OctoDobl_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : DecaDobl_Complex_Polynomials.Poly )
             return DecaDobl_Complex_Polynomials.Poly;
  function Projective_Transformation 
             ( p : HexaDobl_Complex_Polynomials.Poly )
             return HexaDobl_Complex_Polynomials.Poly;
  function Projective_Transformation
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys )
             return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys )
             return PentDobl_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys;
  function Projective_Transformation
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys;

  procedure Projective_Transformation
              ( p : in out Standard_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out DoblDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out TripDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out TripDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out QuadDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out PentDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out PentDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out OctoDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out OctoDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out DecaDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out DecaDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Projective_Transformation
              ( p : in out HexaDobl_Complex_Polynomials.Poly );
  procedure Projective_Transformation
              ( p : in out HexaDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   An additional unknown is added so that each term has the same degree.

  function Projective_Transformation
             ( s : Standard_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : DoblDobl_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : TripDobl_Complex_Solutions.Solution )
             return TripDobl_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : QuadDobl_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : PentDobl_Complex_Solutions.Solution )
             return PentDobl_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : OctoDobl_Complex_Solutions.Solution )
             return OctoDobl_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : DecaDobl_Complex_Solutions.Solution )
             return DecaDobl_Complex_Solutions.Solution;
  function Projective_Transformation
             ( s : HexaDobl_Complex_Solutions.Solution )
             return HexaDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Adds 1 as last component to the solution.

  function Projective_Transformation
             ( sols : Standard_Complex_Solutions.Solution_List )
           return Standard_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : DoblDobl_Complex_Solutions.Solution_List )
           return DoblDobl_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : TripDobl_Complex_Solutions.Solution_List )
           return TripDobl_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : QuadDobl_Complex_Solutions.Solution_List )
           return QuadDobl_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : PentDobl_Complex_Solutions.Solution_List )
           return PentDobl_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : OctoDobl_Complex_Solutions.Solution_List )
           return OctoDobl_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : DecaDobl_Complex_Solutions.Solution_List )
           return DecaDobl_Complex_Solutions.Solution_List;
  function Projective_Transformation
             ( sols : HexaDobl_Complex_Solutions.Solution_List )
           return HexaDobl_Complex_Solutions.Solution_List;
  procedure Projective_Transformation
              ( sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out TripDobl_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out PentDobl_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out OctoDobl_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out DecaDobl_Complex_Solutions.Solution_List );
  procedure Projective_Transformation
              ( sols : in out HexaDobl_Complex_Solutions.Solution_List );

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
