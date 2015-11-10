with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Extrinsic_Diagonal_Solvers is

-- DESCRIPTION :
--   This package offers drivers to compute witness sets on the intersection
--   of positive dimensional components by extrinsic diagonal homotopies.
--   Note that the routines in this package do not solve, but provide the
--   homotopies on which standard continuation can be applied directly.

  procedure Randomize_System;

  -- DESCRIPTION :
  --   Interactive driver to randomize a system so that the number
  --   of equations equals the co-dimension of the solution component
  --   which is to be sampled.

  procedure Build_Diagonal_Cascade;

  -- DESCRIPTION :
  --   Interactive driver to build a homotopy to start a cascade
  --   of homotopies to compute all positive dimensional components
  --   of the intersection of two components.

  procedure Collapse_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim,add2dim : in natural32;
                r : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Collapse_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim,add2dim : in natural32;
                r : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Collapse_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim,add2dim : in natural32;
                r : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Removes the duplicate variables from p and sols,
  --   in standard double, double double, or quad double precision,
  --   embedding the result as a component of dimension dim+addtodim.

  -- ON ENTRY :
  --   p        an embedded polynomial system, as many slacks as dim;
  --   sols     solutions to p, with the proper embedding;
  --   dim      number of slack variables currently in the embedding;
  --   add2dim  number of slack variables that must be added.

  -- ON RETURN :
  --   sols     solutions with dim + add2dim slack variables;
  --   r        system with the diagonal removed, embedded with
  --            as many as dim + add2dim slack variables.

  procedure Collapse_Diagonal_System;

  -- DESCRIPTION :
  --   Once the target system in the diagonal homotopy has been solved,
  --   we want to return to the original set of coordinates.
  --   This is what this interactive driver does.

end Extrinsic_Diagonal_Solvers;
