with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;

package Drivers_to_Breakup_Components is

-- DESCRIPTION :
--   This package offers driver routines to perform a breakup of
--   components of solutions.

  procedure Standard_Breakup_Components
               ( file : in file_type; level,itp : in natural32;
                 skewproj : in boolean;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List );

  procedure Multprec_Breakup_Components
               ( file : in file_type; level,size,itp : in natural32;
                 skewproj : in boolean;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List );

  -- ON ENTRY :
  --   file      to write results and diagnostics on;
  --   level     number of slices added;
  --   size      size of the numbers;
  --   itp       interpolation type :
  --              = 1 : massive interpolation with full grid of points,
  --              = 2 : incremental interpolation, one point after the other,
  --              = 3 : subspace restriction and projection from point;
  --   p         embedded polynomial system;
  --   mp        multi-precision version of original polynomial system;
  --   sols      generic points.

  procedure Breakup_Menu ( itp,size : out natural32; skewproj : out boolean );

  -- DESCRIPTION :
  --   Displays the menu and lets the user select the interpolator type.
  --   If no multi-precision is selected, then size = 0 on return,
  --   otherwise size indicates the size of the numbers to be used.
  --   The skewproj is true on return when skew line projections are used.

  procedure Breakup_with_Interpolation_Filters;

  -- DESCRIPTION :
  --   This is the interactive driver to breakup solution sets into
  --   irreducible components by means of interpolation filters.

  procedure Monodromy_Decomposition
                   ( file : in file_type;
                     p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in Standard_Complex_Solutions.Solution_List;
                     dim : in natural32 );

  procedure Breakup_with_Monodromy_Group_Actions;

  -- DESCRIPTION :
  --   This interactive driver uses the monodromy group actions to
  --   predict the decomposition of an equidimensional solution set
  --   into irreducible components.

end Drivers_to_Breakup_Components;
