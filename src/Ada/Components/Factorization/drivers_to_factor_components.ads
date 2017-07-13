with text_io;                               use text_io;
with Standard_Natural_Numbers;              use Standard_Natural_Numbers;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Drivers_to_Factor_Components is

-- DESCRIPTION :
--   This package collects drivers to different interpolators.

  procedure Call_Standard_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; itp : in natural32 );

  procedure Call_Multprec_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; itp : in natural32 );

  -- DESCRIPTION :
  --   The incremental interpolation scheme is a flexible way to determine
  --   the decomposition of a pure dimensional solution set into irreducible
  --   components of solutions.  However, the use standard arithmetic limits
  --   the application range.

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   The solutions given to this routine are witness points,
  --   defined by an ordinary polynomial system,
  --   computed in standard double precision.
  --   With monodromy we partition the set of witness points according
  --   to the irreducible components of the system p at dimension dim.
  --   The factorization is represented by f on return.

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   The solutions given to this routine are witness points,
  --   defined by a Laurent polynomial system,
  --   computed in standard double precision.
  --   With monodromy we partition the set of witness points according
  --   to the irreducible components of the system p at dimension dim.
  --   The factorization is represented by f on return.

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   The solutions given to this routine are witness points,
  --   defined by an ordinary polynomial system,
  --   computed in double double precision.
  --   With monodromy we partition the set of witness points according
  --   to the irreducible components of the system p at dimension dim.
  --   The factorization is represented by f on return.

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   The solutions given to this routine are witness points,
  --   defined by a Laurent polynomial system,
  --   computed in double double precision.
  --   With monodromy we partition the set of witness points according
  --   to the irreducible components of the system p at dimension dim.
  --   The factorization is represented by f on return.

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   The solutions given to this routine are witness points,
  --   defined by an ordinary polynomial system,
  --   computed in quad double precision.
  --   With monodromy we partition the set of witness points according
  --   to the irreducible components of the system p at dimension dim.
  --   The factorization is represented by f on return.

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   The solutions given to this routine are witness points,
  --   defined by a Laurent polynomial system,
  --   computed in quad double precision.
  --   With monodromy we partition the set of witness points according
  --   to the irreducible components of the system p at dimension dim.
  --   The factorization is represented by f on return.

  procedure Call_Newton_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   The application of Newton interpolation requires the samples to
  --   lie on a structured grid.  Therefore we need a decomposition of
  --   the set of generic points.

  procedure Call_Trace_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   The use of traces leads to another form of interpolating polynomials
  --   through generic points of a solution component.  Also here we need
  --   a predicted decomposition to set up a structured grid of samples.

  procedure Call_Power_Trace_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   The name "power traces" refers to the Newton identities for
  --   transforming values of the elementary symmetric functions into power
  --   sums.  With these identities we need fewer samples to set up the
  --   trace form of the interpolating polynomial.

  procedure Call_Linear_Trace_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   The major advantage of working with traces is that linear ones
  --   suffice to perform the validation process.

end Drivers_to_Factor_Components;
