with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Multitasking_Membership_Tests is

-- DESCRIPTION :
--   A homotopy membership test decides whether a point belongs to
--   the solution set represented by a witness set.
--   This test may involve the tracking of as many paths as the degree
--   of the set.  With multitasking, the test can finish faster.

  function Is_Member ( s : Standard_Complex_Solutions.Solution_List;
                       z : Standard_Complex_Vectors.Vector;
                       tol : double_float ) return natural32;
  function Is_Member ( s : DoblDobl_Complex_Solutions.Solution_List;
                       z : DoblDobl_Complex_Vectors.Vector;
                       tol : double_float ) return natural32;
  function Is_Member ( s : QuadDobl_Complex_Solutions.Solution_List;
                       z : QuadDobl_Complex_Vectors.Vector;
                       tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   Returns 0 if z is not in s with respect to the tolerance tol,
  --   otherwise returns the index of the solution in s which matches z
  --   within the given tolerance.

  function Standard_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Poly_Systems.Poly_Sys;
               gpts : Standard_Complex_Solutions.Solution_List;
               tpnt : Standard_Complex_Vectors.Vector )
             return natural32;
  function Standard_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Laur_Systems.Laur_Sys;
               gpts : Standard_Complex_Solutions.Solution_List;
               tpnt : Standard_Complex_Vectors.Vector )
             return natural32;
  function DoblDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : DoblDobl_Complex_Solutions.Solution_List;
               tpnt : DoblDobl_Complex_Vectors.Vector )
             return natural32;
  function DoblDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : DoblDobl_Complex_Solutions.Solution_List;
               tpnt : DoblDobl_Complex_Vectors.Vector )
             return natural32;
  function QuadDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : QuadDobl_Complex_Solutions.Solution_List;
               tpnt : QuadDobl_Complex_Vectors.Vector )
             return natural32;
  function QuadDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : QuadDobl_Complex_Solutions.Solution_List;
               tpnt : QuadDobl_Complex_Vectors.Vector )
             return natural32;

  -- DESCRIPTION :
  --   Runs a multitasked homotopy membership test 
  --   in double, double double, or quad double precision,
  --   to determine if the point tptn belongs to the set of dimension dim,
  --   represented by start and gpts,
  --   for ordinary and Laurent polynomial systems.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   dim      dimension of the solution set;
  --   tol      tolerance used to compare coordinates;
  --   start    embedded system as start system for the test;
  --   gpts     generic points as start solutions;
  --   tpnt     point to test.

  procedure Standard_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in Standard_Complex_Poly_Systems.Poly_Sys;
               gpts : in Standard_Complex_Solutions.Solution_List;
               tpnt : in Standard_Complex_Vectors.Vector;
               idx : out natural32;
               match : out Standard_Complex_Vectors.Vector );
  procedure Standard_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in Standard_Complex_Laur_Systems.Laur_Sys;
               gpts : in Standard_Complex_Solutions.Solution_List;
               tpnt : in Standard_Complex_Vectors.Vector;
               idx : out natural32;
               match : out Standard_Complex_Vectors.Vector );
  procedure DoblDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : in DoblDobl_Complex_Solutions.Solution_List;
               tpnt : in DoblDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out DoblDobl_Complex_Vectors.Vector );
  procedure DoblDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : in DoblDobl_Complex_Solutions.Solution_List;
               tpnt : in DoblDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out DoblDobl_Complex_Vectors.Vector );
  procedure QuadDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : in QuadDobl_Complex_Solutions.Solution_List;
               tpnt : in QuadDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out QuadDobl_Complex_Vectors.Vector );
  procedure QuadDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : in QuadDobl_Complex_Solutions.Solution_List;
               tpnt : in QuadDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs a multitasked homotopy membership test
  --   in double, double double, or quad double precision,
  --   to determine if the point tptn belongs to the set 
  --   of dimension dim, represented by start and gpts,
  --   for ordinary and Laurent polynomial systems.
  --   Returns also the matching coordinates for testing purposes.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   dim      dimension of the solution set;
  --   tol      tolerance used to compare coordinates;
  --   start    embedded system as start system for the test;
  --   gpts     generic points as start solutions;
  --   tpnt     point to test.

  -- ON RETURN :
  --   idx      index of the matching point in the new list
  --            of generic points;
  --   match    coordinates of the matching point if idx /= 0.

-- FILTERS :

  function Standard_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Poly_Systems.Poly_Sys;
               gpts,sols : Standard_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List;
  function Standard_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Laur_Systems.Laur_Sys;
               gpts,sols : Standard_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List;
  function DoblDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               gpts,sols : DoblDobl_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List;
  function DoblDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               gpts,sols : DoblDobl_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List;
  function QuadDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               gpts,sols : QuadDobl_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List;
  function QuadDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               gpts,sols : QuadDobl_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns those points in sols which fail the membership test.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   dim      dimension of the solution set;
  --   tol      tolerance used to compare coordinates;
  --   start    embedded system as start system for the test;
  --   gpts     generic points as start solutions;
  --   sols     points to test.

-- WITH PREPROCESSING EVALUATION TEST :

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean; nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean; nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean; nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean; nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean; nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean; nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );

  -- DESCRIPTION :
  --   The test point is first evaluated in the system before
  --   running the membership test, as a preprocessing stage.
  --   This preprocessing saves time if the test point is not
  --   a junk point in a cascade and is user given.

  -- ON ENTRY :
  --   verbose  if true, the diagnostic output is written to screen,
  --            if false, then the procedure remain silent;
  --   nt       number of tasks;
  --   ep       embedded polynomial system;
  --   dim      dimension of solution components;
  --   genpts   generic points on the components;
  --   x        coordinates of the point for testing;
  --   restol   tolerance for residual;
  --   homtol   tolerance for homotopy test.

  -- ON RETURN :
  --   success  true when test point satisfies the equations,
  --            false otherwise;
  --   found    true when homotopy test succeeds in finding a matching
  --            generic point, found is false when not success.

end Multitasking_Membership_Tests;
