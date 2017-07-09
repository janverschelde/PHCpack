with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;           use DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;           use QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Homotopy_Membership_Tests is

-- DESCRIPTION :
--   This package offers an implementation of the homotopy test to
--   solve the membership problem for a point to be on a component.

  procedure Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );

  -- DESCRIPTION :
  --   Implementation of the homotopy membership test procedure,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by an ordinary polynomial system.

  -- REQUIRED : the sampling machine is properly initialized and tuned.

  -- ON ENTRY :
  --   verbose  if true, then diagnostic output is written to screen,
  --            if false, then the procedures remain silent;
  --   genpts   list of generic points in a witness set;
  --   x        the test point;
  --   adjsli   coefficients of hyperplane with adjusted constant so that x 
  --            satisfies the equations with coefficients in adjsli;
  --   homtol   tolerance for the homotopy membership test, used to compare
  --            x to the newly computed generic points on the adjusted slices.

  -- ON RETURN :
  --   found    true if the test point belongs to the witness set
  --            represented by the list of generic points, false otherwise.

  procedure Laurent_Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Laurent_Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Laurent_Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );

  -- DESCRIPTION :
  --   Implementation of the homotopy membership test procedure,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by a Laurent polynomial system.

  -- REQUIRED : the sampling machine is properly initialized and tuned.

  -- ON ENTRY :
  --   verbose  if true, then diagnostic output is written to screen,
  --            if false, then the procedures remain silent;
  --   genpts   list of generic points in a witness set;
  --   x        the test point;
  --   adjsli   coefficients of hyperplane with adjusted constant so that x 
  --            satisfies the equations with coefficients in adjsli;
  --   homtol   tolerance for the homotopy membership test, used to compare
  --            x to the newly computed generic points on the adjusted slices.

  -- ON RETURN :
  --   found    true if the test point belongs to the witness set
  --            represented by the list of generic points, false otherwise.

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );

  -- DESCRIPTION :
  --   Implementation of the homotopy membership test procedure,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by an ordinary polynomial system.
  --   We invoke the Sample_with_Stop routine from the Sampling Machine.
  --   Since this routine stops when the given vector x is found, 
  --   we compare x against the last vector in the list of new solutions.
  --   Intermediate output is written to file.

  procedure Laurent_Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Laurent_Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );
  procedure Laurent_Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean );

  -- DESCRIPTION :
  --   Implementation of the homotopy membership test procedure,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by a Laurent polynomial system.
  --   We invoke the Sample_with_Stop routine from the Sampling Machine.
  --   Since this routine stops when the given vector x is found, 
  --   we compare x against the last vector in the list of new solutions.
  --   Intermediate output is written to file.

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );

  -- DESCRIPTION :
  --   This procedure performs the component test for solution s,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by an ordinary polynomial system.

  -- REQUIRED : the sampling machine is initialized and tuned.

  -- ON ENTRY :
  --   verbose  if true, the diagnostic output is written to screen,
  --            if false, then the procedures remain silent;
  --   ep       embedded polynomial system;
  --   dim      dimension of solution components;
  --   sli      slices that cut out the generic points;
  --   genpts   generic points on the components;
  --   x        coordinates of the point for testing;
  --   restol   tolerance for residual;
  --   homtol   tolerance for homotopy test.

  -- ON RETURN :
  --   success  true when test point satisfies the equations,
  --            false otherwise;
  --   found    true when homotopy test succeeds in finding a matching
  --            generic point, found is false when not success.

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean );

  -- DESCRIPTION :
  --   This procedure performs the component test for solution s,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by a Laurent polynomial system.

  -- REQUIRED : the sampling machine is initialized and tuned.

  -- ON ENTRY :
  --   verbose  if true, the diagnostic output is written to screen,
  --            if false, then the procedures remain silent;
  --   ep       embedded polynomial system;
  --   dim      dimension of solution components;
  --   sli      slices that cut out the generic points;
  --   genpts   generic points on the components;
  --   x        coordinates of the point for testing;
  --   restol   tolerance for residual;
  --   homtol   tolerance for homotopy test.

  -- ON RETURN :
  --   success  true when test point satisfies the equations,
  --            false otherwise;
  --   found    true when homotopy test succeeds in finding a matching
  --            generic point, found is false when not success.

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                s : in Standard_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                s : in DoblDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                s : in QuadDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean );

  -- DESCRIPTION :
  --   This procedure performs the component test for solution s,
  --   for a witness set defined by an ordinary polynomial system.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   ep       embedded polynomial system;
  --   dim      dimension of solution components;
  --   sli      slices that cut out the generic points;
  --   genpts   generic points on the components;
  --   s        point up for testing;
  --   restol   tolerance for residual;
  --   homtol   tolerance for homotopy test.

  -- ON RETURN :
  --   success  true when test point satisfies the equations,
  --            false otherwise;
  --   found    true when homotopy test succeeds in finding a matching
  --            generic point, found is false when not success.

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                s : in Standard_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                s : in DoblDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                s : in QuadDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean );

  -- DESCRIPTION :
  --   This procedure performs the component test for solution s,
  --   for a witness set defined by a Laurent polynomial system.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   ep       embedded polynomial system;
  --   dim      dimension of solution components;
  --   sli      slices that cut out the generic points;
  --   genpts   generic points on the components;
  --   s        point up for testing;
  --   restol   tolerance for residual;
  --   homtol   tolerance for homotopy test.

  -- ON RETURN :
  --   success  true when test point satisfies the equations,
  --            false otherwise;
  --   found    true when homotopy test succeeds in finding a matching
  --            generic point, found is false when not success.

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );

  -- DESCRIPTION :
  --   Runs the homotopy membership test on a list of points in sols,
  --   in standard double, double double, or quad double precision,
  --   on a witness sets for an ordinary polynomial system.
  --   This procedure initializes and tunes the sampling machine.

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );

  -- DESCRIPTION :
  --   Runs the homotopy membership test on a list of points in sols,
  --   in standard double, double double, or quad double precision,
  --   on a witness sets for a Laurent polynomial system.
  --   This procedure initializes and tunes the sampling machine.

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );

  -- DESCRIPTION :
  --   Runs the homotopy membership test on a list of points in sols,
  --   in standard double, double double, or quad double precision,
  --   on a witness set defined by an ordinary polynomial systems.

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float );

  -- DESCRIPTION :
  --   Runs the homotopy membership test on a list of points in sols,
  --   in standard double, double double, or quad double precision,
  --   on a witness set defined by a Laurent polynomial systems.

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out Standard_Complex_Solutions.Solution_List );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   This procedure performs the component test for the solutions sols,
  --   for a witness set defined by an ordinary polynomial system.
  --   If provided, the solutions that are not found to be components are
  --   appended to the list isols.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   ep       embedded polynomial system;
  --   dim      number of slices;
  --   genpts   generic points on the slices;
  --   sols     points up for testing, do they lie on the components?
  --   restol   tolerance for residual;
  --   homtol   tolerance for homotopy test.

  -- ON RETURN :
  --   isols    solutions found to be isolated;
  --   isols_last points to the last element in isols.

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out Standard_Complex_Solutions.Solution_List );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   This procedure performs the component test for the solutions sols,
  --   for a witness set defined by a Laurent polynomial system.
  --   If provided, the solutions that are not found to be components are
  --   appended to the list isols.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   ep       embedded polynomial system;
  --   dim      number of slices;
  --   genpts   generic points on the slices;
  --   sols     points up for testing, do they lie on the components?
  --   restol   tolerance for residual;
  --   homtol   tolerance for homotopy test.

  -- ON RETURN :
  --   isols    solutions found to be isolated;
  --   isols_last points to the last element in isols.

end Homotopy_Membership_Tests;
