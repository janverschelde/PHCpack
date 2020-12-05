with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Main_Factorization is

-- DESCRIPTION :
--   The f in factorization applies also to filtering junk points
--   with a homotopy membership test.

  procedure Tune_Member_Tolerances ( restol,homtol : in out double_float );

  -- DESCRIPTION :
  --   Displays the current values for the tolerances to decide whether
  --   a residual is small enough for the point to lie on a hypersurface,
  --   and for a test point to belong to a witness set.

  procedure Standard_Multitasked_Membership_Test
              ( file : in file_type;
                nbtasks,dim : in natural32; homtol : in double_float;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                gpts,sols : in Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   When nt > 0, then the multitasking membership test runs
  --   in double precision.

  procedure DoblDobl_Multitasked_Membership_Test
              ( file : in file_type;
                nbtasks,dim : in natural32; homtol : in double_float;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                gpts,sols : in DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   When nt > 0, then the multitasking membership test runs
  --   in double double precision.

  procedure QuadDobl_Multitasked_Membership_Test
              ( file : in file_type;
                nbtasks,dim : in natural32; homtol : in double_float;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                gpts,sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   When nt > 0, then the multitasking membership test runs
  --   in quad double precision.

  procedure Standard_Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 );
  procedure DoblDobl_Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 );
  procedure QuadDobl_Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   in double, double double, or quad double precision arithmetic.
  --   The number of tasks is in nt and 
  --   vrblvl is the value of the verbose level.

  procedure Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the precision and then runs the membership test
  --   in the selected precision.  The number of tasks is in nt
  --   and vrblvl is the value of the verbose level.

  procedure Trace_Form_Interpolation ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies the linear trace to certify the monodromy breakup.
  --   The value of the verbose level in in vrblvl.

  procedure Incremental_Interpolation ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Interpolates solution sets with multivariate polynomials of
  --   increasing degrees to determine the irreducible components.
  --   The value of the verbose level in in vrblvl.

  procedure Standard_Eliminate
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim : in integer32 );

  -- DESCRIPTION :
  --   Given the embedded system with its solutions, this procedure
  --   initializes the sampler and calls the interpolation routines.
  --   All calculations are done with standard floating-point arithmetic.

  procedure Multprec_Eliminate
                ( file : in file_type;
                  ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim : in integer32; size : in natural32 );

  -- DESCRIPTION :
  --   Given the embedded system with its solutions,
  --   initializes the sampler and calls the interpolation routines.
  --   Sample refinement and interpolation is done with multiprecision
  --   arithmetic.

  procedure Eliminate_Variables ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given an embedding with slices parallel to a given subspace,
  --   with interpolation we eliminate variables not belonging to this
  --   subspace from the system.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -f.  Its main function is to factor a pure dimensional
  --   solution set into irreducible components.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used wherever defined;
  --   infilename     file name for input (embedded system+generic points),
  --                  if empty, then the user will be prompted to supply
  --                  the name of an input file;
  --   outfilename    name of file for output, if empty, then the user will
  --                  be asked to supply the name of an output file;
  --   vrblvl         is the verbose level.

end Main_Factorization;
