with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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

package Homotopy_Membership_Filters is

-- DESCRIPTION :
--   Given a witness set for a higher dimensional solution set
--   and a list of solutions, a homotopy membership filter splits
--   the given list of solutions in two lists:
--   1) the solutions which belong to the solution set; and
--   2) the solutions which do not belong to the solution set.
--   Elements in the list which do not satisfy the given polynomials
--   are remove and in neither of the two lists in the output of the
--   homotopy membership filter.

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Runs a homotopy membership filter on the solutions in totest,
  --   in double, double double, or quad double precision,
  --   for ordinary or Laurent polynomial systems.
  --   No output is written to file.

  -- ON ENTRY :
  --   verbose  if true then output is written to screen,
  --            otherwise the procedure remains silent;
  --   eqs      the equations of the witness set;
  --   pts      generic points in the witness set;
  --   dim      dimension of the witness set;
  --   rcotol   tolerance on the inverse of the condition number estimate,
  --            if the rco of a solution is less than rcotol,
  --            then the solution is considered regular and therefore
  --            it cannot lie on the set represented by the witness set
  --            and this bypasses the homotopy membership test,
  --            if rcotol is set to 0.0, then all points will
  --            go through the homotopy membership test;
  --   restol   tolerance on the residual;
  --   homtol   tolerance for the homotopy membership test;
  --   totest   solutions to test whether or not they belong to the
  --            solution set represented by the witness set.

  -- ON RETURN :
  --   mempts   points in totest which are a member of the solution set,
  --            represented by the witness set;
  --   outpts   outside points in totest not a member of the solution set.

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Runs a homotopy membership filter on the solutions in totest,
  --   in double, double double, or quad double precision,
  --   for ordinary or Laurent polynomial systems.
  --   No output is written to file.

  -- ON ENTRY :
  --   verbose  if true then output is written to screen,
  --            otherwise the procedure remains silent;
  --   nt       number of tasks to be used in the test;
  --   eqs      the equations of the witness set;
  --   pts      generic points in the witness set;
  --   dim      dimension of the witness set;
  --   rcotol   tolerance on the inverse of the condition number estimate,
  --            if the rco of a solution is less than rcotol,
  --            then the solution is considered regular and therefore
  --            it cannot lie on the set represented by the witness set
  --            and this bypasses the homotopy membership test,
  --            if rcotol is set to 0.0, then all points will
  --            go through the homotopy membership test;
  --   restol   tolerance on the residual;
  --   homtol   tolerance for the homotopy membership test;
  --   totest   solutions to test whether or not they belong to the
  --            solution set represented by the witness set.

  -- ON RETURN :
  --   mempts   points in totest which are a member of the solution set,
  --            represented by the witness set;
  --   outpts   outside points in totest not a member of the solution set.

  procedure Filter
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Runs a homotopy membership filter on the solutions in totest,
  --   in double, double double, or quad double precision,
  --   for ordinary or Laurent polynomial systems.
  --   Output is written to file.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   eqs      the equations of the witness set;
  --   pts      generic points in the witness set;
  --   dim      dimension of the witness set;
  --   rcotol   tolerance on the inverse of the condition number estimate,
  --            if the rco of a solution is less than rcotol,
  --            then the solution is considered regular and therefore
  --            it cannot lie on the set represented by the witness set
  --            and this bypasses the homotopy membership test,
  --            if rcotol is set to 0.0, then all points will
  --            go through the homotopy membership test;
  --   restol   tolerance on the residual;
  --   homtol   tolerance for the homotopy membership test;
  --   totest   solutions to test whether or not they belong to the
  --            solution set represented by the witness set.

  -- ON RETURN :
  --   mempts   points in totest which are a member of the solution set,
  --            represented by the witness set;
  --   outpts   outside points in totest not a member of the solution set.

  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Runs a homotopy membership filter on the solutions in totest,
  --   in double, double double, or quad double precision,
  --   for ordinary or Laurent polynomial systems.
  --   Output is written to file.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nt       number of tasks;
  --   eqs      the equations of the witness set;
  --   pts      generic points in the witness set;
  --   dim      dimension of the witness set;
  --   rcotol   tolerance on the inverse of the condition number estimate,
  --            if the rco of a solution is less than rcotol,
  --            then the solution is considered regular and therefore
  --            it cannot lie on the set represented by the witness set
  --            and this bypasses the homotopy membership test,
  --            if rcotol is set to 0.0, then all points will
  --            go through the homotopy membership test;
  --   restol   tolerance on the residual;
  --   homtol   tolerance for the homotopy membership test;
  --   totest   solutions to test whether or not they belong to the
  --            solution set represented by the witness set.

  -- ON RETURN :
  --   mempts   points in totest which are a member of the solution set,
  --            represented by the witness set;
  --   outpts   outside points in totest not a member of the solution set.

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );

  -- DESCRIPTION :
  --   Junk points are removed from the given witness supersets,
  --   running successive homotopy membership filters,
  --   in double, double double, or quad double precision,
  --   for witness sets defined by ordinary or Laurent systems.

  -- ON ENTRY :
  --   verbose  if true then output is written to screen,
  --            otherwise the procedure remains silent;
  --   eqs      sequence of embedded systems, with range 0..topdim;
  --   pts      candidate witness points, for the range 0..topdim,
  --            no junk is assumed in pts(topdim);
  --   topdim   the top dimension of the solution set;
  --   rcotol   tolerance on the inverse of the estimated condition number,
  --            which allows to bypass the homotopy membership if > 0.0;
  --   restol   tolerance on the residual;
  --   homtol   tolerance for the homotopy membership test.

  -- ON RETURN :
  --   pts      points on higher dimensional sets are removed.
  --   filcnt   counts the number of solutions left at each stage
  --            after filtering the junk points, array of range 0..topdim;
  --   times    elapsed CPU user time at each stage;
  --   alltime  the total elapsed CPU user time.

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );

  -- DESCRIPTION :
  --   Junk points are removed from the given witness supersets,
  --   running successive homotopy membership filters,
  --   in double, double double, or quad double precision,
  --   for witness sets defined by ordinary or Laurent systems.

  -- ON ENTRY :
  --   verbose  if true then output is written to screen,
  --            otherwise the procedure remains silent;
  --   nt       number of tasks;
  --   eqs      sequence of embedded systems, with range 0..topdim;
  --   pts      candidate witness points, for the range 0..topdim,
  --            no junk is assumed in pts(topdim);
  --   topdim   the top dimension of the solution set;
  --   rcotol   tolerance on the inverse of the estimated condition number,
  --            which allows to bypass the homotopy membership if > 0.0;
  --   restol   tolerance on the residual;
  --   homtol   tolerance for the homotopy membership test.

  -- ON RETURN :
  --   pts      points on higher dimensional sets are removed.
  --   filcnt   counts the number of solutions left at each stage
  --            after filtering the junk points, array of range 0..topdim;
  --   times    elapsed CPU user time at each stage;
  --   alltime  the total elapsed CPU user time.

end Homotopy_Membership_Filters;
