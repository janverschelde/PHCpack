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

package Cascade_Membership_Filters is

-- DESCRIPTION :
--   Given the witness set for the top dimensional set and a sequence 
--   of super witness sets for the lower dimensional solution sets,
--   a cascade membership filter removes the junk points from the
--   super witness sets and then outputs witness sets for the pure
--   dimensional solution sets in a numerical irreducible decomposition.
--
--   Two types of polynomial systems are supported :
--   (1) ordinary polynomial systems with nonnegative exponents,
--   (2) Laurent polynomial systems may have negative exponents.
--   The filters have three different kinds of output :
--   (1) optional output is written to screen,
--   (2) output to file and optional output to screen,
--   (3) output to callback procedures and optional output to screen.
--   Three different types of precision are supported :
--   (1) standard double precision is the hardware precision,
--   (2) double double precision doubles the hardware precision,
--   (3) quad double precision quadrules the hardware precision.
--   Multitasking is supported :
--   (1) if no number of tasks is provided, then no multitasking is applied,
--   (2) multitasking is applied with a positive number of tasks.
--   Combining all 2x3x3x2 options results in 36 Filter procedures.

-- SINGLE TASKED FILTERS WITH OPTIONAL OUTPUT TO SCREEN :

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

-- SINGLE TASKED FILTERS WITH CALLBACK AND OPTIONAL OUTPUT TO SCREEN :

  procedure Filter_with_Callback
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );

  -- DESCRIPTION :
  --   Junk points are removed from the given witness supersets,
  --   running successive homotopy membership filters,
  --   in double, double double, or quad double precision,
  --   for witness sets defined by ordinary or Laurent systems.
  --   Report_Witness_Set is called for each new filtered set.

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

-- SINGLE TASKED FILTERS WITH OUTPUT TO FILE :

  procedure Filter
              ( file : in file_type; verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean;
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
  --   verbose  if true then also the homotopy membership tests write
  --            to file, otherwise the homotopy membership tests are silent;
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

-- MULTITASKED FILTERS :

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

-- MULTITASKED FILTERS WITH CALLBACK AND OPTIONAL OUTPUT TO SCREEN :

  procedure Filter_with_Callback
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Filter_with_Callback
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );

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

-- MULTITASKED FILTERS WITH OUTPUT TO FILE :

  procedure Filter
              ( file : in file_type; verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Filter
              ( file : in file_type; verbose : in boolean; nt : in natural32;
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
  --   file     for output of diagnostics and results;
  --   verbose  if true then the homotopy membership tests write to file
  --            otherwise the homotopy membership tests remain silent;
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

end Cascade_Membership_Filters;
