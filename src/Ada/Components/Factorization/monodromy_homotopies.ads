with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
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

package Monodromy_Homotopies is

-- DESCRIPTION :
--   Given a witness set for a pure dimensional solution set,
--   a monodromy homotopy computes a partition of the witness points
--   so each set in the partition corresponds to one irreducible factor
--   of the solution set.

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Partitions the witness points according to the irreducible factors,
  --   for ordinary polynomial systems and for Laurent systems,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   verbose  if true then output is written to screen,
  --            otherwise the procedure remains silent;
  --   eqs      the equations of the witness set;
  --   pts      generic points in the witness set;
  --   dim      dimension of the witness set;
  --   nbl      limit on the number of iterations which keep the
  --            factorization unchanged;
  --   tol      tolerance to decide whether points are equal.

  -- ON RETURN :
  --   f        irreducible decomposition of a witness set.

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration );

  -- DESCRIPTION :
  --   Partitions the witness points according to the irreducible factors,
  --   for ordinary polynomial systems and for Laurent systems,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   verbose  if true then output is written to screen,
  --            otherwise the procedure remains silent;
  --   eqs      sequence of embedded systems, of range 0..topdim;
  --   pts      generic points in the witness sets, for range 0..topdim;
  --   topdim   the top dimension of the solution set;
  --   nbl      limit on the number of iterations which keep the
  --            factorization unchanged;
  --   tol      tolerance to decide whether points are equal.

  -- ON RETURN :
  --   f        irreducible decomposition of a solution set,
  --            as an array of range 1..topdim;
  --   times    elapsed CPU user time at each level;
  --   alltime  the total elapsed CPU user time.

end Monodromy_Homotopies;
