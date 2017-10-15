with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_VecVecs;

package Path_Counts_Table is

-- DESCRIPTION :
--   A path counts table contains the counts after the cascade filtering.
--   Every entry in the table is a vector with three numbers,
--   at every dimension, counting the number of paths tracked,
--   and the number of solutions with zero and nonzero slack variables.
--   The range of the path counts table goes from 0 to topdim,
--   where topdim is the top dimension of the solution set.

  procedure Update_Path_Counts
              ( cnts : in out Standard_Natural_VecVecs.VecVec;
                dim,nsols,nsols0,nsols1 : in natural32 );

  -- DESCRIPTION :
  --   Updates the vector of path counts with the result of the
  --   cascade filtering at the dimension dim.

  -- ON ENTRY :
  --   cnts     vector with a range which includes dim;
  --   dim      the current dimension;
  --   nsols    number of solutions filtered;
  --   nsols0   number of solutions with zero slack variables;
  --   nsols1   number of solutions with nonzero slack variables.

  -- ON RETURN :
  --   cnts     cnts(dim) points to the vector with the three numbers
  --            nsols, nsols0, and nsols1.

  procedure Write_Path_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Writes the path counts in cnts to file.

  procedure Write_Path_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec;
                times : in Array_of_Duration; totaltime : in duration );

  -- DESCRIPTION :
  --   Writes the path counts in cnts to file,
  --   with the elapsed CPU user time at each stage,
  --   and the sums for each count.
  --   The totaltime is the total elapsed CPU user time.

  -- REQUIRED : times'range = cnts'range.

  procedure Write_Filter_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec;
                times : in Array_of_Duration; totaltime : in duration );

  -- DESCRIPTION :
  --   Writes the number of points in each stage of the homotopy filter,
  --   for all dimensions.  The counts in cnts are triangular,
  --   as there is no filtering needed at the top dimension, and at
  --   dimension 0, there are as many stages as the top dimension minus one.

  procedure Write_Factor_Counts
              ( file : in file_type;
                deco : in Standard_Natural_VecVecs.Array_of_VecVecs;
                times : in Array_of_Duration; totaltime : in duration );

  -- DESCRIPTION :
  --   Writes the classification of the witness points along the factors
  --   in the irreducible decomposition, for all dimensions.
  --   The range of deco is 1..topdim.
  --   There is no factorization at dimension 0.

end Path_Counts_Table;
