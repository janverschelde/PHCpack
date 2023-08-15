with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;

package Double_Rho_Algorithm is

-- DESCRIPTION :
--   Implements the rho algorithm on double precision real numbers,

  procedure Allocate
              ( tab : out Standard_Floating_VecVecs.VecVec;
                idx : out Standard_Integer_Vectors.Vector;
                dim : in integer32 );

  -- DESCRIPTION :
  --   Allocates every column in tab with dim doubles,
  --   all initialized to zero.
  --   The vector idx will store the last valid index of each column in tab
  --   and is initialized to -1.

  -- REQUIRED : tab'range = idx'range.

  procedure Initialize
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_float;
                verbose : in boolean := true ); 

  -- DESCRIPTION :
  --   Updates the vector of columns tab with the number nbr.
  --   Works only to fill up the first two columns with the
  --   first two numbers of a sequence.
  --   If verbose, then each updated column is printed.

  -- REQUIRED : tab'range = idx'range.

  -- ON ENTRY :
  --   tab      table of extrapolated values;
  --   dim      equals the number of defined columns in the table;
  --   idx      idx(k) is the last index of column k;
  --   nbr      new number to update the table with;
  --   verbose  is the verbose flag.

  -- ON RETURN :
  --   tab      table with up to the first two columns updated;
  --   dim      updated dimension of the table;
  --   idx      updated index for each column.

  procedure Columns 
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_float;
                verbose : in boolean := true ); 

  -- DESCRIPTION :
  --   Updates the defined columns in tab with the number nbr.

  -- REQUIRED : tab'range = idx'range.

  -- ON ENTRY :
  --   tab      table of extrapolated values;
  --   dim      equals the number of defined columns in the table;
  --   idx      idx(k) is the last index of column k;
  --   nbr      new number to update the table with;
  --   verbose  is the verbose flag.

  -- ON RETURN :
  --   tab      table with up to date columns;
  --   idx      updated index for each column.

  procedure New_Column
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ); 

  -- DESCRIPTION :
  --   Adds a new column to tab.
  --   Must be called after the Columns procedure.

  -- REQUIRED : tab'range = idx'range and dim >= 2.

  -- ON ENTRY :
  --   tab      table of extrapolated values;
  --   dim      equals the number of defined columns in the table;
  --   idx      idx(k) is the last index of column k;
  --   verbose  is the verbose flag.

  -- ON RETURN :
  --   tab      table with new columns if possible;
  --   dim      updated dimension of the table;
  --   idx      updated index for each column.

  procedure Extrapolate
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_float;
                verbose : in boolean := true ); 

  -- DESCRIPTION :
  --   Updates the defined columns in tab with the number nbr,
  --   calling Initialize for the first two numbers, and
  --   calling New_Column after Columns for the next numbers.

  -- REQUIRED : tab'range = idx'range.

  -- ON ENTRY :
  --   tab      table of extrapolated values;
  --   dim      equals the number of defined columns in the table;
  --   idx      idx(k) is the last index of column k;
  --   nbr      new number to update the table with;
  --   verbose  is the verbose flag.

  -- ON RETURN :
  --   tab      table with up to date columns;
  --   dim      updated dimension of the table;
  --   idx      updated index for each column.

end Double_Rho_Algorithm;
