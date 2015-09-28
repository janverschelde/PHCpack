with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_Select_Solutions is

-- DESCRIPTION :
--   This package offers routines to select solutions 
--   computed with double double precision arithmetic from file.

  procedure Scan_Solutions
              ( file : in file_type; len,dim : in natural32;
                sv : in Vector; sa : out Solution_Array );

  -- DESCRIPTION :
  --   Scans the file for a selection of solutions.

  -- REQUIRED :
  --   The file must be opened for input and placed after reading
  --   of the dimension information; sv'range = sa'range.

  -- ON ENTRY :
  --   file     input file for the solution list;
  --   len      total number of solutions on the file;
  --   dim      length of each solution vector on the file;
  --   sv       selection of numbers to be scanned from file.

  -- ON RETURN :
  --   sa       contains the selected solutions from file.

  procedure Scan_Solutions
              ( file : in file_type; len,dim : in natural32;
                sv : in Vector; sols : out Solution_List );

  -- DESCRIPTION :
  --   Scans the file for a selection of solutions.

  -- REQUIRED :
  --   The file must be opened for input and placed after reading
  --   of the dimension information; sv'range = sa'range.

  -- ON ENTRY :
  --   file     input file for the solution list;
  --   len      total number of solutions on the file;
  --   dim      length of each solution vector on the file;
  --   sv       selection of numbers to be scanned from file.

  -- ON RETURN :
  --   sols     contains the selected solutions from file.

  function Sort ( v : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the sorted list of numbers in the vector v,
  --   sorted in increasing order.

  function Find ( a : natural32; v : Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns 0 if the value a does not occur in the vector v,
  --   otherwise returns k when v(k) = a.

  procedure Write_Selection
               ( file : in file_type; dim : in natural32;
                 rv,sv : in Vector; sa : in Solution_Array );

  -- DESCRIPTION :
  --   Writes the selected sequence of solutions to file,
  --   in the order as the selected sequence rv.

  -- ON ENTRY :
  --   file      output file opened for output;
  --   dim       dimension of the solutions in sa;
  --   rv        selected sequence of solutions;
  --   sv        sorted sequence of solutions;
  --   sa        the selected sequence sorted as in sv.

  procedure Store_Selection
               ( sols : out Solution_List;
                 rv,sv : in Vector; sa : in Solution_Array );

  -- DESCRIPTION :
  --   Stores the selected sequence of solutions in sa to sols,
  --   in the order as the selected sequence rv.

  -- ON ENTRY :
  --   rv        selected sequence of solutions;
  --   sv        sorted sequence of solutions;
  --   sa        the selected sequence sorted as in sv.

  -- ON RETURN :
  --   sols      list of selected solutions.

  procedure Select_Solutions
              ( sols : in Solution_List;
                rv : in Standard_Natural_Vectors.Vector;
                sv : out Standard_Natural_Vectors.Vector;
                sa : out Solution_Array );

  -- DESCRIPTION :
  --   Selects those solutions in sols that occur at the positions
  --   indicated by the indices in rv.

  -- REQUIRED : rv'range = sv'range = sa'range
  --   for all i in rv'range: rv(i) <= Length_Of(sols).

  -- ON ENTRY :
  --   sols     list of solutions;
  --   rv       indices of elements in sols.

  -- ON RETURN :
  --   sv       sorted vector rv;
  --   sa       selected solutions from sols.

  procedure Select_from_File
              ( file : in file_type;
                rv : in Standard_Natural_Vectors.Vector;
                sa : out Solution_Array; fail : out boolean );
  procedure Select_from_File
              ( file : in file_type;
                rv : in Standard_Natural_Vectors.Vector;
                sv : out Standard_Natural_Vectors.Vector;
                sa : out Solution_Array; fail : out boolean );
  procedure Select_from_File
              ( file : in file_type;
                rv : in Standard_Natural_Vectors.Vector;
                sols : out Solution_List );
  procedure Select_from_File
              ( file : in file_type; bannered : in boolean;
                rv : in Standard_Natural_Vectors.Vector;
                sv : out Standard_Natural_Vectors.Vector;
                sa : out Solution_Array; fail : out boolean );
  procedure Select_from_File
              ( file : in file_type; bannered : in boolean;
                rv : in Standard_Natural_Vectors.Vector;
                sols : out Solution_List );

  -- DESCRIPTION :
  --   Scans the file for solutions and selects those solutions
  --   with indices in rv.  The result is in the list sols.

  -- REQUIRED :
  --   The file should be positioned at the beginning,
  --   prior to the banner.

  -- ON ENTRY :
  --   file     input file where the solutions are;
  --   bannered is true if the solutions are preceeded by a banner,
  --            if omitted, then the user will be prompted;
  --   rv       indices of solutions to select from file.

  -- ON RETURN :
  --   sv       sorted vector rv;
  --   sa       selected solutions as array;
  --   sols     selected solutions as list;
  --   fail     true if something went wrong...

end QuadDobl_Select_Solutions;
