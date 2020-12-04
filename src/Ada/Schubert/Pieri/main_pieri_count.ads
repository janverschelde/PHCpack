with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Main_Pieri_Count is

-- DESCRIPTION :
--   Export the main program to compute a Pieri root count.

  procedure Pieri_Count ( file : in file_type; m,p,q : in natural32 );

  -- DESCRIPTION :
  --   Writes the Pieri root count for all curves of degree q
  --   that produce p-planes that meet m-planes at distinct places.
  
  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -e -b.
  --   Retrieves m, p, and q from the input file
  --   and writes the pieri count to the output file.

  -- ON INPUT :
  --   infilename      name of the file where the input is;
  --   outfilename     name of the output file;
  --   verbose         the verbose level.

end Main_Pieri_Count;
