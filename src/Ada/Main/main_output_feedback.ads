with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Main_Output_Feedback is

-- DESCRIPTION :
--   Encapsulates the computation of output feedback laws,
--   computed by Pieri homotopies, and realized by algorithms
--   written by Yusong Wang.

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -k.
  --   The first two arguments are the names of the input and output files.
  --   The last argument is the verbose level.

end Main_Output_Feedback;
