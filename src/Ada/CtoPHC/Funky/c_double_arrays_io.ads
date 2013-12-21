with text_io;                            use text_io;
with C_Double_Arrays;                    use C_Double_Arrays;

package C_Double_Arrays_io is

-- DESCRIPTION :
--   This package provides output routines for arrays of C doubles.

  procedure put ( n : in natural; v : in C_Double_Array );
  procedure put ( file : in file_type;
                  n : in natural; v : in C_Double_Array );

  -- DESCRIPTION :
  --   Writes the entries of v from 0 to n-1 to standard output or to file.

end C_Double_Arrays_io;
