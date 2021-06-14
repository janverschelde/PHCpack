with text_io;                            use text_io;
with C_Integer_Arrays;                   use C_Integer_Arrays;

package C_Integer_Arrays_io is

-- DESCRIPTION :
--   This package provides output routines for arrays of C integers.

  procedure put ( n : in natural; v : in C_Integer_Array );
  procedure put ( file : in file_type;
                  n : in natural; v : in C_Integer_Array );

  -- DESCRIPTION :
  --   Writes the entries of v from 0 to n-1 to standard output or to file.

end C_Integer_Arrays_io;
