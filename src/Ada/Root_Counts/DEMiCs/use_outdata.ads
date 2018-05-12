with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_outdata ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to DEMiCs_Output_Data.

-- ON ENTRY :
--   job   =  0 : allocates memory for the lifting values,
--                in a[0] are the number of distinct supports,
--                and b is an array of a[0] integers,
--                where b[k] is the number of points in the support k,
--                starting the count k at zero.
--         =  1 : assigns the lifting value of point in support
--                with index a[0] (starting to count at zero) 
--                at position with index b[0] (starting again at zero),
--                to the value given in c[0];
--         =  2 : retrieves the lifting value of point in support
--                with index a[0] (starting to count at zero) 
--                at position with index b[0] (starting again at zero),
--                returns the lifting value in c[0].
--         =  3 : deallocates the memory to store the lifting values;
--         =  4 : appends the string representation of cell indices,
--                in a[0] are the number of characters and
--                in b are the integers corresponding to the characters.
--         =  5 : deallocates the memory occupied by the cell indices.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
