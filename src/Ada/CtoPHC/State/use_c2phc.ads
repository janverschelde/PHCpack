with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_c2phc ( job : integer32;
                     a : C_intarrs.Pointer;
		     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to the operations in PHCpack.
--   See the job codes in the file use_c2phc4c.ads

-- Processing the output data of DEMiCs :
--
--   job  =  834 : allocates memory for the lifting values,
--                 in a[0] are the number of distinct supports,
--                 and b is an array of a[0] integers,
--                 where b[k] is the number of points in the support k,
--                 starting the count k at zero.
--        =  835 : assigns the lifting value of point in support
--                 with index a[0] (starting to count at zero) 
--                 at position with index b[0] (starting again at zero),
--                 to the value given in c[0];
--        =  836 : retrieves the lifting value of point in support
--                 with index a[0] (starting to count at zero) 
--                 at position with index b[0] (starting again at zero),
--                 returns the lifting value in c[0];
--        =  837 : deallocates the memory to store the lifting values;
--        =  838 : appends the string representation of cell indices,
--                 in a[0] are the number of characters and
--                 in b are the integers corresponding to the characters.
--        =  839 : returns the string representation of cell indices,
--                 stored at position a[0], on return in a[0] are the
--                 number of characters in b, the integers representing
--                 the characters in the string;
--        =  840 : deallocates the memory occupied by the cell indices;
--        =  841 : stores the mixed volume given in a;
--        =  842 : returns the stored mixed volume with an assignment to a.
--        =  843 : computes the mixed volume for the system in the systems
--                 container and fills the cells container with the regular
--                 mixed-cell configuration constructed for the mixed volume,
--                 using dynamic enumeration for all mixed cells (demics),
--                 if the standard systems container is empty, then the
--                 system in the standard Laurent systems container is taken;
--        =  844 : computes the mixed volume and the stable mixed volume
--                 for the system in the systems container and fills the 
--                 cells container with the regular mixed-cell configuration
--                 constructed for the stable mixed volume,
--                 using dynamic enumeration for all mixed cells (demics),
--                 if the standard systems container is empty, then the
--                 system in the standard Laurent systems container is taken.
--
-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
