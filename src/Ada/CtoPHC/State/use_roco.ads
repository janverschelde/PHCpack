with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_roco ( job : integer32;
                    a : C_intarrs.Pointer;
                    b : C_intarrs.Pointer;
                    c : C_dblarrs.Pointer;
                    vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to linear-product root counts
--   and the random linear-product systems.

-- ON ENTRY :
--   job     =  0 : constructs a supporting set structure for the system
--                  in the systems container;
--           =  1 : writes the set structure to screen;
--           =  2 : returns in a[0] the root count based on the set structure;
--           =  3 : puts in the system container a random linear-product
--                  system based on the supporting set structure;
--           =  4 : puts in the solution container all solution of
--                  the random linear-product system;
--           =  5 : clears the set structure;
--           =  6 : returns in b the string representation of the current
--                  set structure and in a[0] the number of characters in b.
--           =  7 : parses string in b, with number of characters in a[0]
--                  into a set structure;
--           =  8 : verifies whether the set structure supports the
--                  standard polynomial system in the container,
--                  and return 1 in a[0] if true, 0 in a[0] if false,
--                  if a[0] equals 1 on input, then extra information
--                  is written to screen;
--           = 10 : creates a partition for the system in the container
--                  and returns this partition as a string in b,
--                  a[0] contains the Bezout number and
--                  a[1] the number of characters in the string b;
--           = 11 : given in a[0] the number of characters in the string b
--                  of a partition of the set of unknowns, returns in a[0]
--                  the m-homogeneous Bezout number of the system in the
--                  standard containers with respect to the given partition;
--           = 12 : puts in the systems container a random linear-product
--                  system based on a give partition, apply then job 4
--                  to solve this linear-product system.

-- ON RETURN :
--   0 if operation was successful, otherwise something went wrong...
