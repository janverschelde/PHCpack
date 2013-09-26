with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_roco ( job : integer32;
                    a : C_intarrs.Pointer;
                    b : C_intarrs.Pointer;
                    c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to linear-product root counts
--   and the random linear-product systems.

-- ON ENTRY :
--   job      =  0 : constructs a supporting set structure for the system
--                   in the systems container;
--            =  1 : writes the set structure to screen;
--            =  2 : returns in a the root count based on the set structure;
--            =  3 : puts in the system container a random linear-product
--                   system based on the supporting set structure;
--            =  4 : puts in the solution container all solution of
--                   the random linear-product system;
--            =  5 : clears the set structure.

-- ON RETURN :
--   0 if operation was successful, otherwise something went wrong...
