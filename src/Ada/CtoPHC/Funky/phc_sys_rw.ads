with Standard_Integer_Numbers;       use Standard_Integer_Numbers;
with C_Double_Arrays;                use C_Double_Arrays;

function phc_sys_rw ( rw,size_p : integer32; p : C_DblArrs.Pointer )
                    return C_DblArrs.Pointer;

-- DESCRIPTION :
--   Dual read/write procedure to get polynomial system into C or 
--   to export a concatenated coefficient support representation.

-- ON ENTRY :
--   rw        flag to indicate read (rw = 0) or write (rw /= 0);
--   size_p    size of the buffer;
--   p         concatenated coefficient support representation of
--             a polynomial system when rw /= 0, ignored when rw = 0.

-- ON RETURN :
--   if rw = 0, then on return is a concatenated coefficient support
--   representation of a polynomial system, otherwise ignore return.
