with Interfaces.C.Strings;              use Interfaces.C.Strings;
with C_Double_Arrays;                   use C_Double_Arrays;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

function Pieri_Solver ( m,p,q,nb,output_level : integer32; 
                        points,planes : C_dblarrs.Pointer;
                        filename : chars_ptr ) return integer32;

-- DESCRIPTION :
--   Returns the number of maps of degree q that produce p-planes
--   meeting n given generic m-planes at n specified interpolation points,
--   where n equals m*p + q*(m+p).  The arrays of doubles contain the real
--   and imaginary parts of the complex coefficients.

-- ON ENTRY :
--   m          dimension of the input planes;
--   p          dimension of the output planes;
--   q          degree of the solution maps;
--   nb         number of solution maps that have to be computed,
--              if < 0 then all solution maps will be computed.
--   output_level determines the amount of intermediate output:
--              if 0, then no intermediate output;
--              if 1, only final determinant validation will be done;
--              if 2, all intermediate determinant validations are done;
--              if 3 (or higher), output of all path tracking is shown;
--   points     values of the interpolation points;
--   planes     m-planes which the solution maps have to meet;
--   filename   name of the output file.

-- ON RETURN :
--   result of the combinatorial root count
--   or -1 if the call to the C routine failed.
