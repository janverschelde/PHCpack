with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_c2pieri ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to the Pieri homotopy state machine.
--   The verbose level is given by the value of vrblvl.

-- ON ENTRY :
--   job =  0 : display the menu of all available options;
--   job =  1 : initialize dimensions m,p,q calling use_c2pieri
--              with a[0] = m, a[1] = p, and a[2] = q;
--   job =  2 : initialize input planes calling use_c2pieri
--              with a[0] = m, a[1] = p, and a[2] = q,
--                   b = m*p + q*(m+p), 
--                   c = coefficients of the input planes,
--              in pairs of doubles (real + imaginary part of complex),
--              and as matrices stored row wise;
--   job =  3 : initialize interpolation points calling use_c2pieri
--              with a[0] = m, a[1] = p, and a[2] = q,
--                   b = m*p + q*(m+p), 
--                   c = coefficients of the interpolation points,
--              in pairs of doubles (real + imaginary part of complex);
--   job =  4 : store the start pivots calling use_c2pieri
--              with a[0] = p, and 
--                   b = the 2*p entries of top and bottom pivots;
--   job =  5 : store the target pivots calling use_c2pieri
--              with a[0] = p, and 
--                   b = the 2*p entries of top and bottom pivots;
--   job =  6 : store the start solution curve calling use_c2pieri
--              with a[0] = degree of freedom in start pivots, and
--                   c = coefficients of the complex vector in pairs
--              of doubles (real + imaginary parts);
--   job =  7 : retrieve target solution curve calling use_c2pieri
--              with a[0] = degree of freedom in target pivots,
--              on return are then in c the coefficients of the target;
--   job =  8 : track solution path without intermediate output;
--   job =  9 : track solution path with output diagnostics;
--   job = 10 : verify intersection conditions without output,
--              on return c contains sum of residuals;
--   job = 11 : verify intersection conditions with extra output,
--              on return c contains sum of residuals;
--   job = 12 : destroy the state machine.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
