with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Vectors;          use Multprec_Floating_Vectors;

package Multprec_Winding_Numbers is

-- DESCRIPTION :
--   This package provides facilities to estimate the winding
--   number of a power series in s by extrapolation, given 
--   the logarithms of the absolute values of the power series
--   evaluated at a decreasing sequence of values for s,
--   decreasing geometrically at a fixed rate.
--   The procedures work with multiprecision arithmetic.

  function Consecutive_Errors ( logx : in Vector ) return Vector;

  -- DESCRIPTION :
  --   Computes the vector of errors on consecutive differences
  --   of the logarithms in logx.  If logx has range 0..r,
  --   then the vector on return has range 0..r-2.

  procedure Write_Extrapolation_Errors
              ( file : in file_type; e : in Vector;
                log10h : in Floating_Number; ew : in integer32 );

  -- DESCRIPTION :
  --   Writes the estimate for m along with the errors
  --   at the end of the extrapolation process.

-- DRIVER PROCEDURES :

  procedure Extrapolate_on_Errors
              ( logx : in Vector; h : in Floating_Number;
                ew : out integer32; error : out Floating_Number );
  procedure Extrapolate_on_Errors
              ( file : in file_type;
                logx : in Vector; h : in Floating_Number;
                ew : out integer32; error : out Floating_Number );

  -- ON ENTRY :
  --   file     for extra output and diagnostics,
  --            if omitted, then no output will be written;
  --   logx     logarithms of absolute values of power series
  --            evaluated at a decreasing sequences of s values;
  --   h        ratio at which values for s are decreasing.

  -- ON RETURN 
  --   ew       estimate for the winding number;
  --   error    error on the estimate for the winding number.

end Multprec_Winding_Numbers;
