with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Vectors;                use Quad_Double_Vectors;

package QuadDobl_Winding_Numbers is

-- DESCRIPTION :
--   This package provides facilities to estimate the winding
--   number of a power series in s by extrapolation, given 
--   the logarithms of the absolute values of the power series
--   evaluated at a decreasing sequence of values for s,
--   decreasing geometrically at a fixed rate.
--   All calculations are done with quad double arithmetic.

  function Consecutive_Differences ( logx : in Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the vector of consecutive differences of the logarithms
  --   in logx.  The vector on return has range logx'first..logx'last-1.

  function Consecutive_Errors ( difs : in Vector ) return Vector;

  -- DESCRIPTION :
  --   Given the vector of consecutive differences, on return is
  --   the vector of errors on these consecutive differences.
  --   The vector on return has range difs'first..difs'last-1.

  procedure Write_Extrapolation_Errors
              ( file : in file_type; e : in Vector;
                log10h : in quad_double; ew : in integer32 );

  -- DESCRIPTION :
  --   Writes the estimate for m along with the errors
  --   at the end of the extrapolation process.

-- DRIVER PROCEDURES :

  procedure Extrapolate_on_Errors_full
              ( logx : in Vector; h,log10h : in quad_double;
                csc_dif,csc_err,ext_err : out Vector; 
                ew : out integer32; error : out quad_double );
  procedure Extrapolate_on_Errors_full
              ( file : in file_type;
                logx : in Vector; h,log10h : in quad_double;
                csc_dif,csc_err,ext_err : out Vector;
                ew : out integer32; error : out quad_double );

  -- DESCRIPTION :
  --   Does the full extrapolation on the complete vector logx.
  --   So the size of the vector logx determines the order of
  --   the extrapolation method.

  -- REQUIRED :
  --   csc_dif'range = logx'first..logx'last-1,
  --   csc_err'range = logx'first..logx'last-2,
  --   est_err'range = logx'first..logx'last-3.

  -- ON ENTRY :
  --   file     for extra output and diagnostics,
  --            if omitted, then no output will be written;
  --   logx     logarithms of absolute values of power series
  --            evaluated at a decreasing sequences of s values;
  --   h        ratio at which values for s are decreasing;
  --   log10h   decimal logarithm of h.

  -- ON RETURN 
  --   csc_dif  consecutive differences of logarithms in logx;
  --   csc_err  consecutive errors, differences of differences of
  --            consecutive logarithms in logx,
  --   ext_err  extrapolated errors;
  --   ew       estimate for the winding number;
  --   error    error on the estimate for the winding number.

end QuadDobl_Winding_Numbers;
