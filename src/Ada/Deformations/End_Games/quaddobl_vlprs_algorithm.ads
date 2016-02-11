with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Vectors;                use Quad_Double_Vectors;
with Quad_Double_Matrices;               use Quad_Double_Matrices;

package QuadDobl_vLpRs_Algorithm is

-- DESCRIPTION :
--   This package provides several routines for an extrapolation method
--   on the power series x(s) = a s^w ( 1 + O(s) ),
--   using quad double floating-point arithmetic.
--   Because the aim is to estimate the leading power w, the algorithm
--   takes the logs as input.  So it considers the sequence
--     log(|x(s_k)|) = log(|a|) + w log(s_k) + log(1 + O(s_k)),
--   for a sequence s_0 > s_1 > .. > s_m > 0, with m >= the order r.

-- USAGE :
--   The first two routines are interesting to fill up the pipe, until
--   sufficient data points are available to match the order.
--   From then on, the fully incremental version can be used.

-- START UP FROM SCRATCH :

  procedure vLpRs_full
                ( r : in integer32; s,logs,logx : in Vector;
                  srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix );

  -- DESCRIPTION :
  --   Calls the full version of the vLpRs-Algorithm.

  -- REQUIRED :
  --   s'range = logs'range = logx'range, s'last >= r,
  --   srp'range = dsp'range = 1..r-1 = rt1'range(*) = rt2'range(*),
  --   p'range = 0..r-1, L'range = v'range = 0..r.

  -- ON ENTRY :
  --   r          order of the extrapolation method;
  --   s          strictly decreasing sequence of positive numbers;
  --   logs       logarithms of the s-values;
  --   logx       logarithms of |x(s_k)|.

  -- ON RETURN :
  --   srp        last row of powers of consecutive s-values;
  --   dsp        last row of differences of consecutive s-powers;
  --   p          last row of p-factors;
  --   L          last row in the L-table, with error O(s^r).
  --   v          last row in the v-table, with error O(s^r).
  --   rt1,rt2    last consecutive R-tables.

  procedure vLpRs_pipe
                ( r : in integer32; s,logs,logx : in Vector;
                  srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix );
  procedure vLpRs_pipe
                ( file : in file_type;
                  r : in integer32; s,logs,logx : in Vector;
                  srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix );

  -- DESCRIPTION :
  --   Constructs the extrapolation table in an incremental way.

  -- REQUIRED :
  --   s'range = logs'range = logx'range, s'last >= r,
  --   srp'range = dsp'range = 1..r-1 = rt1'range(*) = rt2'range(*),
  --   p'range = 0..r-1, L'range = v'range = 0..r.

  -- ON ENTRY :
  --   file       to write error table on;
  --   r          order of the extrapolation method;
  --   s          strictly decreasing sequence of positive numbers;
  --   logs       logarithms of the s-values;
  --   logx       logarithms of |x(s_k)|.

  -- ON RETURN :
  --   srp        last row of powers of consecutive s-values;
  --   dsp        last row of differences of consecutive s-powers;
  --   p          last row of p-factors;
  --   L          last row in the L-table, with error O(s^r).
  --   v          last row in the v-table, with error O(s^r).
  --   rt1,rt2    last consecutive R-tables.

-- INCREMENTAL UPDATE :

  procedure vLpRs_pipe
                ( s,logs,logx : in quad_double;
                  srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix );
  
  procedure vLpRs_pipe
                ( file : in file_type; s,logs,logx : in quad_double;
                  srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix );

  -- DESCRIPTION :
  --   One additional row of every table is computed.
  
  -- REQUIRED :
  --   srp'range = dsp'range = 1..r-1 = rt1'range(*) = rt2'range(*),
  --   p'range = 0..r-1, L'range = v'range = 0..r.

  -- ON ENTRY :
  --   file       to write table with errors on;
  --   s          new s-value, must be smaller than latest one and nonzero;
  --   logs       logarithm of s;
  --   logx       logarithm of the absolute value of |x(s)|;
  --   srp        last row of powers of consecutive s-values;
  --   dsp        last row of differences of consecutive s-powers;
  --   p          last row of p-factors;
  --   L          last row of extrapolated logarithms of s-values;
  --   v          last row of extrapolated logarithms of data points;
  --   rt1,rt2    last consecutive R-tables.

  -- ON RETURN :
  --   srp        updated row of powers of consecutive s-values;
  --   dsp        updated row of differences of consecutive s-powers;
  --   p          updated row of p-factors;
  --   L          updated row of extrapolated logarithms of s-values;
  --   v          updated row of extrapolated logarithms of data points;
  --   rt1,rt2    updated consecutive R-tables.
  
end QuadDobl_vLpRs_Algorithm;
