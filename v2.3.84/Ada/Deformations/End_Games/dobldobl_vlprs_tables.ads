with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Vectors;              use Double_Double_Vectors;
with Double_Double_Matrices;             use Double_Double_Matrices;

package DoblDobl_vLpRs_Tables is

-- DESCRIPTION :
--   This package implements the tables for a general r-order 
--   extrapolation method, using double double floating-point arithmetic.
--   For all tables, there is a full and a pipe-lining version.
--   For efficiency reasons the full v- and L-table are combined.

-- I. The v-table : extrapolated values of log(|x_i(s_k)|).

  procedure v_pipe ( v : in out Vector; p : in Vector; vr : in double_double );
  procedure v_pipe ( v,p : in Vector; vr : in double_double;
                     vrp : in out Vector );

  -- DESCRIPTION : computes one additional row of the v-table.

  -- REQUIRED : v'range = vrp'range = 0..r, p'range = 0..r-1.

  -- ON ENTRY :
  --   v             previous row of the v-table;
  --   p             current vector of p-factors;
  --   vr            equals log(|x_i(s_r)|).

  -- ON RETURN :
  --   v or vrp      last row in the v-table: vrp(r) should have error O(s^r).

-- II. The L-table : linear combinations of the logarithms log(s_k).

  procedure L_pipe ( l : in out Vector; p : in Vector; lr : in double_double );
  procedure L_pipe ( l,p : in Vector; lr : in double_double;
                     lrp : in out Vector );

  -- DESCRIPTION : computes one additional row of the L-table.

  -- REQUIRED : l'range = lrp'range = 0..r, p'range = 0..r-1.

  -- ON ENTRY :
  --   l             previous row of the v-table;
  --   p             current vector of p'factors;
  --   lr            equals log(s_r).

  -- ON RETURN :
  --   l or lrp      last row in the L-table, vrp(r)/lpr(r) = w_i + O(s^r).

-- The full computation of both v- and L-table :

  procedure vL_full ( s,l,v : in Vector; srp,dsp,p,lrp,vrp : out Vector;
                      rt1,rt2 : in out Matrix );

  -- DESCRIPTION : computes the last row of v- and L-table from s-values.

  -- REQUIRED :
  --   s'range = l'range = v'range = 0..r = lrp'range = vrp'range, 
  --   srp'range = dsp'range = 1..r-1, p'range = 0..r-1,
  --   rt1'range(*) = 1..r-1 = rt2'range(*).

  -- ON ENTRY :
  --   s             consecutive s-values: s(0) > s(1) > .. > s(r) > 0;
  --   l             logs of the s-values: l(k) = log(s(k));
  --   v             points at the s-values: v(k) = log(|x_i(s(k))|).

  -- ON RETURN :
  --   srp           consecutive powers of s(r): srp(l) = s(r)**l;
  --   dsp           differences of powers: dsp(l) = srp(l) - s(r-1)**l.
  --   p             last row used in the p-table;
  --   vrp           last row in the v-table;
  --   lrp           last row in the L-table, vrp(r)/lpr(r) = w_i + O(s^r);
  --   rt1           previous instance of the r-table;
  --   rt2           last instance of the r-table.

-- III. The p-table : extrapolation factors needed for v-table and L-table.

  procedure p_full ( s : in Vector; srp,dsp,p : out Vector;
                     rt1,rt2 : in out Matrix );

  -- DESCRIPTION : computes the last row of p-table from s-values.

  -- REQUIRED :
  --   s'range = 0..r, srp'range = dsp'range = 1..r-1, p'range = 0..r-1,
  --   rt1'range(*) = 1..r-1 = rt2'range(*).

  -- NOTE : p_full = R_full followed by p_pipe.

  procedure p_pipe ( rt1,rt2 : in Matrix; p : out Vector );

  -- DESCRIPTION : update of one row of the p-table.

  -- REQUIRED : rt1'range(*) = rt2'range(*) = 1..r-1, p'range = 0..r-1.

  -- ON ENTRY :
  --   rt1           previous instance of the R-table;
  --   rt2           last instance of the R-table.

  -- ON RETURN :
  --   p             p-factors, with p(0) = 1, p(i) = rt1(i,i)/rt2(i,i).

-- IV. The R-table : propagated error factors.

  procedure R_full ( s : in Vector; srp,dsp,p : out Vector;
                     rt1,rt2 : in out Matrix );

  procedure RR_full ( s : in Vector; srp,dsp,p : out Vector;
                      rt1,rt2 : in out Matrix );

  -- DESCRIPTION : computation of the r-table from a sequence of s-values.
  --   The RR_full computes the subdiagonal as well.

  -- REQUIRED :
  --   s'range = 0..r, srp'range = dsp'range = 1..r-1, p'range = 0..r-1,
  --   rt1'range(*) = 1..r-1 = rt2'range(*).

  -- ON ENTRY :
  --   s             sequence of consecutive s-values.

  -- ON RETURN :
  --   srp,dsp       see the output of s_full;
  --   p             last row used in the p-table;
  --   rt1           previous instance of the r-table;
  --   rt2           last instance of the r-table.

  procedure R_pipe ( rt1 : in Matrix; s,p : in Vector; rt2 : in out Matrix );
  procedure RR_pipe ( rt1 : in Matrix; s,p : in Vector; rt2 : in out Matrix );

  -- DESCRIPTION : update of the r-table.  RR_pipe updates also subdiagonal.

  -- REQUIRED :
  --   rt1'range(*) = rt2'range(*) = s'range = 1..r-1, p'range = 0..r-1.

  -- ON ENTRY :
  --   rt1           last instance of the R-table;
  --   s             current row of differences of powers of s-values;
  --   p             last row of the current p-table.

  -- ON RETURN :
  --   rt2           new instance of the R-table.

-- V. The s-table : consecutive s-values with differences of their powers.

  procedure s_full ( s : in Vector; srp,dsp : out Vector );

  -- DESCRIPTION : computes the s-table in full.

  -- REQUIRED : s'range = 0..r, srp'range = dsp'range = 1..r-1.

  -- ON ENTRY :
  --   s             the s-values, with s(0) > s(1) > .. > s(r) > 0.

  -- ON RETURN :
  --   srp           consecutive powers of s(r): srp(l) = s(r)**l;
  --   dsp           differences of powers: dsp(l) = srp(l) - s(r-1)**l.

  procedure s_pipe ( srp : in out Vector; sr : in double_double;
                     dsp : out Vector );
  procedure s_pipe ( sr1 : in Vector; sr : in double_double;
                     srp,dsp : out Vector );

  -- DESCRIPTION : computes an additional row of the s-table.

  -- REQUIRED : sr1'range = srp'range = dsp'range = 1..r-1.

  -- ON ENTRY :
  --   sr1 or srp    consecutive powers of s(r-1): srp(l) = s(r-1)**l.
  --   sr            new last value for s(r).

  -- ON RETURN : same as s_full.

end DoblDobl_vLpRs_Tables;
