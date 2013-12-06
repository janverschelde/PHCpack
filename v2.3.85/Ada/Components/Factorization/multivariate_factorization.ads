with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Multivariate_Factorization is

-- DESCRIPTION :
--   This package offers some routines to factor multivariate polynomials
--   with complex coefficients.  There are three stages :
--     1) with monodromy we group witness points on same factor;
--     2) certification of monodromy groupings with linear traces;
--     3) interpolation at the witness points of the factors.
--   Complementary to stages 1) and 2) is the "Trace_Factor", this
--   combinatorially exploits linear traces to find the factorization.

  procedure Normalize ( p : in out Poly );
  procedure Normalize ( p : in out Poly_Sys );

  -- DESCRIPTION :
  --   A polynomial is normalized if its leading coefficient is one.
  --   This procedure normalizes the polynomial dividing the polynomial
  --   by its leading coefficient.

  procedure Factor 
              ( p : in Poly; n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Link_to_Poly_Sys; fail : out boolean );

  procedure Factor 
              ( p : in Poly; n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );

  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in Poly; n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Link_to_Poly_Sys; fail : out boolean );

  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in Poly; n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );

  -- DESCRIPTION :
  --   Applies monodromy to factor a multivariate polynomial p of degree d.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   output   if true, then continuation will produce extra output,
  --            otherwise the path tracking will run silently;
  --   p        polynomial in n variables of degree d;
  --   n        number of variables;
  --   d        degree of the polynomial.

  -- ON RETURN :
  --   factors  grouping of the labels of the generic points;
  --   mf       mf(i) is the multiplicity of factors(i);
  --   b,v      two vectors of range 1..n, defining a line b+t*v;
  --   wp       d witness points, counted with multiplicity,
  --            with duplicates removed length of vector <= d,
  --            every entry in wp is a t-value for the line b+t*v.
  --   mw       mw(i) is the multiplicity of wp(i);
  --   rdp      random derivatives of p, up to the maximal multiplicity;
  --   rad      rad(i) is distance to closest other point in cluster,
  --            when wp(i) is single, then rad(i) is zero;
  --   dst      dst(i) is distance to closest other point outside cluster;
  --   fail     calculation of witness points failed.

  procedure Certify
              ( p : in Poly;                       
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Link_to_Poly_Sys;
                maxdif : out double_float );
  procedure Certify
              ( file : in file_type; p : in Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Link_to_Poly_Sys;
                maxdif : out double_float );

  -- DESCRITPION :
  --   Certifies the factorization of p using linear traces.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   p        polynomial in n variables;
  --   b        offset vector of a random line b + t*v;
  --   v        direction of a random line b + t*v;
  --   wp       witness points on the line satisfying p(x) = 0;
  --   mw       mw(i) is the multiplicity of w(i);
  --   f        groups witness points according to factors;
  --   mf       mf(i) is the multiplicity of factor f(i);
  --   rdp      random derivatives of p, output of Factor.

  -- ON RETURN :
  --   maxdif   maximal difference between value at linear trace
  --            and the sum of one component at all solutions.

  procedure Interpolate
              ( p : in Poly;                       
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Link_to_Poly_Sys;
                factors : out Link_to_Poly_Sys );
  procedure Interpolate
              ( file : in file_type; p : in Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Link_to_Poly_Sys;
                factors : out Link_to_Poly_Sys );

  -- DESCRITPION :
  --   Certifies the factorization of p using linear traces.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   p        polynomial in n variables;
  --   b        offset vector of a random line b + t*v;
  --   v        direction of a random line b + t*v;
  --   wp       witness points on the line satisfying p(x) = 0;
  --   mw       mw(i) is the multiplicity of w(i);
  --   f        groups witness points according to factors;
  --   mf       mf(i) is the multiplicity of factor f(i);
  --   rdp      random derivatives of p, output of Factor.

  -- ON RETURN :
  --   factors  polynomials interpolating at the factors,
  --            the polynomials on return are normalized.

  function Multiply ( factors : Poly_Sys;
                      mu : Standard_Natural_Vectors.Vector ) return Poly;

  -- DESCRIPTION :
  --   Returns the result of the multiplication of the factors,
  --   with their multiplicities in the vector mu.

  procedure Trace_Factor
              ( p : in Poly; n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Trace_Factor
              ( file : in file_type; p : in Poly; n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );

  -- DESCRIPTION :
  --   Combinatorial factorization using linear traces only.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   p        polynomial in n variables to be factored;
  --   n        number of variables of the polynomial p;
  --   d        degree of the polynomial p.

  -- ON RETURN :
  --   factors  groups witness points according to factors;
  --   mf       mf(i) is the multiplicity of factors(i);
  --   b        offset vector of a random line b + t*v;
  --   v        direction of a random line b + t*v;
  --   wp       witness points on the line satisfying p(x) = 0;
  --   mw       mw(i) is multiplicity of witness point wp(i);
  --   rdp      random derivatives of p, up to the maximal multiplicity;
  --   rad      rad(i) is distance to closest other point in cluster,
  --            when wp(i) is single, then rad(i) is zero;
  --   dst      dst(i) is distance to closest other point outside cluster;
  --   fail     true if failed to compute witness points.

end Multivariate_Factorization;
