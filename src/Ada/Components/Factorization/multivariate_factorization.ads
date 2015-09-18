with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Systems;

package Multivariate_Factorization is

-- DESCRIPTION :
--   This package offers some routines to factor multivariate polynomials
--   with complex coefficients.  There are three stages :
--     (1) with monodromy we group witness points on same factor;
--     (2) certification of monodromy groupings with linear traces;
--     (3) interpolation at the witness points of the factors.
--   Stage (2) is also applied in the "Trace_Factor", in the combinatorial
--   exploration with linear traces to find the factorization.

-- FACTORING WITH GIVEN GENERIC POINTS :

  procedure Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;     
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;     
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;     
                 mw : out Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Performs the factorization when multiplicities > 1, 
  --   without any intermediate output.

  -- ON ENTRY :
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   ep        nested Horner form of the polynomial;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

  procedure Factor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Factor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Factor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Performs the factorization when multiplicities > 1,
  --   with intermediate output to file.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   output    if output needed during continuation;
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   ep        nested Horner form of the polynomial;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

  procedure Trace_Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Trace_Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Trace_Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the combinatorial factorization with linear traces
  --   to the case with witness points of higher multiplicity,
  --   in standard double, double double, or quad double precision,
  --   without any intermediate output.

  -- ON ENTRY :
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

  procedure Trace_Factor_with_Multiplicities
               ( file : in file_type; n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Trace_Factor_with_Multiplicities
               ( file : in file_type; n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Trace_Factor_with_Multiplicities
               ( file : in file_type; n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the combinatorial factorization with linear traces
  --   to the case with witness points of higher multiplicity,
  --   with intermediate output to file.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   ep        nested Horner form of the polynomial;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

-- FACTORING WITH MONODROMY :

  procedure Factor 
              ( p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                fail : out boolean );
  procedure Factor 
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                fail : out boolean );
  procedure Factor 
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                fail : out boolean );

  procedure Factor 
              ( p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Factor 
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Factor 
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );

  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                fail : out boolean );
  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                fail : out boolean );
  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                fail : out boolean );

  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Factor
              ( file : in file_type; output : in boolean;
                p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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
              ( p : in Standard_Complex_Polynomials.Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                maxdif : out double_float );
  procedure Certify
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                maxdif : out double_float );
  procedure Certify
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                maxdif : out double_float );
  procedure Certify
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                maxdif : out double_float );
  procedure Certify
              ( file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly;
                b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                maxdif : out double_float );
  procedure Certify
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly;
                b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                maxdif : out double_float );

  -- DESCRITPION :
  --   Certifies the factorization of p using linear traces,
  --   in standard double, double double, or quad double precision,
  --   without intermediate output or with output written to file.

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
              ( p : in Standard_Complex_Polynomials.Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Interpolate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Interpolate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Interpolate
              ( file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly;
                b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Interpolate
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly;
                b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRITPION :
  --   Certifies the factorization of p using linear traces,
  --   in standard double, double double, or quad double precision,
  --   without intermediate output or with output written to file.

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

  function Multiply ( factors : Standard_Complex_Poly_Systems.Poly_Sys;
                      mu : Standard_Natural_Vectors.Vector )
                    return Standard_Complex_Polynomials.Poly;
  function Multiply ( factors : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                      mu : Standard_Natural_Vectors.Vector )
                    return DoblDobl_Complex_Polynomials.Poly;
  function Multiply ( factors : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                      mu : Standard_Natural_Vectors.Vector )
                    return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the result of the multiplication of the factors,
  --   in standard double, double double, or quad double precision,
  --   with their multiplicities in the vector mu.

-- COMBINATORIAL EXPLORATION :

  procedure Trace_Factor
              ( p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Trace_Factor
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Trace_Factor
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Trace_Factor
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Trace_Factor
              ( file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean );
  procedure Trace_Factor
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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
