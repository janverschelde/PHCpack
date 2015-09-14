with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;

package Standard_Lined_Hypersurfaces is

-- DESCRIPTION :
--   A generic point on a hypersurface is the intersection of the hypersurface
--   with a general line.  A 'lined hypersurface' is a general line with
--   generic points on that line and the hypersurface.
--   The routines in this package produce generic points in standard double
--   precision by varying the general line.

-- PART I : computing generic points on a random line.

  procedure Generic_Points
                ( file : in file_type;
                  p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean );
  procedure Generic_Points
                ( file : in file_type;
                  p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean;
                  m : out Standard_Natural_Vectors.Vector );
  procedure Generic_Points
                ( file : in file_type;
                  p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean;
                  m : out Standard_Natural_Vectors.Vector;
                  rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Generic_Points
                ( file : in file_type;
                  p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean;
                  m : out Standard_Natural_Vectors.Vector;
                  rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  rad,dst : out Standard_Floating_Vectors.Vector );
  procedure Generic_Points
                ( p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean );
  procedure Generic_Points
                ( p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean;
                  m : out Standard_Natural_Vectors.Vector );
  procedure Generic_Points
                ( p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean;
                  m : out Standard_Natural_Vectors.Vector;
                  rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Generic_Points
                ( p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  d : in natural32;
                  b,v : in Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_Vectors.Vector;
                  fail : out boolean;
                  m : out Standard_Natural_Vectors.Vector;
                  rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  rad,dst : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Intersecting a degree d polynomial p with a general line will
  --   lead to d generic points on the hypersurface p(x) = 0.

  -- ON ENTRY :
  --   file       for diagnostics and intermediate output;
  --   p          polynomial is several variables;
  --   ep         polynomial p in nested Horner form for evaluation;
  --   d          degree of the polynomial p;
  --   b          base point for general line x(t) = b + t*v;
  --   v          direction of general line;
  --   eps        accuracy requirement;
  --   maxit      maximal number of iterations.

  -- ON RETURN :
  --   t          vector of range 1..d with t values for generic points :
  --                p(b + t(i)*v) = 0, for i in 1..d;
  --   fail       true if failed to converge, false otherwise;
  --   m          m(i) is the multiplicity for the m-th root;
  --   rdp        rdp(m(i)) contains the (m(i)-1)-th random derivative
  --              of the original polynomial p;
  --   rad        distance of multiple point to nearest closest other
  --              point in the cluster, is zero for single point;
  --   dst        distance of point to nearest point outside cluster.

  procedure Cluster_Analysis
	        ( x : in Standard_Complex_Vectors.Vector;
                  tol : in double_float;
                  radius,distance : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs a cluster analysis on the points in x.  Points of x belong
  --   to the same cluster if they lie within distance tol from each other.
  --   For multiple points, the radius of the cluster and distance to the
  --   other points in x is computed.

-- PART II : sampling and refining along moving lines

  procedure Silent_Refiner
                ( p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                  b,v : in Standard_Complex_Vectors.Vector;
                  t : in out Standard_Complex_Vectors.Vector;
                  ft,dt : out Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32 );
  procedure Silent_Refiner
                ( p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                  b,v : in Standard_Complex_Vectors.Vector;
                  m : in Standard_Natural_Vectors.Vector;
                  t : in out Standard_Complex_Vectors.Vector;
                  ft,dt : out Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32 );
  procedure Reporting_Refiner
                ( file : in file_type;
                  p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                  b,v : in Standard_Complex_Vectors.Vector;
                  t : in out Standard_Complex_Vectors.Vector;
                  ft,dt : out Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32 );
  procedure Reporting_Refiner
                ( file : in file_type;
                  p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                  b,v : in Standard_Complex_Vectors.Vector;
                  m : in Standard_Natural_Vectors.Vector;
                  t : in out Standard_Complex_Vectors.Vector;
                  ft,dt : out Standard_Complex_Vectors.Vector;
                  eps : in double_float; maxit : in natural32 );

  -- DESCRIPTION :
  --   Applies at most four Newton iterations to each root,
  --   eventually using multiplicities in a modified Newton.

  -- ON ENTRY :
  --   file       to write a conclusion line for each root;
  --   p          vector of Horner forms of multivariate polynomials :
  --                p(0) defines the equation of the hypersurface,
  --                p(i) is the i-th partial derivative;
  --   (b,v)      base point and direction of the line b + t*v;
  --   m          multiplicity information, the i-th root occurs m(i) times;
  --   t          approximations for p(b + t(i)*v) = 0, i=1,2,..Degree(p);
  --   eps        accuracy requirement;
  --   maxit      maximal number of iterations allowed.

  -- ON RETURN :
  --   t          refined approximations;
  --   ft         ft(i) is the function evaluated at t(i);
  --   dt         dt(i) is the last correction term.

  procedure Silent_Hypersurface_Sampler
                ( p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                  b0,v0,b1,v1 : in Standard_Complex_Vectors.Vector;
                  t : in out Standard_Complex_Vectors.Vector );
  procedure Reporting_Hypersurface_Sampler
                ( file : in file_type;
                  p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                  b0,v0,b1,v1 : in Standard_Complex_Vectors.Vector;
                  output : in boolean;
                  t : in out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes new samples of the hypersurface defined by the polynomial p,
  --   of degree d.

  -- ON ENTRY :
  --   file       for intermediate results and diagnostics;
  --   p          vector of Horner forms of multivariate polynomials :
  --                p(0) is the defining equation of the hypersurface,
  --                p(i) is the i-th partial derivative;
  --   (b0,v0)    base point and direction of the start line b0 + t*v0;
  --   (b1,v1)    base point and direction of the target line b1 + t*v1;
  --   output     if true, then predictor-corrector will produce output,
  --              otherwise, the predictor-corrector will be silent;
  --   t          vector with roots p(b0 + t(i)*v0) = 0, for i=1,2,..,d.

  -- ON RETURN :
  --   t          vector with roots p(b1 + t(i)*v1) = 0, for i=1,2,..,d.

end Standard_Lined_Hypersurfaces;
