with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;          use Multprec_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;

package Multprec_Subspace_Restrictions is

-- DESCRIPTION :
--   This package provides constructors for the equations of the container
--   subspaces and utilities to restrict solutions to these subspaces.

  procedure Container_Dimension
                ( file : in file_type; k,n,size : in natural32;
                  samples : in Multprec_Complex_VecVecs.VecVec;
                  tol : in double_float; mat,vec : in out Matrix;
                  dim : out natural32 );

  -- DESCRIPTION :
  --   Determines the dimension of the space spanned by the samples.

  -- ON ENTRY :
  --   file       to write diagnostics;
  --   k          the number of vectors used in the vector configuration;
  --   n          dimension of the original space, before embedding;
  --   size       size of the numbers;
  --   samples    array of at least n+1 points;
  --   tol        tolerance to determine whether number is zero or not;

  -- ON RETURN :
  --   mat        matrix that contains points in its rows, range 1..n+1,1..n;
  --   vec        triangulated vector configuration, range 1..n,1..n;
  --   dim        dimension of the space spanned by the samples.

  procedure Container_Subspace
                ( file : in file_type; k,n,level,dim : in natural32;
                  p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  embp : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mat,vec : in Matrix; tol : in double_float;
                  samples : in Multprec_Complex_VecVecs.VecVec;
                  restsamp : out Multprec_Complex_VecVecs.VecVec;
                  lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                  kerpols,restp
                    : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                  restembp
                    : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Constructs equations for the linear subspace spanned by points.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   n          original number of variables;
  --   level      number of added variables in the embedding;
  --   dim        dimension of the container subspace;
  --   p          original polynomial system of dimension n;
  --   embp       embedded polynomial system of dimension n+level;
  --   mat        point configuration that spans the subspace;
  --   vec        vector configuration in upper triangular form;
  --   tol        tolerance to decide whether number is zero;
  --   samples    current set of samples.

  -- ON RETURN :
  --   restsamp   are the samples restricted to the subspace;
  --   lpiv       link to pivots, remaining variables after restriction;
  --   kerpols    linear polynomials describing the subspace;
  --   restp      polynomial system p restricted to the subspace;
  --   restembp   restricted embedded polynomial system.

  procedure Subspace_Restriction
                ( file : in file_type;
                  orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                  ind,k,n,level,size : in natural32;
                  samples : in Multprec_Complex_VecVecs.Array_of_VecVecs;
                  restsamp : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
                  lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                  kerpols,restorgsys
                    : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                  restembsys 
                    : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  dim : out natural32 );

  -- DESCRIPTION :
  --   Determines the dimension of the subspace spanned by the samples
  --   and restricts the polynomial system to this subspace if dim < n.

  -- ON ENTRY :
  --   file       output file to write diagnostics and intermediate results;
  --   orgsys     original polynomial system;
  --   ind        index to the current solution samples;
  --   k          size of vector configuration;
  --   n          dimension of the original space before the embedding;
  --   level      number of slices added in the embedding;
  --   size       size of the numbers;
  --   samples    samples(ind) must contain at least k+1 points.

  -- ON RETURN :
  --   restsamp   samples restricted to the subspace;
  --   lpiv       link to the pivots, indicating the remaining variables;
  --   kerpols    defining polynomials of the subspace, if dim < k;
  --   restorgsys is the original system restricted to subspace, if dim < k;
  --   dim        dimension of the subspace.

  function Collapse_Equations
                ( p : Multprec_Complex_Poly_Systems.Poly_Sys; dim : natural32 )
                return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns a system of range 1..dim with random combinations from p.

  function Embed_Collapsed_Equations
               ( restorgsys : Multprec_Complex_Poly_Systems.Poly_Sys;
                 restembsys : Standard_Complex_Poly_Systems.Poly_Sys;
                 dim,level : natural32 )
               return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the embedding of the collapsed system restorgsys (which is
  --   the output of "Collapse_Equations" above).  The last "level" equations
  --   of restembsys (restricted embedded system) are copied to the result.

end Multprec_Subspace_Restrictions;
