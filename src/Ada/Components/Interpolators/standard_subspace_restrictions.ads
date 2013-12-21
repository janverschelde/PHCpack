with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

package Standard_Subspace_Restrictions is

-- DESCRIPTION :
--   This package provides constructors for the equations of the container
--   subspaces and utilities to restrict solutions to these subspaces.

  procedure Container_Dimension
                ( file : in file_type; k,n : in natural32;
                  samples : in Standard_Complex_VecVecs.VecVec;
                  tol : in double_float; mat,vec : out Matrix;
                  dim : out natural32 );

  -- DESCRIPTION :
  --   Determines the dimension of the space spanned by the samples.

  -- ON ENTRY :
  --   file       to write diagnostics;
  --   k          the number of vectors used in the vector configuration;
  --   n          dimension of the original space, before embedding;
  --   samples    array of at least k+1 points;
  --   tol        tolerance to determine whether number is zero or not;

  -- ON RETURN :
  --   mat        matrix that contains points in its rows, range 1..k+1,1..n;
  --   vec        triangulated vector configuration, range 1..k,1..n;
  --   dim        dimension of the space spanned by the samples.

  procedure Container_Subspace
                ( file : in file_type; k,n,level,dim : in natural32;
                  p : in Poly_Sys; mat,vec : in Matrix; tol : in double_float;
                  samples : in Standard_Complex_VecVecs.VecVec;
                  restsamp : out Standard_Complex_VecVecs.VecVec;
                  lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                  kerpols,restp : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Constructs equations for the linear subspace spanned by points
  --   and restricts the polynomial system to this subspace.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   k          number of vectors used in vec, which implies that the
  --              useful number of rows in vec is 1..k, and the
  --              useful number of rows in mat is 1..k+1;
  --   n          original number of variables;
  --   level      number of added variables in the embedding;
  --   dim        dimension of the container subspace;
  --   p          embedded polynomial system of dimension = n+level;
  --   mat        point configuration that spans the subspace;
  --   vec        vector configuration in upper triangular form;
  --   tol        tolerance to decide whether number is zero;
  --   samples    current set of samples.

  -- ON RETURN :
  --   restsamp   are the samples restricted to the subspace;
  --   lpiv       link to pivot vector;
  --   kerpols    linear polynomials describing the subspace;
  --   restp      polynomial system p restricted to the subspace.

  procedure Subspace_Restriction
                ( file : in file_type; embsys : in Poly_Sys;
                  ind,k,n,level : in natural32;
                  samples : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  restsamp : in out Standard_Complex_VecVecs.Array_of_VecVecs;
                  lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                  kerpols,restembsys : out Link_to_Poly_Sys;
                  dim : out natural32 );

  -- DESCRIPTION :
  --   Determines the dimension of the subspace spanned by the samples
  --   and restricts the polynomial system to this subspace if dim < n.

  -- ON ENTRY :
  --   file       output file to write diagnostics and intermediate results;
  --   embsys     embedded polynomial system;
  --   ind        index to the current solution samples;
  --   k          size of vector configuration;
  --   n          dimension of the original space before the embedding;
  --   level      number of slices added in the embedding;
  --   samples    samples(ind) must contain at least k+1 points.

  -- ON RETURN :
  --   restsamp   samples restricted to the subspace;
  --   lpiv       link to the pivots, indicating the remaining variables;
  --   kerpols    defining polynomials of the subspace, if dim < k;
  --   restembsys is the embedded system restricted to subspace, if dim < k;
  --   dim        dimension of the subspace.

  function Collapse_Equations
                ( p : Poly_Sys; dim,level : natural32 ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns a random combination equations from p to get a square system.
  --   All nonzero equations in p have as many variables as dim.
  --   The last "level" equations are copied.

end Standard_Subspace_Restrictions;
