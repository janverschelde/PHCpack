with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Matrices;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Matrices;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;

package Rectangular_Sample_Grids is

-- DESCRIPTION :
--   A rectangular sample grid has its samples on parallel slices.
--   This package provides creators, diagnostics and selectors.

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine in the
  --   packages Sample_Points and Sample_Point_Lists.

-- UTILITIES :

  function Extended_Random
             ( v : Standard_Complex_Vectors.Vector; size : natural32 )
             return Multprec_Complex_Vectors.Vector;
  function Extended_Random
             ( v : Standard_Complex_VecVecs.VecVec; size : natural32 )
             return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   The vectors on return are extended with random numbers,
  --   up to the given size.

-- CREATORS :

  function Create1 ( sps : Standard_Sample_List; m : natural32 )
                   return Array_of_Standard_Sample_Lists;
  function Create1 ( sps : Standard_Sample_List; m,size : natural32 )
                   return Array_of_Multprec_Sample_Lists;

  -- DESCRIPTION :
  --   Returns a grid of range 0..m, where every grid(i) contains as
  --   many samples as the list sps on parallel slices.
  --   The constant terms of the hyperplane section are taken to be
  --   roots of unity on the complex unit circle.
  --   For a "square" grid, set m = Length_Of(sps).

  -- ASSUMED : The component is of dimension one.

  -- REQUIRED :
  --   The sampling machine is initialized and properly tuned.
  --   All samples in sps have the same hyperplane sections.

  procedure Create1 ( sps : in Standard_Sample_List; m,size : natural32;
                      stgrid : out Array_of_Standard_Sample_Lists;
                      mpgrid : out Array_of_Multprec_Sample_Lists );

  -- DESCRIPTION :
  --   Returns both the original and refined grid of sample lists,
  --   with the given size of the multi-precision numbers.
  --   The same assumptions and requirments hold as above.

  function Triangular_Create1
             ( sps : Standard_Sample_List; m : natural32 )
             return Array_of_Standard_Sample_Lists;
  function Triangular_Create1
             ( sps : Standard_Sample_List; m,size : natural32 )
             return Array_of_Multprec_Sample_Lists;


  -- DESCRIPTION :
  --   Returns a "triangular" grid of samples as array of range 0..m,
  --   where the i-th list contains d-i+1 samples, for i >= 1 and
  --   with d = Length_Of(sps).

-- DIAGNOSTICS :

  function Errors ( grid : Array_of_Standard_Sample_Lists )
                  return Standard_Floating_Matrices.Matrix;
  function Errors ( grid : Array_of_Multprec_Sample_Lists )
                  return Multprec_Floating_Matrices.Matrix;

  -- DESCRIPTION :
  --   The accuracy of the grid is determined by the errors on the samples.
  --   Rows in the matrix of errors correspond to the lists in the grid.

  function Maximal_Error ( grid_errors : Standard_Floating_Matrices.Matrix )
                         return double_float;
  function Maximal_Error ( grid_errors : Multprec_Floating_Matrices.Matrix )
                         return Floating_Number;
  function Maximal_Error ( grid : Array_of_Standard_Sample_Lists )
                         return double_float;
  function Maximal_Error ( grid : Array_of_Multprec_Sample_Lists )
                         return Floating_Number;

  -- DESCRIPTION :
  --   Returns the maximal errors on the samples in the grid.

  function Distance ( spt1,spt2 : Standard_Sample ) return double_float;
  function Distance ( spt1,spt2 : Multprec_Sample ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the max norm of the difference between the solution vectors
  --   of the two samples.

  function Distance ( sps : Standard_Sample_List; i : natural32; 
                      spt : Standard_Sample ) return double_float;
  function Distance ( sps : Multprec_Sample_List; i : natural32; 
                      spt : Multprec_Sample ) return Floating_Number;

  -- DESCRIPTION :
  --   The distance between the i-th sample spt in the list sps and
  --   all other samples in the list is the minimum of the distances
  --   between that sample and all other samples in the list sps.

  function Distances ( grid : Array_of_Standard_Sample_Lists )
                     return Standard_Floating_Matrices.Matrix;
  function Distances ( grid : Array_of_Multprec_Sample_Lists )
                     return Multprec_Floating_Matrices.Matrix;

  -- DESCRIPTION :
  --   The accuracy of the interpolator is determined by the distance
  --   between the samples in the grid.  The (i,j)-th entry in the matrix
  --   on return lists the minimal distance of the j-th sample in the i-th
  --   list of the grid between all the other samples in that i-th list.
  --   So the number of rows in the matrix on return is grid'length and
  --   the numer of columns equals the number of elements in each list.

  function Minimal_Distance ( grid_dist : Standard_Floating_Matrices.Matrix )
                            return double_float;
  function Minimal_Distance ( grid_dist : Multprec_Floating_Matrices.Matrix )
                            return Floating_Number;
  function Minimal_Distance ( grid : Array_of_Standard_Sample_Lists )
                            return double_float;
  function Minimal_Distance ( grid : Array_of_Multprec_Sample_Lists )
                            return Floating_Number;

  -- DESCRIPTION :
  --   Returns the minimal distance between the samples in the lists
  --   of the grid.

-- SELECTORS :

  function Abscisses
             ( grid : Array_of_Standard_Sample_Lists; i : natural32 )
             return Standard_Complex_Vectors.Vector;
  function Abscisses
             ( grid : Array_of_Multprec_Sample_Lists; i : natural32 )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns -c(0) for the slices in the range 0..i in the grid.
  --   These are the abscisses x_k for the points in the grid:
  --   for every x_k there are d values y(x_k) = y_kl, l=1..d.

  function Extract_Samples
                ( grid : Array_of_Standard_Sample_Lists )
                return Standard_Complex_VecVecs.VecVec;
  function Extract_Samples
                ( grid : Array_of_Multprec_Sample_Lists )
                return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   The vector of sample points consists of a well distributed
  --   selection of points on all slices.
  --   Only the first two coordinates of every solution are returned.

  function Rotate_and_Project ( v,x : Standard_Complex_Vectors.Vector )
                              return Standard_Complex_Vectors.Vector;
  function Rotate_and_Project ( v,x : Multprec_Complex_Vectors.Vector )
                              return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Performs a linear coordinate transformation so that the
  --   first component of the vector on return equals v*x.
  --   The general formula for this linear transformation is
  --     y(1) := v(x'range)*x;
  --     y(i) := v(i)*x(1) - v(1)*x(i), for i in 2..x'last.
  --   Only the first two coordinates of y are needed and returned.

  function Rotate_and_Project2 ( v,x : Standard_Complex_Vectors.Vector )
                               return Standard_Complex_Numbers.Complex_Number;
  function Rotate_and_Project2 ( v,x : Multprec_Complex_Vectors.Vector )
                               return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns only the second coordinate of the rotation, as the first
  --   coordinate in the full Rotate_and_Project is actually minus the
  --   constant term of the slicing hyperplane.

  function Inverse_Rotate ( v,z : Standard_Complex_Vectors.Vector )
                          return Standard_Complex_Vectors.Vector;
  function Inverse_Rotate ( v,z : Multprec_Complex_Vectors.Vector )
                          return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Rotate_and_Project(v,Inverse_Rotate(v,z)) = z.

  function Rotate_Samples ( d,k : natural32;
                            rot : Standard_Complex_Vectors.Vector;
                            grid : Array_of_Standard_Sample_Lists )
                          return Standard_Complex_Matrices.Matrix;
  function Rotate_Samples ( d,k : natural32;
                            rot : Multprec_Complex_Vectors.Vector;
                            grid : Array_of_Multprec_Sample_Lists )
                          return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   The matrix on return contains the 2nd coordinate of the
  --   rotated samples.  If ya is the name of the matrix, ya(i,j)
  --   contains the second coordinate of the i-th sample on
  --   the j-th slice, for i from 1 to d, and j from 0 to d.

  -- REQUIRED : grid'range is 0..d.

  -- ON ENTRY :
  --   d         degree of the grid;
  --   k         index to last slice to be considered;
  --   rot       rotation vector;
  --   grid      grid of samples.

end Rectangular_Sample_Grids;
