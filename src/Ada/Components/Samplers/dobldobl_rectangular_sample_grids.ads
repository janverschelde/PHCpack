with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with Double_Double_Matrices;
with DoblDobl_Complex_Matrices;
with DoblDobl_Sample_Points;             use DoblDobl_Sample_Points;
with DoblDobl_Sample_Lists;              use DoblDobl_Sample_Lists;

package DoblDobl_Rectangular_Sample_Grids is

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

-- CREATORS :

  function Create1 ( sps : DoblDobl_Sample_List; m : natural32 )
                   return Array_of_DoblDobl_Sample_Lists;

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

  function Triangular_Create1
             ( sps : DoblDobl_Sample_List; m : natural32 )
             return Array_of_DoblDobl_Sample_Lists;

  -- DESCRIPTION :
  --   Returns a "triangular" grid of samples as array of range 0..m,
  --   where the i-th list contains d-i+1 samples, for i >= 1 and
  --   with d = Length_Of(sps).

-- DIAGNOSTICS :

  function Errors ( grid : Array_of_DoblDobl_Sample_Lists )
                  return Double_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   The accuracy of the grid is determined by the errors on the samples.
  --   Rows in the matrix of errors correspond to the lists in the grid.

  function Maximal_Error ( grid_errors : Double_Double_Matrices.Matrix )
                         return double_double;
  function Maximal_Error ( grid : Array_of_DoblDobl_Sample_Lists )
                         return double_double;

  -- DESCRIPTION :
  --   Returns the maximal errors on the samples in the grid.

  function Distance ( spt1,spt2 : DoblDobl_Sample ) return double_double;

  -- DESCRIPTION :
  --   Returns the max norm of the difference between the solution vectors
  --   of the two samples.

  function Distance ( sps : DoblDobl_Sample_List; i : natural32; 
                      spt : DoblDobl_Sample ) return double_double;

  -- DESCRIPTION :
  --   The distance between the i-th sample spt in the list sps and
  --   all other samples in the list is the minimum of the distances
  --   between that sample and all other samples in the list sps.

  function Distances ( grid : Array_of_DoblDobl_Sample_Lists )
                     return Double_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   The accuracy of the interpolator is determined by the distance
  --   between the samples in the grid.  The (i,j)-th entry in the matrix
  --   on return lists the minimal distance of the j-th sample in the i-th
  --   list of the grid between all the other samples in that i-th list.
  --   So the number of rows in the matrix on return is grid'length and
  --   the numer of columns equals the number of elements in each list.

  function Minimal_Distance ( grid_dist : Double_Double_Matrices.Matrix )
                            return double_double;
  function Minimal_Distance ( grid : Array_of_DoblDobl_Sample_Lists )
                            return double_double;

  -- DESCRIPTION :
  --   Returns the minimal distance between the samples in the lists
  --   of the grid.

-- SELECTORS :

  function Abscisses
             ( grid : Array_of_DoblDobl_Sample_Lists; i : natural32 )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns -c(0) for the slices in the range 0..i in the grid.
  --   These are the abscisses x_k for the points in the grid:
  --   for every x_k there are d values y(x_k) = y_kl, l=1..d.

  function Extract_Samples
                ( grid : Array_of_DoblDobl_Sample_Lists )
                return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   The vector of sample points consists of a well distributed
  --   selection of points on all slices.
  --   Only the first two coordinates of every solution are returned.

  function Rotate_and_Project ( v,x : DoblDobl_Complex_Vectors.Vector )
                              return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Performs a linear coordinate transformation so that the
  --   first component of the vector on return equals v*x.
  --   The general formula for this linear transformation is
  --     y(1) := v(x'range)*x;
  --     y(i) := v(i)*x(1) - v(1)*x(i), for i in 2..x'last.
  --   Only the first two coordinates of y are needed and returned.

  function Rotate_and_Project2 ( v,x : DoblDobl_Complex_Vectors.Vector )
                               return DoblDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns only the second coordinate of the rotation, as the first
  --   coordinate in the full Rotate_and_Project is actually minus the
  --   constant term of the slicing hyperplane.

  function Inverse_Rotate ( v,z : DoblDobl_Complex_Vectors.Vector )
                          return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Rotate_and_Project(v,Inverse_Rotate(v,z)) = z.

  function Rotate_Samples ( d,k : natural32;
                            rot : DoblDobl_Complex_Vectors.Vector;
                            grid : Array_of_DoblDobl_Sample_Lists )
                          return DoblDobl_Complex_Matrices.Matrix;

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

end DoblDobl_Rectangular_Sample_Grids;
