with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Lists_of_Integer_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package m_homogeneous_permanent_factors is

-- DESCRIPTION :
--   Solving m-homogeneous linear-product start systems mirrors 
--   the row expansion for the permanent of a degree matrix.
--   The factors that contribute to the permanent identify the
--   linear systems that define the start solutions.

  procedure Permanent
              ( row : in integer32;
                deg : in Standard_Integer_Matrices.Matrix;
                cols,crd : in out Standard_Integer_Vectors.Vector;
                per : in out integer32;
                ind,ind_last : in out Lists_of_integer_Vectors.List );

  -- DESCRIPTION :
  --   Row expansion for the permanent of an integer matrix.

  -- ON ENTRY :
  --   row      current row index, initialize with deg'first(1);
  --   deg      degree matrix, with as many columns as crd'length;
  --   cols     selected column indices of the degree matrix,
  --            the range of cols must be deg'range(2);
  --   crd      cardinalities of the sets in the partition,
  --            crd(k) equals the number of elements to be taken
  --            of the k-th column of deg;
  --   per      initialize with zero.

  -- ON RETURN :
  --   per      the permanent of the matrix deg;
  --   ind      list of all column indices used in the permanent;
  --   ind_last is the pointer the to last element in the list ind.

  procedure Split_Indices
              ( ind : integer32;
                accu : in out Standard_Integer_Vectors.Vector;
                deg : in Standard_Integer_Matrices.Matrix;
                cols,base : in Standard_Integer_Vectors.Vector;
                res,res_last : in out Lists_of_Integer_Vectors.List );

  -- DESCRIPTION :
  --   Splits the column indices in cols with respect to the base,
  --   appending to the list res with its pointer to last rest_last.
  --   Initialize ind with accu'first, pointing to the first element
  --   in the accumulator.

  function Split_Column_Indices
             ( deg : Standard_Integer_Matrices.Matrix;
               ind : Lists_of_Integer_Vectors.List )
             return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   For every vector of column indices to the degree matrix deg,
  --   returns a vector of indices to the equations in the corresponding
  --   linear-product system.  In particular, for every vector v in the
  --   list on return we have j = v(i) equals the j-th linear factor in
  --   the i-th equation of the linear-product start system.

  procedure Permanent_Factors
              ( p : in Poly_Sys; z : in Partition;
                sols : out Lists_of_Integer_Vectors.List );

  -- DESCRIPTION :
  --   Computes the permanent of the degree table for the partition z
  --   of the set of unknowns of the system p.
  --   Returns in sols those indices of the factors in the linear equations
  --   of the linear-product start system with respect to the partition z.
  --   The number of elements in sols equals the m-homogeneous Bezout number
  --   with respect to the partition z.

  function Solve_Linear_System
             ( ind : Standard_Integer_Vectors.Vector ) return Solution;

  -- DESCRIPTION :
  --   Solves the system defined by the linear factors where
  --   the factor of the k-th equation is defined by ind(k),
  --   for k in ind'range.

  -- REQUIRED : Standard_Linear_Product_System is initialized well.

  procedure Solve_m_Homogeneous_Start_System
              ( ind_sols : in Lists_of_Integer_Vectors.List;
                q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Solves the linear-product start system stored in 
  --   Standard_LInear_Product_System, using the list of indices.

end m_homogeneous_permanent_factors;
