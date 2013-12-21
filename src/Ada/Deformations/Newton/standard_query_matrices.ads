with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Complex_Matrices;        use Standard_Complex_Matrices;

package Standard_Query_Matrices is

-- DESCRIPTION :
--   This package offers routine to let the user selective query the elements
--   of a matrix.  For testing operations in large matrices, it is no longer
--   practical to write the whole matrix to screen, instead we use the
--   operations in this package.

  procedure Show_Element ( A : in character; i,j : in integer32;
                           aij : in Complex_Number );

  -- DESCRIPTION :
  --   Shows the (i,j)-th element of the matrix a, in the format
  --   A(i,j) = aij, showing values for i, j, and aij.

  procedure Write_Matrix ( a : in Matrix );

  -- DESCRIPTION :
  --   Writes the elements of the matrix a line by line.

  procedure Show_Matrix ( a : in Matrix );

  -- DESCRIPTION :
  --   Shows the elements in a matrix, allowing the user to pause.

  procedure Query_Matrix ( a : in Matrix );

  -- DESCRIPTION :
  --   Interactive routines to query the elements of a matrix.

  procedure Query_Matrices ( a,b : in Matrix );

  -- DESCRIPTION :
  --   Interactive routines to query the elements of the matrices a and b.

  function Difference_of_Matrices ( a,b : Matrix ) return double_float;

  -- DESCRIPTION :
  --   Returns the absolute value of the sum of all componentwise
  --   differences between corresponding elements in a and b.

  procedure Show_Difference_of_Matrices ( a,b : in Matrix );

  -- DESCRIPTION :
  --   Shows the user the difference between the matrices a and b.

  procedure Show_Differences_in_Matrices ( a,b : in Matrix );

  -- DESCRIPTION :
  --   Only shows the user the different entries between a and b.

end Standard_Query_Matrices;
