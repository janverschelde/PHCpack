with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;

package Facet_Vertex_Enumeration is

-- DESCRIPTION :
--   The following two problems are dual to each other:
--   1) given vertices of a polytope, enumerate all facet inequalities;
--   2) given facet inequalities of a polytope, enumerate all vertices.
--   This package provides data encapsulation routines.

-- DATA MANIPULATION :

  function List_to_Matrix
               ( n : integer32; L : Lists_of_Floating_Vectors.List )
               return Standard_Floating_Matrices.Matrix;

  -- DESCRIPTION :
  --   The columns of the matrix on return contain the first n entries
  --   of the vectors in the given list L.

  function List_to_Vector
               ( i : integer32; L : Lists_of_Floating_Vectors.List )
               return Standard_Floating_Vectors.Vector;

  -- DESCRPTION :
  --   The vector on return contains the ith entry of every vector
  --   in the given list L.

  function Matrix_to_List
               ( n : integer32; m : Standard_Floating_Matrices.Matrix )
               return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION :
  --   The list on return contains as vectors the first n entries of
  --   the columns of the given matrix m.

-- ENUMERATORS OF VERTICES :

  procedure Enumerate_Vertices
               ( cff : in Standard_Floating_Matrices.Matrix;
                 rhs : in Standard_Floating_Vectors.Vector;
                 tol : in double_float;
                 points : out Lists_of_Floating_Vectors.List;
                 labels : out Lists_of_Integer_Vectors.List;
                 fail,infty : out boolean );

  -- DESCRIPTION :
  --   Returns the list of all vertex points x satisfying the inequalities:
  --     cff(1,i)*x(1) + cff(2,i)*x(2) + .. + cff(n,i)*x(n) >= rhs(i),
  --   for i running along all columns of the coefficient matrix cff,
  --   and where "=" is reached for at least n columns of cff.
  --   Also the labels of the equations where "=" is returned.

  -- REQUIRED : cff'range(2) = rhs'range.

  -- ON ENTRY :
  --   cff       coefficients of the inequalities;
  --   rhs       right-hand size vector;
  --   tol       tolerance to decide whether a number equals zero or not.

  -- ON RETURN :
  --   points    coordinates of the vertex points;
  --   labels    labels to the defining equalities of the vertex points;
  --   fail      the system of inequalities is infeasible;
  --   infty     optimization problem is unbounded.

  function Enumerate_Vertex_Points
               ( cff : Standard_Floating_Matrices.Matrix;
                 rhs : Standard_Floating_Vectors.Vector; tol : double_float )
               return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of all vertex points x satisfying the inequalities:
  --     cff(1,i)*x(1) + cff(2,i)*x(2) + .. + cff(n,i)*x(n) >= rhs(i),
  --   for i running along all columns of the coefficient matrix cff,
  --   and where "=" is reached for at least n columns of cff.

  -- REQUIRED : cff'range(2) = rhs'range.

  -- ON ENTRY :
  --   cff       coefficients of the inequalities;
  --   rhs       right-hand side vector;
  --   tol       tolerance to decide whether a number equals zero or not.

  function Enumerate_Vertex_Labels
               ( cff : Standard_Floating_Matrices.Matrix;
                 rhs : Standard_Floating_Vectors.Vector; tol : double_float )
               return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of all vertex points x satisfying the inequalities:
  --     cff(1,i)*x(1) + cff(2,i)*x(2) + .. + cff(n,i)*x(n) >= rhs(i),
  --   for i running along all columns of the coefficient matrix cff,
  --   and where "=" is reached for at least n columns of cff.
  --   The list on return contains the labels of the equalities that
  --   define the vertex points.

  -- REQUIRED : cff'range(2) = rhs'range.

  -- ON ENTRY :
  --   cff       coefficients of the inequalities;
  --   rhs       right-hand size vector;
  --   tol       tolerance to decide whether a number equals zero or not.

-- ENUMERATORS OF FACETS :

  procedure Enumerate_Facets
               ( pts : in Standard_Floating_Matrices.Matrix;
                 tol : in double_float;
                 facets : out Lists_of_Floating_Vectors.List;
                 labels : out Lists_of_Integer_Vectors.List;
                 fail,infty : out boolean );

  -- DESCRIPTION :
  --   Returns the list of facet inequalities x = (x(1..n),x(n+1)):
  --     pts(1,i)*x(1) + pts(2,i)*x(2) + .. + pts(n,i)*x(n) >= x(n+1),
  --   for all i running along all columns of the matrix of vertices,
  --   and where "=" is reached for at least n vertex points.

  -- ON ENTRY :
  --   pts       columns of the matrix contain points that span polytope;
  --   tol       tolerance to decide whether a number equals zero or not.

  -- ON RETURN :
  --   facets    facet inequalities, in the format as x above;
  --   labels    labels to the vertex points that span the facet;
  --   fail      the system of inequalities is infeasible;
  --   infty     optimization problem is unbounded.

  function Enumerate_Facet_Inequalities
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of facet inequalities x = (x(1..n),x(n+1)):
  --     pts(1,i)*x(1) + pts(2,i)*x(2) + .. + pts(n,i)*x(n) >= x(n+1),
  --   for all i running along all columns of the matrix of vertices,
  --   and where "=" is reached for at least n vertex points.

  function Enumerate_Facet_Labels
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Computes the list of facet inequalities x = (x(1..n),x(n+1)):
  --     pts(1,i)*x(1) + pts(2,i)*x(2) + .. + pts(n,i)*x(n) >= x(n+1),
  --   for all i running along all columns of the matrix of vertices,
  --   and where "=" is reached for at least n vertex points.
  --   The list on return contains the labels of the points that span
  --   the facets.

-- ENUMERATORS OF LOWER FACETS :

  procedure Enumerate_Lower_Facets
               ( pts : in Standard_Floating_Matrices.Matrix;
                 tol : in double_float;
                 facets : out Lists_of_Floating_Vectors.List;
                 labels : out Lists_of_Integer_Vectors.List;
                 fail,infty : out boolean );

  -- DESCRIPTION :
  --   Returns the facet inequalities of all lower facets, 
  --   where x(n) = 1 for all vectors x in the list on return.
  --   The parameters have the same meaning as in the facet
  --   enumeration analogue above.

  function Enumerate_Lower_Facet_Inequalities
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION :
  --   Returns the facet inequalities of all lower facets, 
  --   where x(n) = 1 for all vectors x in the list on return.
  --   The format of the input is the facet enumeration analogue above.

  function Enumerate_Lower_Facet_Labels
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns the facet inequalities of all lower facets, 
  --   where x(n) = 1 for all vectors x in the list on return.
  --   The format of the input is the facet enumeration analogue above.

end Facet_Vertex_Enumeration;
