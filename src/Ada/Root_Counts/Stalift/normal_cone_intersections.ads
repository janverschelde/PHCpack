with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;

package Normal_Cone_Intersections is

-- DESCRIPTION :
--   This package provides a data abstraction to represent the intersections
--   of the generators of a normal cone with a tuple of normal cone complexes.

-- DATA STRUCTURE :

  type Intersection_Matrix ( n,m,nc : integer32 ) is record
    sv : Vector(1..n);
    im : Matrix(0..m,1..nc);
  end record; 

  -- The aim of an intersection matrix is to answer the question :
  --   does the ith generator of the normal cone of the point x
  --   belong to the normal cone of the kth point of the jth support?

  -- The parameters of the three-dimensional intersection matrix are
  --   n  : number of supports to consider the point x to;
  --   m  : number of generators of the normal cone to x;
  --   nc : total number of normal cones that need to be considered.

  -- The data of the intersection matrix are
  --   sv(1..n) a vector whose entries indicates the starting position of the
  --     the supports in the matrix im, i.e., sv(i) gives the first column in
  --     the matrix im that collects the data of the normal cone of the first
  --     points in the ith support;
  --   im(0..m,1..nc) is a matrix that contains the answers to the question:
  --     im(i,sv(j)+k) equals 0 or 1, 0 when the ith generator does not belong
  --     to the normal cone of the kth points in the jth support, 1 otherwise;
  --     im(0,sv(j)+k) is the sum of im(i,sv(j)+k) for all i in 1..m.

-- CONSTRUCTORS :

  function Number_of_Cones
             ( L : Array_of_Lists; i : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Computes the number of cones, i.e.: returns the sum of the lengths of 
  --   all lists of l, without consider the ith one.
  --   This is an auxiliary for determining the third dimension of the 
  --   intersection matrix.

  function Lengths ( L : Array_of_Lists; i : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of dimension equal to l that accumulates the lengths
  --   of the lists of l, without considering the ith component.
  --   More precisely, the jth component of the vector on return equals
  --   one plus the sum of all lengths of the lists in l(l'first..j),
  --   minus the length of l(i), if i < j.  So the vector on return can serve
  --   as the vector sv, except for the last component that equals one plus
  --   the total number of cones to consider.

  function Create ( L : Array_of_Lists; g : List; i : integer32 ) 
                  return Intersection_Matrix;

  -- DESCRIPTION :
  --   Returns the intersection matrix of the list of generators of the normal
  --   cone of a point that belongs to the list l(i).

-- ELEMENTARY SELECTORS :

  function Is_In ( ima : Intersection_Matrix;
                   i,j,k : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true when the ith generator belongs to the normal cone of the
  --   kth point of the jth support list.

  function Maximal_Column ( ima : Intersection_Matrix ) return integer32;

  -- DESCRIPTION :
  --   Returns the index to the column in the intersection matrix with
  --   the maximal column sum.

  function Component ( ima : Intersection_Matrix; column : integer32 ) 
                     return integer32;

  -- DESCRIPTION :
  --   Returns the number of the component of the intersection matrix the
  --   given column index corresponds to.  The number on return equals the
  --   index of the support of the corresponding column.

  function Length ( ima : Intersection_Matrix;
                    i : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the length of the ith component of the matrix, i.e., returns
  --   the length of the ith support list.

  function Row_Sum ( ima : Intersection_Matrix;
                     i,j : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the sum of the elements on the ith row, for all normal cones
  --   of the points of the jth support.

-- ENUMERATING COMPLEMENTARY COLUMS :

  generic

    with procedure Process ( cols : in Vector; continue : out boolean );

    -- DESCRIPTION :
    --   This procedure is invoked each time a set of complementary columns
    --   has been found.  The vector on return has the following meaning:
    --     cols(i) = 0 means that no normal cone of the ith component is taken,
    --     cols(i) = j means that the jth normal cone of the ith component
    --                 belongs to the complementary columns.
    --   Note that the range of cols is 1..n-1, with n = #supports.

  procedure Complementary_Columns ( ima : in Intersection_Matrix );

  -- DESCRIPTION :
  --   This procedure enumerates all complementary columns.
  --   A set of columns, at most one of each component, is said to be
  --   complementary if its union contains all generators of a normal cone.

  function Partition ( ima : Intersection_Matrix; cols : Vector; g : List )
                     return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns a partition of the set of generators g w.r.t. the columns cols
  --   in the intersection matrix ima.  More precisely: the ith list on return
  --   contains those generators that belong to the normal cone cols(i) of the
  --   ith component.
  --   If the same generator belongs to several cones, it will be contained
  --   only in the list with smallest index.

  -- REQUIRED : cols is a set of complementary columns.

  function Partition_in_Union ( partg,points : Array_of_Lists; i : integer32;
                                cols : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the set of generators of a normal cones of the ith
  --   component of the point lists belongs to the union of normal cones,
  --   as indicated by the set of complementary columns.

  function Contained_in_Union
             ( l : Array_of_Lists; i : integer32; g : List;
               ima : Intersection_Matrix; cols : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the list of generators g is contained in the normal
  --   cones to the point selected by the columns of the intersection matrix.

-- FINAL TARGET ROUTINE :

  function Contained_in_Union
             ( l : Array_of_Lists; i : integer32; g : List;
               ima : Intersection_Matrix ) return boolean;

  -- DESCRIPTION :
  --   Enumerates all sets of complementary columns, until its union
  --   contains the convex cone spanned by the list of generators,
  --   until all possibilities are exhausted.
  --   In the first case, true is returned, otherwise false is returned.

end Normal_Cone_Intersections;
