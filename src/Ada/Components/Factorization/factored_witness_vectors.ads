with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;

package Factored_Witness_Vectors is

-- DESCRIPTION :
--   The factors of a multivariate polynomial are represented by
--   vectors of generic points, representing witness sets.
--   The functions and procedures defined in this package manipulate
--   these vectors of witness points.
--   Three different levels of precision are supported:
--   standard double, double double, and quad double precision.

  procedure Swap ( m : in out Standard_Natural_Vectors.Vector;
                   i,j : in integer32 );

  -- DESCRIPTION :
  --   Swaps the elements m(i) and m(j).

  procedure Swap ( v : in out Standard_Complex_Vectors.Vector;
                   i,j : in integer32 );
  procedure Swap ( v : in out DoblDobl_Complex_Vectors.Vector;
                   i,j : in integer32 );
  procedure Swap ( v : in out QuadDobl_Complex_Vectors.Vector;
                   i,j : in integer32 );

  -- DESCRIPTION :
  --   Swaps the elements v(i) and v(j).

  procedure Sort ( m : in out Standard_Natural_Vectors.Vector;
                   w : in out Standard_Complex_Vectors.Vector );
  procedure Sort ( m : in out Standard_Natural_Vectors.Vector;
                   w : in out DoblDobl_Complex_Vectors.Vector );
  procedure Sort ( m : in out Standard_Natural_Vectors.Vector;
                   w : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sorts the generic points in w according to the multiplicities in m
  --   in ascending order.

  function Countmu ( m : Standard_Natural_Vectors.Vector; mu : natural32 )
                   return natural32;

  -- DESCRIPTION :
  --   Returns the number of entries i in m for which m(i) = mu.

  function Select_Multiple_Factors
               ( m : Standard_Natural_Vectors.Vector;
                 w : Standard_Complex_Vectors.Vector; mu : natural32 ) 
               return Standard_Complex_Vectors.Vector;
  function Select_Multiple_Factors
               ( m : Standard_Natural_Vectors.Vector;
                 w : DoblDobl_Complex_Vectors.Vector; mu : natural32 ) 
               return DoblDobl_Complex_Vectors.Vector;
  function Select_Multiple_Factors
               ( m : Standard_Natural_Vectors.Vector;
                 w : QuadDobl_Complex_Vectors.Vector; mu : natural32 ) 
               return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Selects those witness points in w with multiplicity mu.

  -- ASSUMED :
  --   There is at least one witness point of multiplicity mu.

  function Is_In ( v : Standard_Complex_Vectors.Vector;
                   x : Standard_Complex_Numbers.Complex_Number;
                   tol : double_float ) return boolean;
  function Is_In ( v : DoblDobl_Complex_Vectors.Vector;
                   x : DoblDobl_Complex_Numbers.Complex_Number;
                   tol : double_float ) return boolean;
  function Is_In ( v : QuadDobl_Complex_Vectors.Vector;
                   x : QuadDobl_Complex_Numbers.Complex_Number;
                   tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there is an entry i in v such that |v(i) - x| <= tol.

  function Position ( v : Standard_Complex_Vectors.Vector;
                      x : Standard_Complex_Numbers.Complex_Number;
                      tol : double_float ) return integer32;
  function Position ( v : DoblDobl_Complex_Vectors.Vector;
                      x : DoblDobl_Complex_Numbers.Complex_Number;
                      tol : double_float ) return integer32;
  function Position ( v : QuadDobl_Complex_Vectors.Vector;
                      x : QuadDobl_Complex_Numbers.Complex_Number;
                      tol : double_float ) return integer32;

  -- DESCRIPTION :
  --   Returns the position of the point x in the vector v.

  function Positions ( v,x : Standard_Complex_Vectors.Vector;
                       tol : double_float )
                     return Standard_Natural_Vectors.Vector;
  function Positions ( v,x : DoblDobl_Complex_Vectors.Vector;
                       tol : double_float )
                     return Standard_Natural_Vectors.Vector;
  function Positions ( v,x : QuadDobl_Complex_Vectors.Vector;
                       tol : double_float )
                     return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of positions of the elements of x in v.

  function Remove_Duplicates
               ( w : Standard_Complex_Vectors.Vector; tol : double_float )
               return Standard_Complex_Vectors.Vector;
  function Remove_Duplicates
               ( w : DoblDobl_Complex_Vectors.Vector; tol : double_float )
               return DoblDobl_Complex_Vectors.Vector;
  function Remove_Duplicates
               ( w : QuadDobl_Complex_Vectors.Vector; tol : double_float )
               return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   The vector on return does not contain elements that are as close
  --   to other elements as indicated by the given tolerance tol.

  function Remove_Duplicates
               ( w : Standard_Complex_Vectors.Vector; tol : double_float;
                 m : Standard_Natural_Vectors.Vector )
               return Standard_Natural_Vectors.Vector;
  function Remove_Duplicates
               ( w : DoblDobl_Complex_Vectors.Vector; tol : double_float;
                 m : Standard_Natural_Vectors.Vector )
               return Standard_Natural_Vectors.Vector;
  function Remove_Duplicates
               ( w : QuadDobl_Complex_Vectors.Vector; tol : double_float;
                 m : Standard_Natural_Vectors.Vector )
               return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Only the multiplicity of the first occurrence of a multiple root
  --   will be retained.

  -- ON ENTRY :
  --   w         witness points, with multiple occurrences;
  --   tol       tolerance to decide whether two numbers are equal;
  --   m         multiplicities of the witness points in w.

  function Select_Points ( w : in Standard_Complex_Vectors.Vector;
                           k : in Standard_Natural_Vectors.Vector )
                         return Standard_Complex_Vectors.Vector;
  function Select_Points ( w : in DoblDobl_Complex_Vectors.Vector;
                           k : in Standard_Natural_Vectors.Vector )
                         return DoblDobl_Complex_Vectors.Vector;
  function Select_Points ( w : in QuadDobl_Complex_Vectors.Vector;
                           k : in Standard_Natural_Vectors.Vector )
                         return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns those points in w whose entries are in k.

  function Multiplicity_of_Factors
              ( factors : Standard_Natural_VecVecs.VecVec;
                m : in Standard_Natural_Vectors.Vector )
              return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the multiplicity of each factor.

end Factored_Witness_Vectors;
