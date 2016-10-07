with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;

package Integer_Lifting_Utilities is

-- DESCRIPTION :
--   This package provides some utilities for dealing with lifting functions.

  function Adaptive_Lifting ( L : Array_of_Lists ) return Vector;

  -- DESCRIPTION :
  --   Returns upper bounds for a random lifting, depending on the lengths
  --   of the lists in L.

--  function Select_Subsystem ( p : Laur_Sys; mix : Vector; mic : Mixed_Cell )
--                            return Laur_Sys;

  -- DESCRIPTION :
  --   Given a Laurent polynomial system and a mixed cell,
  --   the corresponding subsystem will be returned.

  -- ON ENTRY :
  --   p          a Laurent polynomial system;
  --   mix        type of mixture: occurencies of the supports;
  --   mic        a mixed cell.

  -- REQUIRED :
  --   The polynomials in p must be ordered according to the type of mixture.

  function Perform_Lifting
               ( n : integer32; L : List; p : Poly ) return Poly;
  function Perform_Lifting
               ( n : integer32; mix : Vector; lifted : Array_of_Lists;
                 p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Construction of the polyhedral homotopy, given the lifted supports.
  --   The p will be lifted according to the lifted points, i.e. each
  --   monomial whose exponent vector corresponds to a vector in the
  --   lifted points (except for the last entry of course), will be
  --   extended with an additional unknown, with as exponent, the
  --   lifting value of that vector.

  function Copy_Lifting ( lifted : List; pt : Link_to_Vector )
                        return Link_to_Vector;

  -- DESCRIPTION :
  --   Searches the corresponding point in the list lifted and returns
  --   the lifted point.  If the corresponding point has not been found,
  --   then the original point pt will be returned.

  function Copy_Lifting ( lifted,pts : List ) return List;

  -- DESCRIPTION :
  --   Copies the lifting on the points lifted to the points in pts,
  --   i.e., each point in pts will get the same lifting as the corresponding
  --   lifted point in the list lifted.

  procedure Search_Lifting ( L : in List; pt : in Vector;
                             found : out boolean; lif : out integer32 );

  -- DESCRIPTION :
  --   Searches the lifting of the point in the lifted list l.
  --   If found, then lif equals the lifting, otherwise lif is meaningless.

  function Search_and_Lift ( L : List; pt : Vector ) return Vector;

  -- DESCRIPTION :
  --   Given a lifted list of points and a unlifted vector, the function
  --   either returns the corresponding lifted vector from the list, or
  --   the same point, when there is no lifted point in l whose projection
  --   equals the given point pt.

  function Induced_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                             pt : Vector ) return Vector;
  function Induced_Lifting
               ( n : integer32; mix : Vector; points : Array_of_Lists;
                 mixsub : Mixed_Subdivision ) return Array_of_Lists;

  -- DESCRIPTION :
  --   Given a mixed subdivision for a tuple of supports,
  --   then the lifted points will be returned as induced by the
  --   subdivision.   When points do not occur in the mixed subdivision,
  --   they will be lifted conservatively.

  procedure Constant_Lifting
               ( L : in List; liftval : in integer32;
                 lifted,lifted_last : in out List );

  -- DESCRIPTION :
  --   Gives all points in l the constant lifting liftval,
  --   and appends them to the list lifted, where lifted_last is
  --   a pointer to the last element of lifted.

  procedure Constant_Lifting
               ( al : in Array_of_Lists; liftval : in integer32;
                 lifted,lifted_last : in out Array_of_Lists );

  -- DESCRIPTION :
  --   Gives all points in al(i) the constant lifting liftval,
  --   and appends them to the lists lifted(i), where lifted_last(i) is
  --   a pointer to the last element of lifted(i).

  function Conservative_Lifting 
               ( mic : Mixed_Cell; k : integer32; point : Vector )
               return integer32;
  function Conservative_Lifting 
               ( mixsub : Mixed_Subdivision; k : integer32; point : Vector )
               return integer32;
  
  -- DESCRIPTION :
  --   Returns the value of the conservative lifting function of the point
  --   to be considered w.r.t. the kth polytope.

  -- REQUIRED :
  --   The given point must already be in the lifted space and its last
  --   coordinate must contain already a lower bound for the lifting value.

  function Lower_Lifting ( mic : Mixed_Cell; k : integer32; point : Vector )
                         return integer32;

  -- DESCRIPTION :
  --   Returns a lower bound on the lifting value.  In case that the point
  --   belongs to the kth component of the cell, the value zero is returned.
  --   Otherwise, the conservative lifting value of the point w.r.t. the
  --   cell is returned.

  -- REQUIRED :
  --   The point must be already in lifted space.

  function Lower_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                           point : Vector ) return integer32;

  -- DESCRIPTION :
  --   Applies the lower lifting function to all cells in the subdivision.
  --   Stops when the lower bound equals the lifting value of the point
  --   equals the lifting of the point given on entry.

  function Lifted_Supports ( r : integer32; mixsub : Mixed_Subdivision )
                           return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns an array of lists of range 1..r with the lifted supports
  --   for all points in the given mixed subdivision.

  -- REQUIRED : r = c.pts'last for all cells c in mixsub.

end Integer_Lifting_Utilities;
