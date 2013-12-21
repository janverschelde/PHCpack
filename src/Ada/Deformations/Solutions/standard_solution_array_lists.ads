with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Generic_Lists;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_Solution_Array_Lists is

-- DESCRIPTION :
--   When path tracking we convert a solution list into an array.
--   For huge solution lists this causes indexing problems.
--   Therefore we use lists of solution arrays.

  type Link_to_Solution_Array is access Solution_Array;
  package Lists_of_Solution_Arrays is
    new Generic_Lists(Link_to_Solution_Array);
  type Solution_Array_List is new Lists_of_Solution_Arrays.List;

-- CREATORS :

  function Create ( s : Solution_Array ) return Solution_Array_List;
  function Create ( s : Link_to_Solution_Array ) return Solution_Array_List;

  -- DESCRIPTION :
  --   Returns a list with one element: the given array s.

  function Create ( s : Solution_List; size : natural32 )
                  return Solution_Array_List;

  -- DESCRIPTION :
  --   Returns a list of solution arrays where every array contains
  --   as many solutions as the given size, except perhaps for the 
  --   last array in the list which may contain fewer.

  procedure Append ( first,last : in out Solution_Array_List;
                     sa : in Solution_Array );

  -- DESCRIPTION :
  --   Appends the solution array to the list of solution arrays,
  --   with head in first and last element pointed to by last.

-- CONVERTOR : 

  function Concat ( s : Solution_Array_List ) return Solution_List;

  -- DESCRIPTION :
  --   Concatenates the solution arrays in s into one big list.
  --   The pointers of the solutions are copied.

-- DESTRUCTORS :

  procedure Deep_Clear ( s : in out Link_to_Solution_Array );
  procedure Shallow_Clear ( s : in out Link_to_Solution_Array );
  procedure Deep_Clear ( s : in out Solution_Array_List );
  procedure Shallow_Clear ( s : in out Solution_Array_List );

  -- DESCRIPTION :
  --   A shallow clear releases only the pointers, while a deep clear
  --   destroys also all content and nested pointer structures.

end Standard_Solution_Array_Lists;
