with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Sets_of_Unknowns is

-- DESCRIPTION :
--   This package provides a data abstraction for dealing with
--   sets with a limited number of unknowns.

  type Set is private;

-- CREATORS :

  function Create ( n : natural32 ) return Set;
  
  -- DESCRIPTION :
  --   Creates a set which can contain n unknowns.
  --   After this operation Dimension(s) = n.

  function Create ( s : Set ) return Set;

  -- DESCRIPTION :
  --   Returns a new set with the same contents as the given one.
  --   Note that s1 := s2 does only create a new name of the same set,
  --   but does not create a new set.

  function Universe ( n : natural32 ) return Set;

  -- DESCRIPTION : Returns the set of all n unknowns.

-- CONSTRUCTORS :

  procedure Add ( s : in out Set; i : in natural32 );

  -- DESCRIPTION :
  --   Adds the ith unknown to the set.
  --   The set must be created and 1 <= i <= Dimension(s).

  procedure Union ( s1 : in out Set; s2 : in Set );
  function  Union ( s1,s2 : Set ) return Set;

  -- DESCRIPTION :
  --   Constructs the union of two sets of equal dimension.
  --   Either the result will be contained in s1, or a new set
  --   will be constructed.

  procedure Remove ( s : in out Set; i : in natural32 );

  -- DESCRIPTION :
  --   Removes the ith unknown from the set.
  --   The set must be created and 1 <= i <= Dimension(s);

  procedure Difference ( s1 : in out Set; s2 : in Set );
  function  Difference ( s1,s2 : Set ) return Set;

  -- DESCRIPTION :
  --   Constructs the difference s1\s2 of two sets of equal dimension.
  --   Either the result will be contained in s1, or a new set
  --   will be constructed.

  procedure Intersection ( s1 : in out Set; s2 : in Set );
  function  Intersection ( s1,s2 : Set ) return Set;

  -- DESCRIPTION :
  --   Construct the intersection of two sets of equal dimension.
  --   Either the result will be contained in s1, or a new set
  --   will be constructed.

-- SELECTORS :

  function Dimension ( s : Set ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of unknowns the set can contain.
  --   For an empty set, Dimension(s) = 0.

  function Extent_Of ( s : Set ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of unknowns the set contains.

  function Is_In ( s : Set; i : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the ith unknown belongs to the set.

  function Is_Subset ( s1,s2 : Set ) return boolean;

  -- DESCRIPTION :
  --   Returns true when the set s1 is a subset of s2.

  function Is_Equal ( s1,s2 : Set ) return boolean;

  -- DESCRIPTION :
  --   Returns true when the two sets have the same dimension
  --   and the same elements.

  generic
    with procedure Process ( sub : in Set; continue : out boolean );
  procedure Generate_Subsets ( s : in Set; k : in natural32 );

  -- DESCRIPTION :
  --   Generates all proper subsets with k elements of the set s.
  --   Each time a new subset is found, Process is called, with the
  --   subset as input parameter.  When the output parameter continue
  --   is set to false, the iteration stops.  Otherwise it continues.

  -- NOTE :
  --   The same subset will be used over and over.  If necessary,
  --   copies should be taken when processing the subset.

  generic
    with procedure Process ( sub : in Set; continue : out boolean );
  procedure Generate_All_Subsets ( s : in Set );

  -- DESCRIPTION :
  --   Generates all nonempty subsets of the set s, the set s included.
  --   As above, the procedure Process allows the processing of a subset.

-- DESTRUCTOR :

  procedure Clear ( s : in out Set );

  -- DESCRIPTION :
  --   Frees all occupied memory of the given set.
  --   After Clear(s), Dimension(s) = 0.

private

  type Set_Rep;
  type Set is access Set_Rep;

end Sets_of_Unknowns;
