with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;

package Monodromy_Group_Actions is

-- DESCRIPTION :
--   This package records the actions of the monodromy group to
--   decompose the k-dimensional solution set of a polynomial system
--   into irreducible components.

  type Irreducible_Components is private;

-- CONSTRUCTORS :

  function Create ( n : integer32 ) return Irreducible_Components;

  -- DESCRIPTION :
  --   The decomposition on return consists of n linear components.
  --   Until the structure on return is cleared, it sum of degrees
  --   will always equal n.

  procedure Add ( ic : in out Irreducible_Components;
                  i : in integer32; j : in natural32 );

  -- DESCRIPTION :
  --   Adds j to the i-th component.

  procedure Merge ( ic : in out Irreducible_Components;
                    i,j : in integer32 );

  -- DESCRIPTION :
  --   Merges the j-th set with the set corresponding to the i-th variable. 

  procedure Act ( ic : in out Irreducible_Components; map : in Vector );

  -- DESCRIPTION :
  --   Records the action of the map on the decomposition where
  --   map(i) = j indicates i and j belong to the same component.

-- SELECTORS :

  function Empty ( ic : Irreducible_Components; i : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the i-th set in the irreducible components is empty;
  --   otherwise, false is returned

  function Cardinality ( ic : Irreducible_Components; i : integer32 )
                       return natural32;

  -- DESCRIPTION :
  --   Returns the cardinality of the i-th set in the irreducible components.
  --   If the set is empty, then 0 is returned.

  function Component ( ic : Irreducible_Components; i : integer32 )
                     return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..Cardinality(ic,i) with labels to those
  --   points that belong to the same component.

  function Is_In ( ic : Irreducible_Components;
                   i : integer32; j : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if j belongs to the i-th component.

  function Sum_of_Degrees ( ic : Irreducible_Components ) return integer32;

  -- DESCRIPTION :
  --   Returns the sum of the degrees of the components.
  --   This is also an upper bound on the number of components.

  function Empty ( ic : Irreducible_Components ) return boolean;

  -- DESCRIPTION :
  ---  Returns true if there are no irreducible components. 

  function Cardinality ( ic : Irreducible_Components ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of irreducible components.

  function Degrees ( ic : Irreducible_Components ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector with the degrees of the components,
  --   of range 1..Cardinality(ic), sorted in increasing order.

  function Nonempty_Sets ( ic : Irreducible_Components ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..Cardinality(ic) with in its entries
  --   the label of the representative generic point.

  function Representatives ( ic : Irreducible_Components; k : natural32 )
                           return Vector;

  -- DESCRIPTION :
  --   Returns the k-th representative in every nonempty set, where
  --   k is taken mod the cardinality of the set, plus 1 if the result
  --   of the mod operation is zero.

  function Representatives ( ic : Irreducible_Components; k,j : integer32 )
                           return Vector;

  -- DESCRIPTION :
  --   Returns representatives in every nonempty set, starting at the
  --   k-th element (eventually modulo the cardinality) and all other
  --   elements within jump space from the starting element.
  --   In particular, if j = 1, then all elements will be taken;
  --   if j = 2, then only every other element is representative;
  --   for j = Sum_of_Degree(ic) or j = 0, then every set is only 
  --   represented once.

-- DESTRUCTOR :

  procedure Clear ( ic : in out Irreducible_Components );

  -- DESCRIPTION :
  --   All occupied memory is released.

private

  type Irreducible_Components_Rep;
  type Irreducible_Components is access Irreducible_Components_Rep;

end Monodromy_Group_Actions;
