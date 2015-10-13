with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;

package Brackets is

-- DESCRIPTION :
--   A bracket is a selection of rows and represents a maximal minor
--   which is use in the Pluecker embedding of the Grassmannian.

-- DATA STRUCTURES :

  type Bracket is new Standard_Natural_Vectors.Vector;
                                                 -- typically of range 1..d
  type Link_to_Bracket is access Bracket;

  type Array_of_Brackets is array ( integer32 range <> ) of Link_to_Bracket;

-- CONSTRUCTORS : puts the bracket in "normal form" with all entries sorted.

  procedure Create ( v : in Standard_Natural_Vectors.Vector;
                     b : out Bracket; sign : out integer32 );

  -- DESCRIPTION :
  --   Sorts the entries in the rows in ascending order.
  --   Eventually the sign of the bracket coefficient changes,
  --   when the number of transpositions is odd.

  procedure Create ( v : in Standard_Natural_Vectors.Vector;
                     perm : out Standard_Natural_Vectors.Vector;
                     b : out Bracket; sign : out integer32 );

  -- DESCRIPTION :
  --   Returns also the permutation "perm" that is used in the sorting.

  function Modulo ( b : Bracket; n : natural32 ) return Bracket;

  -- DESCRIPTION :
  --   Returns the modular form of the bracket: every entry mod n.

  procedure Modulo ( b : in Bracket; n : in natural32;
                     perm : out Standard_Natural_Vectors.Vector;
                     mb : out Bracket );

  -- DESCRIPTION :
  --   mb = Modulo(b,n) and perm is the permutation used to sort mb.

-- SELECTORS : all require that the brackets are sorted.

  function Is_Zero ( b : Bracket ) return boolean;

  -- DESCRIPTION :
  --   The bracket is zero when it has repeated entries.

  function Is_Equal ( b1,b2 : Bracket ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both brackets are the same.

  function "<" ( b1,b2 : Bracket ) return boolean;

  -- DESCRIPTION :
  --   Returns true if b1 precedes b2 lexicographically.

  function ">" ( b1,b2 : Bracket ) return boolean;

  -- DESCRIPTION :
  --   Returns true if b1 is larger than b2, lexicographically.

  function Is_Standard ( b1,b2 : Bracket ) return natural32;

  -- DESCRIPTION :
  --   Let s = Is_Standard(b1,b2).  If s = 0, then the pair b1*b2 constitutes
  --   a Standard bracket monomial, otherwise, s indicates the first entry
  --   where b1(s) > b2(s).

  -- REQUIRED : b1 < b2, b1 precedes b2 lexicographically.

  function "<=" ( alpha,beta : Bracket ) return boolean;

  -- DESCRIPTION :
  --   We have that alpha <= beta if for alpha(i) <= beta(i)
  --   for all i in alpha'range = beta'range.

-- DESTRUCTORS :

  procedure Clear ( b : in out Link_to_Bracket );
  procedure Clear ( b : in out Array_of_Brackets );

  -- DESCRIPTION :
  --   Destroys the allocate space.

end Brackets;
