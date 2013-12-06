with Generic_Lists;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Brackets;                          use Brackets;

package Bracket_Monomials is

-- DESCRIPTION :
--   This package provides a data abstraction to deal with
--   bracket monomials as well as operations to manipulate them.

  type Bracket_Monomial is private;
  type Array_of_Bracket_Monomials is 
    array ( integer32 range <> ) of Bracket_Monomial;

-- CONSTRUCTORS : a bracket monomial is the product of brackets.

  function Create ( b : Bracket ) return Bracket_Monomial;

  -- DESCRIPTION :
  --   Creates a bracket monomial, consisting out of one single bracket.

  procedure Multiply ( bm : in out Bracket_Monomial; b : in Bracket );

  -- DESCRIPTION :
  --   The bracket monomial is extended with another bracket.

  procedure Append ( bm : in out Bracket_Monomial; b : in Bracket );

  -- DESCRIPTION :
  --   Appends the bracket b to the bracket monomial bm.

  procedure Copy_Multiply
              ( bm1 : in Bracket_Monomial; bm2 : in out Bracket_Monomial );
  procedure Copy_Append
              ( bm1 : in Bracket_Monomial; bm2 : in out Bracket_Monomial );

  -- DESCRIPTION :
  --   Copies the first bracket monomial to the second one, either by
  --   multiplication (which also sorts) or by simply append (no sort).
  --   Note that bm2 := bm1 leads to sharing data and side effects.

  function Create ( bm : Bracket_Monomial ) return Array_of_Brackets;

  -- DESCRIPTION :
  --   Returns the array of pointers to the brackets in the bracket monomial.
  --   There is sharing of data.

-- OPERATIONS :

  function "*" ( b1,b2 : Bracket ) return Bracket_Monomial;

  -- DESCRIPTION :
  --   The resulting bracket monomial is the product of b1 and b2.

  function "*" ( bm : Bracket_Monomial; b : Bracket ) return Bracket_Monomial;
  function "*" ( b : Bracket; bm : Bracket_Monomial ) return Bracket_Monomial;

  -- DESCRIPTION :
  --   Multiplies a bracket monomial with a bracket.

  function "*" ( bm1,bm2 : Bracket_Monomial ) return Bracket_Monomial;
  procedure Multiply ( bm1 : in out Bracket_Monomial;
                       bm2 : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Multiplication of two bracket monomials.

  function Is_Null ( bm : Bracket_Monomial ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the bracket monomial is null.

  function Is_Equal ( bm1,bm2 : Bracket_Monomial ) return boolean;

  -- DESCRIPTION :
  --   Returns true when both monomials contains the same brackets.

  function Dimension ( bm : Bracket_Monomial ) return natural32;

  -- DESCRIPTION :
  --   Returns 0 if Is_Null(bm), otherwise it returns the dimension
  --   of the first bracket in the monomial.

  function "<" ( bm1,bm2 : Bracket_Monomial ) return boolean;

  -- DESCRIPTION :
  --   Let bm1 = b11*b12*..*b1k and bm2 = b21*b22*..*b2l.
  --   Then bm1 < bm2 if b1i < b2i, for some i.

  function ">" ( bm1,bm2 : Bracket_Monomial ) return boolean;

  -- DESCRIPTION :
  --   Let bm1 = b11*b12*..*b1k and bm2 = b21*b22*..*b2l.
  --   Then bm1 > bm2 if b1i > b2i, for some i.

  function Divisible ( bm : Bracket_Monomial; b : Bracket ) return boolean;

  -- DESCRIPTION :
  --   Returns true if bm is divisible by b, i.e.: if b occurs in bm.

-- ITERATORS OVER THE BRACKETS :

  function Number_of_Brackets ( bm : Bracket_Monomial ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of brackets in the monomial bm,
  --   or equivalently the degree of the bracket monomial.

  generic
    with procedure Process ( b : in Bracket; continue : out boolean );
  procedure Enumerate_Brackets ( bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Enumerate the brackets in the bracket monomial.
  --   The brackets appear lexicographically in ascending order
  --   as argument of the procedure Process.
  --   Enumeration stops when continue is set to false.

-- DESTRUCTOR :

  procedure Clear ( bm : in out Bracket_Monomial );

  -- DESCRIPTION :
  --   Deallocates the occupied memory space.

private

  package Lists_of_Brackets is new Generic_Lists(Link_to_Bracket);
  type Bracket_Monomial is new Lists_of_Brackets.List;

end Bracket_Monomials;
