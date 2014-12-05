with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Lists;
with Generic_Polynomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Polynomials is new Generic_Polynomials(Ring);

package Generic_Lists_of_Terms is 

-- DESCRIPTION :
--   This package provides polynomial systems as an array of polynomials
--   with some arithmetic operations.

  use Ring,Polynomials;

-- DATA STRUCTURE :

  package List_of_Terms is new Generic_Lists(Term);
  type Term_List is new List_of_Terms.List;
  type Array_of_Term_Lists is array ( integer32 range <> ) of Term_List;

-- CONSTRUCTORS :

  function Create ( p : Poly ) return Term_List;

  -- DESCRIPTION :
  --   Returns the terms that appear in p with nonzero coefficient.

  function Create ( t : Term_List ) return Poly;

  -- DESCRIPTION :
  --   Returns the polynomial as the sum of the terms in t.

  procedure Append ( first,last : in out Term_List; t : in Term );

  -- DESCRIPTION :
  --   Appends a copy of the term t to the list first, where last is 
  --   a pointer to the last element of the list.

  procedure Merge_Append ( first,last : in out Term_List; t : in Term );

  -- DESCRIPTION :
  --   Appends a copy of the term t to the list first, where last is
  --   a pointer to the last element of the list in case first does
  --   not contain a term with matching exponents as in t.dg.
  --   Otherwise, the coefficient of t is added to the coefficient of
  --   the term with matching exponents in first.

  procedure Concat ( first,last : in out Term_List; t : in Term_List );

  -- DESCRIPTION :
  --   Concatenates the terms in t to the list first, where last points
  --   to the last term in the list first.

  procedure Merge_Concat ( first,last : in out Term_List; t : in Term_List );

  -- DESCRIPTION :
  --   Concatenates the terms in t to the list first, where last points
  --   to the last term in the list first, each time calling Merge_Append
  --   to take into account matching exponents.

-- COPYING :

  procedure Copy ( p : in Term_List; q,q_last : in out Term_List );

  -- DESCRIPTION :
  --   Copies the list p to q, where q_last on return points to the
  --   last element of the list q.

-- SELECTORS :

  function Is_In ( p : in Term_List; t : in Term ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the term t belongs to the list p.

  generic
    with procedure process ( t : in Term; continue : out boolean );
  procedure Iterator ( p : in Term_List );

  -- DESCRIPTION :
  --   Applies the procedure process to each term in the list t
  --   as long as the continue remains true.

-- DESTRUCTORS :

  procedure Clear ( p : in out Term_List );
  procedure Clear ( p : in out Array_of_Term_Lists );
  procedure Shallow_Clear ( p : in out Term_List );
  procedure Shallow_Clear ( p : in out Array_of_Term_Lists );

  -- DESCRIPTION :
  --   Deallocation of memory, a shallow clear only releases the pointer
  --   to the list; the terms themselves are not destroyed.

end Generic_Lists_of_Terms;
