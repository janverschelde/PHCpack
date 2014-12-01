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

-- CONSTRUCTORS :

  function Create ( p : Poly ) return Term_List;

  -- DESCRIPTION :
  --   Returns the terms that appear in p with nonzero coefficient.

  procedure Append ( first,last : in out Term_List; t : in Term );

  -- DESCRIPTION :
  --   Appends a copy of the term t to the list first, where last is 
  --   a pointer to the last element of the list.

-- COPYING :

  procedure Copy ( p : in Term_List; q : in out Term_List );

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
  procedure Shallow_Clear ( p : in out Term_List );

  -- DESCRIPTION :
  --   Deallocation of memory, a shallow clear only releases the pointer
  --   to the list; the terms themselves are not destroyed.

end Generic_Lists_of_Terms;
