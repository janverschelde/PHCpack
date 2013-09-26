with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with generic_lists;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;

package Standard_Complex_Poly_Lists is

-- DESCRIPTION :
--   A product of polynomials is stored as a list of polynomials.

  package Lists_of_Standard_Complex_Polynomials is new generic_lists(Poly);
  type Prod_Poly is new Lists_of_Standard_Complex_Polynomials.List;

  function Create ( p : Poly ) return Prod_Poly;

  -- DESCRIPTION :
  --   Returns a list with one element: the given polynomial p.

  function Number_of_Factors ( p : Prod_Poly ) return natural32;

  -- DESCRIPTION :
  --   Returns the length of the list of factors.

  function Number_of_Unknowns ( p : Prod_Poly ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of unknowns in the first factor,
  --   zero if there are no factors.
  --   Note that the number of unknowns is expected to be the same
  --   in every factor of the product.

  function Expand ( p : Prod_Poly ) return Poly;

  -- DESCRIPTION :
  --   Returns the expanded form of the product polynomial p.

  procedure Shallow_Clear ( p : in out Prod_Poly );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the list structure,
  --   but the deallocation is shallow: the polynomials inside the
  --   list are not cleared.

  procedure Deep_Clear ( p : in out Prod_Poly );

  -- DESCRIPTION :
  --   Deallocates all memory occupied by the list, including the
  --   memory occupied by the polynomials in the list.

end Standard_Complex_Poly_Lists;
