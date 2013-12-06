with Generic_Lists;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Standard_Complex_Laur_Lists is

-- DESCRIPTION :
--   This package provides lists of Laurent polynomials
--   with standard complex coefficients.

  package List_of_Laurentials is new Generic_Lists(Poly);
  type Poly_List is new List_of_Laurentials.List;

-- CONSTRUCTORS :

  function Create ( p : Poly ) return Poly_List;

  -- DESCRIPTION :
  --   Returns a list with as element p.

  function Create_Copy ( p : Poly ) return Poly_List;

  -- DESCRIPTION :
  --   Returns a list with as element a copy of the polynomial p.

  function Create ( p : Laur_Sys ) return Poly_List;

  -- DESCRIPTION :
  --   Returns a shallow copy of the polynomials in p.

  function Create ( p : Poly_List ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns a shallow copy of the polynomials in p.

  function Create_Copy ( p : Laur_Sys ) return Poly_List;

  -- DESCRIPTION :
  --   Every polynomial in the system p is copied into
  --   the list on return.

  function Create_Copy ( p : Poly_List ) return Laur_Sys;
 
  -- DESCRIPTION :
  --   Makes a copy of every polynomial in the list p
  --   and places the copies in the equations of the Laurent system.

  procedure Append ( first,last : in out Poly_List; p : in Poly );

  -- DESCRIPTION :
  --   Appends the polynomial p to the list headed by first
  --   and with last element in last.

  procedure Append_Copy ( first,last : in out Poly_List; p : in Poly );

  -- DESCRIPTION :
  --   Makes a copy of the polynomial p and appends this copy
  --   to the list headed by first and with last element in last.

  procedure Concatenate ( first,last : in out Poly_List; pl : in Poly_List );

  -- DESCRIPTION :
  --   Concatenates the polynomials in the list pl to the list first.

-- DESTRUCTORS :

  procedure Shallow_Clear ( pl : in out Poly_List );

  -- DESCRIPTION :
  --   This shallow clear deallocates the list structure,
  --   but not the polynomials inside the list.

  procedure Deep_Clear ( pl : in out Poly_List );

  -- DESCRIPTION :
  --   Does a deep clear of the list of Laurent polynomials in pl:
  --   space occupied by all polynomials in the list is released.

end Standard_Complex_Laur_Lists;
