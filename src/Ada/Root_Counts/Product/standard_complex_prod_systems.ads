with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Lists;        use Standard_Complex_Poly_Lists;

package Standard_Complex_Prod_Systems is

-- DESCRIPTION :
--   A product system is a systems of product polynomials.

  type Prod_Sys is array ( integer32 range <> ) of Prod_Poly;
  type Link_to_Prod_Sys is access Prod_Sys;

  function Expand ( p : Prod_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the expanded form of the product system.

-- DESTRUCTORS :

  procedure Shallow_Clear ( p : in out Prod_Sys );
  procedure Shallow_Clear ( p : in out Link_to_Prod_Sys );

  -- DESCRIPTION :
  --   Deallocation of the list structures in p,
  --   but not of the actual factors in p.

  procedure Deep_Clear ( p : in out Prod_Sys );
  procedure Deep_Clear ( p : in out Link_to_Prod_Sys );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the product system.

end Standard_Complex_Prod_Systems; 
