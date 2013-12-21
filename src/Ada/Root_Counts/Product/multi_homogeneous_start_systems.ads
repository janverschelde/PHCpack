with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Multi_Homogeneous_Start_Systems is

-- DESCRIPTION :
--   This package contains a routine for the construction
--   of random product start systems.

  procedure RPQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List; nl : out natural32 );

  -- DESCRIPTION :
  --   This routine constructs a start system q with the same structure
  --   as the system p.
  --   Solving q happens by factoring all possible linear systems.

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   q          a suitable start system for p;
  --   sols       the solutions of the start system;
  --   nl         the number of matrices that are factored.

  procedure GBQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List );
  procedure GBQ ( p : in Poly_Sys; q : out Poly_Sys; rq : out Prod_Sys;
                  sols : in out Solution_List );

  -- DESCRIPTION :
  --   This routine constructs a start system q with the same structure
  --   as the system p.
  --   Solving q happens by factoring only these linear systems that
  --   correspond to admissible products.

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   q          a suitable start system for p;
  --   rq         start system in product format;
  --   sols       the solutions of the start system.

end Multi_Homogeneous_Start_Systems;
