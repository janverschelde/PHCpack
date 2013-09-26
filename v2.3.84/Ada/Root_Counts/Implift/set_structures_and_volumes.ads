with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Trees_of_Vectors;                   use Trees_of_Vectors;

package Set_Structures_and_Volumes is

-- DESCRIPTION :
--   This package contains routines for implementing the
--   hybrid approach of working with set structures and
--   computing mixed volumes.

  procedure Build_RPS ( k,n : in natural32; p : in Poly_Sys );

  -- DESCRIPTION :
  --   If the set structure is empty, then a set structure will
  --   be constructed for the first k polynomials of p.
  --   This allows the construction of random product polynomials.

  -- ON ENTRY :
  --   k        number of polynomials for which a random product structure
  --            needs to be constructed, 0 < k <= n;
  --   n        total number of polynomials in the system;
  --   p        polynomial system.

  -- ON RETURN :
  --   The result of this operation can be found in
  --   the data of the package Random_Product_System. 

  function Bezout_BKK ( k,n : natural32; p : Poly_Sys ) return natural32;

  function Bezout_BKK ( k,n : natural32; p : Poly_Sys; tv : Tree_of_Vectors )
		      return natural32;

  procedure Bezout_BKK ( k,n : in natural32; p : in Poly_Sys;
			 tv : in out Tree_of_Vectors; bb : out natural32 );

  -- DESCRIPTION :
  --   If the set structure is empty, then a set structure will be
  --   constructed for the first k equations.
  --   This set structure will be used to eliminate k unknowns,
  --   for the last equations of p, the BKK bound will be computed.

  -- ON ENTRY :
  --   k        number of polynomials for which a random product
  --            structure needs to be constructed, 0 <= k <= n;

  -- ON RETURN :
  --   tv       tree of useful directions used to compute mixed volume,
  --            once tv has been constructed, the mixed volume
  --            can be computed more efficiently;
  --   bb       if k = 0, then the BKK bound is returned,
  --            if 0 < k < n, then for the first k polynomials a product
  --            structure has been used to eliminate k unknowns,
  --            if k = n, then the Bezout number based on the set structure
  --            will be returned;
 
  procedure Mixed_Solve 
               ( file : in file_type; k,n : in natural32; p : in Poly_Sys;
                 bb : out natural32;
                 g : in out Poly_Sys; sols : in out Solution_List );

  procedure Mixed_Solve
               ( file : in file_type; k,n : in natural32; p : in Poly_Sys;
                 tv : in Tree_of_Vectors;
                 bb : out natural32; g : in out Poly_Sys;
                 sols : in out Solution_List );

  -- DESCRIPTION :
  --   Constructs a random product coefficient start system for p.

  -- ON ENTRY :
  --   file      needed to write intermediate results on;
  --   k         a number between 0 and n;
  --   n         the number of variables and equations;
  --   p         a polynomial system;
  --   tv        the tree of useful directions.

  -- ON RETURN :
  --   g         a random product coefficient start system;
  --             if k = 0, then a random coefficient start system,
  --             if k = n, then a random product start system;
  --   bb        the Bezout BKK bound for the system p;
  --             if k = 0, then bb = BKK bound,
  --             if k = n, then bb = Bezout number based on set structure;
  --   sols      the solutions of g.

end Set_Structures_and_Volumes;
