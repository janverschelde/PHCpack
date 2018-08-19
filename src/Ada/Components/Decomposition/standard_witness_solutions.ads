with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;

package Standard_Witness_Solutions is

-- DESCRIPTION :
--   This package stores the witness set representation for the
--   numerical irreducible decomposition of the solution set of
--   an ordinary or a Laurent polynomial system.
--   The data stored by this package is computed by the procedures
--   in the embeddings_and_cascades package, which defines the
--   blackbox solver for polynomial systems, as called by phc -B,
--   computed in standard double precision.

-- CONSTRUCTOR :

  procedure Initialize ( nbvar,topdim : in natural32 );

  -- DESCRIPTION :
  --   Initialized the package with the number of variables in nbvar
  --   and the top dimension in topdim.

  procedure Save_Embedded_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in natural32 );
  procedure Save_Embedded_System
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                k : in natural32 );

  -- DESCRIPTION :
  --   Stores the embedded system p at position k in the sequence
  --   of the systems in the cascade used to compute all witness points.

-- SELECTORS :

  function Load_Embedded_System
              ( k : natural32 )
              return Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  function Load_Embedded_System
              ( k : natural32 )
              return Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  -- DESCRIPTION :
  --   Returns the embedded system at position k in the sequence
  --   of the systems in the cascade used to compute all witness points.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Deallocates the data stored in the package.

end Standard_Witness_Solutions;
