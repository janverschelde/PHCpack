with text_io;                            use text_io;
with Multprec_Irreducible_Decomp;        use Multprec_Irreducible_Decomp;

package Multprec_Irreducible_Decomp_io is

-- DESCRIPTION :
--   This package provides output routines of hypersurfaces
--   and solution components.

  procedure put ( h : in Hypersurfaces );
  procedure put ( file : in file_type; h : in Hypersurfaces );

  procedure put ( s : in Solution_Components );
  procedure put ( file : in file_type; s : in Solution_Components );

end Multprec_Irreducible_Decomp_io;
