with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Main_SAGBI_Homotopies is

-- DESCRIPTION :
--   This package offers two interactive procedures to launch
--   the SAGBI homotopies to intersect m-planes in n-space, m = n-d.

  procedure Main ( n,d : in natural32 );
  procedure Main ( file : in file_type; n,d : in natural32 );

  -- DESCRIPTION :
  --   Computes the intersection of (n-d)-planes in n-space.
  --   To obtain a complete intersection, (n-d)*d planes are needed.

  -- ON ENTRY :
  --   file      output file (is optional) must be opened for output;
  --   n         dimension of the ambient space;
  --   d         co-dimension of the planes we intersect.

end Main_SAGBI_Homotopies;
