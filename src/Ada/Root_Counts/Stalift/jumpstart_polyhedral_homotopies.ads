with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

package Jumpstart_Polyhedral_Homotopies is

-- DESCRIPTION :
--   This package provides drivers to jumpstart polyhedral homotopies.

  function Mixed_Volume ( file : file_type; n,r : integer32;
                          mix : Standard_Integer_Vectors.Link_to_Vector )
                        return natural32;

  -- DESCRIPTION :
  --   Returns the mixed volume for a labeled mixed cell configuration.

  -- REQUIRED :
  --   The dimensions (n,r,mix) have been read from file,
  --   so the file is positioned just before the lifted supports.

  -- ON ENTRY :
  --   file     input file for a regular mixed-cell configuration,
  --            in labeled representation;
  --   n        ambient dimension before the lifting;
  --   r        number of different supports;
  --   mix      counts the number of occurrences of each support.

  -- ON RETURN :
  --   sum of all volumes of the mixed cells on file.

  procedure Retrieve_Mixed_Volume
              ( file : in out file_type; n,r : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                mv : out natural32 );

  -- DESCRIPTION :
  --   Asks first if the user knows the mixed volume,
  --   if not the cells will be written from file and the mixed volume
  --   will be computed from the labeled mixed-cell configuration.

  -- REQUIRED :
  --   The file is open for input and the dimensions have been read.

  -- ON ENTRY :
  --   file     input file with labeled mixed-cell configuration,
  --            positioned after reading of dimensions (n,r,mix);
  --   n        ambient dimension before the lifting;
  --   r        number of different supports;
  --   mix      counts number of occurrences of each support.

  -- ON RETURN :
  --   file     may have been reset to the same position;
  --   mv       mixed volume either entered by user or computed.

  procedure Jumpstart_Polyhedral_Continuation ( p : in Poly_Sys );

  -- DESCRIPTION :
  --   Asks the user for a labeled representation of a mixed-cell
  --   configuration, tunes the parameters, and launches the trackers.

  procedure Read_Cells_and_Track 
              ( p : in Poly_Sys;
                cfile : in out file_type; ofile : in file_type );

  -- DESCRIPTION :
  --   Solves a random coefficient system to solve the system p,
  --   using a regular mixed-cell configuration on cfile.

  -- ON ENTRY :
  --   p        a target system;
  --   cfile    regular mixed-cell configuration in labeled format
  --            used to compute the mixed volume for p;
  --   ofile    output file to write start system and solutions.

  -- ON RETURN :
  --   cfile    is in/out because it may have been reset.

end Jumpstart_Polyhedral_Homotopies;
