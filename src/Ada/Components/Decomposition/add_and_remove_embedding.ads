with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Add_and_Remove_Embedding is

-- DESCRIPTION :
--   An embedding of a polynomial system is applied to help computing
--   witness set representations of polynomial systems.
--   The procedures in this package remove the embedding.

  procedure Standard_Square_and_Embed ( iptname,optname : in string );
  procedure DoblDobl_Square_and_Embed ( iptname,optname : in string );
  procedure QuadDobl_Square_and_Embed ( iptname,optname : in string );

  -- DESCRIPTION :
  --   A system is made square by
  --      adding extra hyperplanes, if it is underdetermined; or
  --      adding slack variables, if it is overdetermined.
  --   The embedding to capture k dimensional components is executed
  --   by adding k random hyperplanes and k slack variables.
  --   This interactive routine reads in a system and creates a new
  --   square and embedded polynomial system,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   iptname   name of the input file, passed as command line argument;
  --   optname   name of the output file, given as command line argument.

  procedure Driver_to_Square_and_Embed ( iptname,optname : in string );

  -- DESCRIPTION :
  --   Prompts the user for the level of precision and then calls
  --   the proper driver procedure.

  -- ON ENTRY :
  --   iptname   name of the input file, passed as command line argument;
  --   optname   name of the output file, given as command line argument.

  procedure Standard_Remove_Embedding
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k,ns : in natural32 );
  procedure Standard_Remove_Embedding
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k,ns : in natural32 );
  procedure DoblDobl_Remove_Embedding
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 );
  procedure DoblDobl_Remove_Embedding
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 );
  procedure QuadDobl_Remove_Embedding
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 );
  procedure QuadDobl_Remove_Embedding
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 );

  -- DESCRIPTION :
  --   Removes the embedding from the system p and solutions sols,
  --   writing the resulting system and solutions to file.

  -- ON ENTRY :
  --   file     file opened for output;
  --   p        embedded system, with slack variables;
  --   sols     solutions of the embedded system p;
  --   k        number of embed variables, equals the dimension;
  --   ns       number of slack variables to make the system square.

  procedure Standard_Remove_Embedding ( inpname,outname : in string );
  procedure DoblDobl_Remove_Embedding ( inpname,outname : in string );
  procedure QuadDobl_Remove_Embedding ( inpname,outname : in string );

  -- DESCRIPTION :
  --   Removes embed and slack variables from the system read on input.
  --   This operation undoes the squaring and embedding,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   inpname   name of the input file, passed as command line argument;
  --   outname   name of the output file, given as command line argument.

  procedure Driver_to_Remove_Embedding ( inpname,outname : in string );

  -- DESCRIPTION :
  --   Prompts the user for the level of precision and then calls
  --   the proper driver procedure.

  -- ON ENTRY :
  --   inpname   name of the input file, passed as command line argument;
  --   outname   name of the output file, given as command line argument.

end Add_and_Remove_Embedding;
