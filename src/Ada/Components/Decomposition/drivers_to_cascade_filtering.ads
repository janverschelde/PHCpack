with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;

package Drivers_to_Cascade_Filtering is

-- DESCRIPTION :
--   This package contains procedures to set up the sequence of
--   homotopies to compute generic points on all components.

  procedure Driver_to_Square_and_Embed;

  -- DESCRIPTION :
  --   A system is made square by
  --      adding extra hyperplanes, if it is underdetermined; or
  --      adding slack variables, if it is overdetermined.
  --   The embedding to capture k dimensional components is executed
  --   by adding k random hyperplanes and k slack variables.
  --   This interactive routine reads in a system and creates a new
  --   square and embedded polynomial system.

  procedure Driver_to_Remove_Embedding;

  -- DESCRIPTION :
  --   Removes embed and slack variables from the system read on input.
  --   This operation undoes the squaring and embedding.

  procedure Witness_Generate
               ( outfile,resfile : in file_type;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 k : in natural32; zerotol : in double_float );

  -- DESCRIPTION :
  --   Calculates candidate witness points on every component,
  --   starting at the component of dimension k.

  -- ON ENTRY :
  --   outfile   file for intermediate results and diagnostics;
  --   resfile   file for final results:
  --               1) system at each level, for k down to 0; and
  --               2) solutions with zero slack variables;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   k         number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   zerotol   tolerance to decide whether a number is zero or not.

  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 k : in natural32; zerotol : in double_float );

  -- DESCRIPTION :
  --   This witness generate writes the witness supersets to files.

  -- ON ENTRY :
  --   name      file name for the top embedded system;
  --   outfile   file for all intermediate and final results;  
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   k         number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   zerotol   tolerance to decide whether a number is zero or not.

  procedure Driver_to_Witness_Generate;

  -- DESCRIPTION :
  --   Interactive driver to call the Witness_Generate procedure.

  procedure Black_Box_Solver
               ( file : in file_type;
                 sys : in Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean;
                 sols : out Standard_Complex_Solutions.Solution_List;
                 rc : out natural32; totaltime : out duration );
 
  -- DESCRIPTION :
  --   The black-box solver consists of the following stages:
  --     1) root-counting, when deg is true then degrees are used
  --     2) construction of start system; and
  --     3) path following to the actual target system.
  --   Several diagnostics and intermediate results are written to file.

  -- ON ENTRY :
  --   file      output file for diagnostics and intermediate results;
  --   sys       a square polynomial system;
  --   deg       flag to indicate when degrees or polytopes are used
  --             to count roots and construct the start system;
 
  -- ON RETURN :
  --   rc        root count equals the number of paths followed;
  --   sols      solutions found at the end of the paths;
  --   totaltime is the total cpu time used in computing the results.

  procedure Driver_for_Cascade_Filter
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32 );
  
  -- DESCRIPTION :
  --   Performs a cascade of homotopies to find generic points on
  --   the solution set of the system p.  k is the top dimension.
  --   All results are written to the file.

  procedure Interactive_Embed_Square_System 
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 );

  -- DESCRIPTION :
  --   Prompts the user to enter the expected top dimension, 
  --   which is returned in topdim,
  --   creates the embedded system and writes it on file.
  --   This procedure is called by Interactive_Square_and_Embed,
  --   in case the given polynomial system is square.

  procedure Embed_Square_System 
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Noninteractive version of the previous procedure.

  -- ON ENTRY :
  --   p        system with as many equations as variables;
  --   topdim   the topdimension.
 
  -- ON RETURN :
  --   embsys   system with as many extra linear equations as topdim
  --            and with the symbol table enlarged with as many
  --            symbols for the slack variables as the value of topdim.

  function Full_Embed_Nonsquare_System
              ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                nq,nv,k : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Constructs an embedding of a nonsquare system,
  --   using slices not restricted to any particular subspace.

  -- ON ENTRY :
  --   p        nonsquare polynomial system;
  --   nq       number of equations;
  --   nv       number of variables;
  --   k        number of slices to be added to the system.

  -- ON RETURN :
  --   Square polynomial system with k additional linear equations.

  procedure Interactive_Embed_Nonsquare_System
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk : in natural32;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 );

  -- DESCRIPTION :
  --   Constructs an embedding of a nonsquare system with number of
  --   equations in nbequ and number of unknowns in nbunk.
  --   The user is prompted for the expected top dimension.
  --   Slack variables are added for overdetermined systems.
  --   Dummy variables are added for underdetermined systems.
  --   The embedded system is written to file.

  procedure Embed_Nonsquare_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk,topdim : in natural32;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Noninteractive version of the above procedure.

  procedure Interactive_Square_and_Embed
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                ep : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                k : out natural32 );

  -- DESCRIPTION :
  --   The embedding of nonsquare systems involves the addition of
  --   extra slack variables (in case the system is overdetermined)
  --   or the use of dummy variables (for underdetermined systems).
  --   Therefore the embedding for square systems is treated separately
  --   from the embedding of the nonsquare systems.

  -- ON ENTRY :
  --   file     for output;
  --   p        system in any number of equations and variables.
 
  -- ON RETURN :
  --   ep       embedded system with extra slack variables;
  --   k        the dimension as entered by the user.

  procedure Square_and_Embed
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                ep : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Noninteractive version of the procedure above,
  --   works for systems in any number of equations and variables
  --   and takes care of the symbol table adjustment.

  procedure Embed_and_Cascade;

  -- DESCRIPTION :
  --   Does the embedding of the top dimension and runs the cascade.

end Drivers_to_Cascade_Filtering;
