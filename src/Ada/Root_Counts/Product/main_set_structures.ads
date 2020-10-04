with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
-- with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Main_Set_Structures is

-- DESCRIPTION :
--   The main interactive procedure to compute Bezout numbers for
--   general linear-product start systems are defined by this package.

  procedure Set_Structure_Info;

  -- DESCRIPTION :
  --   Displays information on set structures on screen.

  procedure Read_Set_Structure ( n : in natural32 );

  -- DESCRIPTION :
  --   Prompts for a set structure for a system of n equations.

  procedure Menu_Prompt ( choice : out character; bb : in natural32 );

  -- DESCRIPTION :
  --   Prompts for a choice after displaying a menu.

  procedure Menu_Handler ( file : in file_type; choice : in character;
                           p : in Poly_Sys; bb : in out natural32 );

  -- DESCRIPTION :
  --   Depending on whether the user wants the computer to generate
  --   a set structure or whether a set structure needs to be read,
  --   a generalized Bezout number is computed and returned in bb.
  --   The permanent computation based on the bipartite matching problem
  --   is more efficient than the one based on set unions.

  procedure Compute_Bezout_Number
              ( file : in file_type;
                p : in Poly_Sys; bb : in out natural32 );

  -- DESCRIPTION :
  --   Interactive calculation of a generalized Bezout number bb,
  --   based on a supporting set structure for the polynomial system p.
  --   The result is written to file, which much be opened for output.

  procedure Construct_Start_System
              ( file : in file_type; p : in Poly_Sys;
                bb : in natural32;
                q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Constructs a linear-product start system q for p,
  --   with its solutions in qsols, based on a supporting
  --   set structure for p with Bezout number in bb.
  --   Output is written to file.

  procedure Write_Results ( file : in file_type; bb : in natural32 );

  -- DESCRIPTION :
  --   Writes the generalized Bezout number with its corresponding
  --   set structure to file.

  procedure Save_Results ( q : in Prod_Sys; qsols : in Solution_List );

  -- DESCRIPTION :
  --   Prompts for the name of an output file and writes the start
  --   system q and its solution to that file.

  procedure Main ( file : in file_type; p : in Poly_Sys;
                   b : in out natural32; -- lpos : in out List;
                   q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Allows the interactive computation of a generalized Bezout number,
  --   with an optional construction of a start system.

  -- ON ENTRY :
  --   file      output file;
  --   p         a polynomial system.

  -- ON RETURN :
  --   b         a bound based on the set structure;
  --   lpos      a list of positions indicating the acceptable classes;
  --   q         a random product start system;
  --   qsols     the solutions of q.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a polynomial system, an output file,
  --   and then computes.

end Main_Set_Structures;
