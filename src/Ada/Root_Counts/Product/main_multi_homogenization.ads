with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Main_Multi_Homogenization is

-- DESCRIPTION :
--   The main interactive procedure to compute multi-homogeneous 
--   Bezout numbers are defined by this package.

  procedure Multi_Homogenization_Info;

  -- DESCRIPTION :
  --   Displays information about multi-homogenization on screen.

  procedure Menu_Prompt ( choice : out character; gb : in natural32 );

  -- DESCRIPTION :
  --   Prompts for a choice after showing a menu.
  --   The generalized Bezout number is given in gb.

  procedure Menu_Handle
              ( file : in file_type; choice : in character;
                p : in Poly_Sys; gb : in out natural32 );

  -- DESCRIPTION :
  --   Handles the menu choice of the Menu_Prompt to compute
  --   a generalized Bezout number gb for the system p.

  procedure Define_Partitions
              ( file : in file_type; p : in Poly_Sys;
                gb : in out natural32; b : out natural32 );

  -- DESCRIPTION :
  --   Defines the partitions for the set of variables of the polynomials
  --   in the system b to compute a generalized Bezout number.
  --   Output is written to file.  The gb is the working variable
  --   for the generalized Bezout number, its ends value is in b.

  procedure Construct_Start_System
              ( file : in file_type; p : in Poly_Sys;
                q : out Poly_Sys; rq : out Prod_Sys;
                qsols : out Solution_List );

  -- DESCRIPTION :
  --   Constructs a multi-homogeneous start system q,
  --   represented in product form by rq, with solutions in qsols,
  --   for the system p.  Output is written to file.

  procedure Write_Results
              ( file : in file_type; p : in Poly_Sys; gb : in natural32 );

  -- DESCRIPTION :
  --   Writes the multi-homogeneous Bezout number gb
  --   for the system p to file.

  procedure Save_Results ( q : in Prod_Sys; sols : in Solution_List );

  -- DESCRIPTION :
  --   For a nonzero list of solutions, the user is prompted to give
  --   a file name to where to write the start system and its solutions.

  procedure Main ( file : in file_type; p : in Poly_Sys;
                   b : in out natural32;
                   q : out Poly_Sys; rq : out Prod_Sys;
                   qsols : out Solution_List );

  -- DESCRIPTION :
  --   The main interactive procedure for multi-homogenization.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   p         a polynomial system.

  -- ON RETURN :
  --   b         a bound based on the degree structure;
  --   q         a random product start system;
  --   rq        product format of the start system q;
  --   qsols     the solutions of q.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a polynomial system, an output file,
  --   and then computes.

end Main_Multi_Homogenization;
