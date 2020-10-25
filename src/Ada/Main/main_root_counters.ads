with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Main_Root_Counters is

-- DESCRIPTION :
--   Interactive procedures to count roots for polynomial systems
--   or for Laurent systems which may have negative exponents,
--   as called in phc -r.

  procedure Polynomial_Main
               ( file : in file_type; nt : in natural32;
                 p,q : in out Poly_Sys; own : in boolean;
                 qsols : in out Solution_List; roco : out natural32;
                 verbose : in integer32 := 0 );
  procedure Laurent_Main
               ( file : in file_type; nt : in natural32;
                 p,q : in out Laur_Sys;
                 qsols : in out Solution_List; roco : out natural32;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This procedure implements an interactive driver for several
  --   root counting methods.  If the system is a Laurent system,
  --   then only the polyhedral methods are available.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   nt        number of tasks, 0 for no multitasking;
  --   p         a polynomial system;
  --   own       if true, then the user has the possibility to give
  --             an own start system;
  --   verbose   is the verbose level.

  -- ON RETURN :
  --   p         has eventually been made homogeneous;
  --   q         a start system based on the chosen root count;
  --   qsols     the solutions of q;
  --   roco      the root count.

  procedure Read_System
              ( filename : in string; lp : out Link_to_Laur_Sys );

  -- DESCRIPTION :
  --   If the string on input is not empty,
  --   then the file will be opened for reading.
  --   An exception will be handled if something goes wrong
  --   with the reading of the system on file.
  --   On return is a pointer to a polynomial system,
  --   or to null if something went wrong.

  function Is_Square ( p : Laur_Sys ) return boolean;

  -- DESRIPTOIN :
  --   Returns true if the system p has 
  --   as many variables as equations.
  --   Writes an error message and returns false
  --   if the system is not square.

  procedure Count_Roots
              ( nt : in natural32; outfilename : in string;
                p : in out Laur_Sys; v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Counts roots for a square Laurent system p,
  --   using nt tasks if nt is nonzero and writing output to the file
  --   with name in outfilename.
  --   The verbose level is in the value of v.

  procedure Count_Roots
               ( nt : in natural32; outfilename : in string;
                 p : in out Poly_Sys; v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Counts the roots of a square polynomial system p,
  --   using nt tasks if nt is nonzero and writing output to the file
  --   with name in outfilename.
  --   The verbose level is in the value of v.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Handles the input and output, checks whether the system
  --   is square and genuinely Laurent, before calling the root counters.
 
  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   verbose        the verbose level.

end Main_Root_Counters;
