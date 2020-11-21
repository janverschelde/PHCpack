with text_io;                            use text_io;
with Ada.Calendar;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with String_Splitters;                   use String_Splitters;
with Standard_Integer_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials; 
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Monomial_Maps;             use Standard_Monomial_Maps;

package Standard_Blackbox_Solvers is

-- DESCRIPTION :
--   This package provides procedure to solve polynomial systems
--   in blackbox mode, in double precision.

  procedure Write_Toric_Binomial_Solutions
              ( file : in file_type; d : in natural32;
                M : in Standard_Integer_Matrices.Matrix;
                c : in Solution_List );

  -- DESCRIPTION :
  --   Writes the solutions of the binomial system to file.

  procedure Append_Toric_Binomial_Solutions_to_Input_File
              ( name : in string; d : in natural32;
                M : in Standard_Integer_Matrices.Matrix;
                c : in Solution_List );

  -- DESCRIPTION :
  --   Appends solutions to the file with the given name.

  procedure Append_Affine_Binomial_Solutions_to_Input_File
              ( name : in string;
                c : in Link_to_Array_of_Monomial_Map_Lists );

  -- DESCRIPTION :
  --   Appends solution maps to the file with the given name.

  procedure Toric_Binomial_Solver
              ( nt : in natural32; start_moment : in Ada.Calendar.Time;
                p : in Laur_Sys; append_sols : in boolean;
                infilename,outfilename : in string; outfile : out file_type;
                to_file,fail : out boolean; v : in integer32 := 0 );

  -- DESCRIPTION :
  --   If p is a binomial system, then it is solved.

  procedure Affine_Binomial_Solver
              ( nt : in natural32; start_moment : in Ada.Calendar.Time;
                p : in Laur_Sys; append_sols : in boolean;
                infilename,outfilename : in string; outfile : out file_type;
                outnewname : out Link_to_String;
                to_file,fail : out boolean; v : in integer32 := 0 );

  -- DESCRIPTION :
  --   If p is a binomial system, then it is solved.

  procedure Toric_Binomial_Solver
              ( nt : in natural32; start_moment : in Ada.Calendar.Time;
                p : in Poly_Sys; append_sols : in boolean;
                infilename,outfilename : in string; outfile : out file_type;
                outnewname : out Link_to_String;
                to_file,fail : out boolean; v : in integer32 := 0 );

  -- DESCRIPTION :
  --   If p is a binomial system, then it is solved.

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Link_to_Poly_Sys; append_sols : in boolean;
                    v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the blackbox solver for a polynomial system.

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Link_to_Laur_Sys; append_sols : in boolean;
                    v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the blackbox solver for a Laurent polynomial system.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Note about the difference between Laurent and ordinary polynomials:
  --   For Laurent binomial systems (the genuine ones with negative powers),
  --   a stable mixed volume or an affine solution set does not make sense.
  --   For ordinary binomial systems of positive dimension, we compute all
  --   affine maps, or a stable mixed volume is computed.

end Standard_Blackbox_Solvers;
