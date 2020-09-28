with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Black_Box_Mixed_Volumes is

-- DESCRIPTION :
--   If the focus in the black box solver is on polyhedral homotopies,
--   then the black box root counters are specialized to polyhedral methods.

  procedure Mixed_Volume
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 mivo,stmv : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 );
  procedure Mixed_Volume
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 mivo,stmv : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the mixed volume and the stable mixed volume for p.

  -- ON ENTRY :
  --   file      to write mixed volumes and timing information;
  --   p         a polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   mivo      mixed volume;
  --   stmv      stable mixed volume;
  --   lifsup    lifted supports of the system;
  --   stlb      lifting of the artificial origin;
  --   mix       type of mixture;
  --   perm      permutation of the equations in p;
  --   iprm      induced permutation of the blackbox mixed volume calculator;
  --   orgmcc    regular mixed-cell configuration to compute mivo;
  --   stbmcc    extra stable mixed cells that contribute to stmv;
  --   rocotime  is the time it took to compute the root count.

  procedure Construct_Start_System
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : in Link_to_Vector;
                 stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                 orgmcc,stbmcc : in Mixed_Subdivision;
                 q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                 hocotime : out duration; verbose : in integer32 := 0 );
  procedure Construct_Start_System
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : in Link_to_Vector;
                 stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                 orgmcc,stbmcc : in Mixed_Subdivision;
                 q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                 hocotime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs polyhedral homotopies to solve a random coefficient system.

  -- ON ENTRY :
  --   file      to write start system and timing information;
  --   p         polynomial system;
  --   mix       type of mixture of the supports;
  --   stlb      lifting for the artificial origin;
  --   lifted    lifted supports;
  --   orgmcc    regular mixed-cell configuration to compute mv;
  --   stbmcc    extra stable mixed cells that contributed to smv;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   q         start system;
  --   qsols     solutions of q;
  --   qsols0    solutions of q with zero components;
  --   hocotime  is the time it took to construct the start system.

  procedure Black_Box_Polyhedral_Homotopies
               ( silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the mixed volume and solves a random coefficient start system
  --   with polyhedral homotopies.

  -- ON ENTRY :
  --   silent    if silent, then the computed root counts will not be shown
  --             on screen, otherwise, the root counter remains silent;
  --   p         a polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, mixed volume or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Polyhedral_Homotopies
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the mixed volume and solves a random coefficient start system
  --   with polyhedral homotopies.
  --   Returns a string with mixed volume and stable mixed volume.

  -- ON ENTRY :
  --   p         a polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, mixed volume or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   rocos     string with the information about the root counts,
  --             with the same information as the above procedures
  --             when silent equals false;
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Polyhedral_Homotopies
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the mixed volume and solves a random coefficient start system
  --   with polyhedral homotopies.
  --   Writes mixed volume, stable mixed volume, start system, 
  --   start solutions, and timing information to file.

  -- ON ENTRY :
  --   p         a polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, mixed volume or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   rocos     string with the information about the root counts,
  --             with the same information as the above procedures
  --             when silent equals false;
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Polyhedral_Homotopies
               ( nt : in natural32; silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the mixed volume and solves a random coefficient start system
  --   with polyhedral homotopies, with multitasking.

  -- ON ENTRY :
  --   nt        number of tasks, must be at least 2;
  --   silent    if not silent, then writes the mixed volume
  --             and the stable mixed volume as root counts;
  --   p         a polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, mixed volume or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Polyhedral_Homotopies
               ( nt : in natural32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the mixed volume and solves a random coefficient start system
  --   with polyhedral homotopies, with multitasking.

  -- ON ENTRY :
  --   nt        number of tasks, must be at least 2;
  --   p         a polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, mixed volume or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   rocos     string with the information about the root counts,
  --             with the same information as the above procedures
  --             when silent equals false;
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Polyhedral_Homotopies
               ( file : in file_type; nt : in natural32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the mixed volume and solves a random coefficient start system
  --   with polyhedral homotopies, with multitasking.
  --   Writes mixed volume, stable mixed volume, start system, 
  --   start solutions, and timing information to file.

  -- ON ENTRY :
  --   nt        number of tasks, must be at least 2;
  --   p         a polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, mixed volume or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

end Black_Box_Mixed_Volumes;
