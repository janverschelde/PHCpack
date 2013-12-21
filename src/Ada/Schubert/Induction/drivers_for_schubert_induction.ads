with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Bracket_Monomials;                  use Bracket_Monomials;
with Intersection_Posets;                use Intersection_Posets;

package Drivers_for_Schubert_Induction is

-- DESCRIPTION :
--   Offers drivers to the Littlewood-Richardson homotopies.

  function Prompt_for_Bracket_Monomial return Bracket_Monomial;

  -- DESCRIPTION :
  --   Prompts the user for a product of brackets,
  --   i.e.: for a Schubert intersection problem.

  function Is_Isolated 
              ( b : Standard_Natural_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   A bracket represents an isolated solution if b(i) = i
  --   for all i in b'range.

  function Number_of_Isolated_Solutions
              ( ips : Intersection_Poset ) return Natural_Number;

  -- DECRIPTION :
  --   Returns the number of isolated solutions at the top level
  --   of an intersection poset of a resolved intersection condition.

  procedure Resolve_Intersection_Condition ( n : in natural32 );

  -- DESCRIPTION :
  --   Interactive procedure to resolve a Schubert condition in n-space.
  --   Prompts the user for a condition and shows its resolution on screen.

  function Get_Intersection_Conditions 
              ( k : natural32 ) return Standard_Natural_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Prompts the user for the number of intersection conditions.

  procedure Read_Intersection_Conditions
              ( ip : in Standard_Natural_Vectors.Vector;
                rows,cols : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Interactive routine to read two brackets representing a
  --   Schubert intersection condition.  As long as they do not form
  --   a happy configuration, the users is prompted to re-enter.

  -- ON ENTRY :
  --   ip       the identity permutation.

  -- ON RETURN :
  --   rows     position of rows of white checkers;
  --   cols     position of columns of white checkers.

  procedure Write_Results
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec;
                vfs : in Standard_Complex_VecMats.VecMat;
                sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the polynomial system and the solutions to file.

  -- ON ENTRY :
  --   file     output file, must be opened for output;
  --   n        ambient dimension;
  --   k        dimension of the solution planes;
  --   q        permutation defines the location of the black checkers;
  --   rows     row positions for white checkers
  --   cols     columns of white checkers of resolved condition;
  --   cnds     conditions kept fixed during flag continuation;
  --   vfs      fixed flags, vfs'range = cnds'range;
  --   sols     solution k-planes.

  procedure Run_Moving_Flag_Continuation ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts the users first for rows, columns, and the third fixed
  --   intersection condition on k-planes and n-space before calling
  --   the other Run_Moving_Flag_Continuation below.

  procedure Run_Moving_Flag_Continuation
              ( n,k : in integer32;
                rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Runs the Littlewood-Richardson homotopies to compute k-planes
  --   meeting generic flags in n-space along a specific attitude.
  --   All output is written to a solution file.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution planes;
  --   rows     row positions for white checkers
  --   cols     columns of white checkers of resolved condition;
  --   cnds     conditions kept fixed during flag continuation.

  procedure Set_Symbol_Table
              ( n,k : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Initializes the symbol table for a k-plane in n-space, for
  --   black and white checkers respectively defined by p and (rows,cols).

  -- IMPORTANT :
  --   This initialization must be done for running cheater's homotopy,
  --   using start solutions of generic problems read from file.

  procedure Scan_for_Start_Schubert_Problem
              ( file : in file_type; n : in integer32;
                vf : out Standard_Complex_VecMats.VecMat;
                p : out Link_to_Poly_Sys; sols : out Solution_List;
                fail : out boolean );

  -- DESCRIPTION :
  --   A Schubert problem consists of a set of fixed flags in n-space,
  --   a polynomial system and a list of solutions.
  --   If fail is false on return, then (vf,p,sols) define a start
  --   Schubert problem stored as the output of Littlewood-Richardson
  --   homotopies on the given file.

  procedure Run_Cheater_Flag_Homotopy
              ( n,k : in integer32;
                rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec;
                inpt : in boolean );

  -- DESCRIPTION :
  --   Runs a cheater's homotopy from a generic instance to either other
  --   random flags or to specific flags entered by the user.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution planes;
  --   rows     row positions for white checkers
  --   cols     columns of white checkers of resolved condition;
  --   cnds     conditions kept fixed during flag continuation;
  --   inpt     if true, then user will give fixed flags,
  --            if false, then user will generate random flags.

  procedure Solve_Schubert_Problems ( n : in integer32 );

  -- DESCRIPTION :
  --   Interactive procedure to compute solutions to Schubert problems
  --   in n-space.  Prompts the user for data.

end Drivers_for_Schubert_Induction;
