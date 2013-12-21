with text_io;                            use text_io;
with Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Jumpstart_Diagonal_Homotopies is

-- DESCRIPTION :
--   This package provides jumpstarting of diagonal homotopies.

  function Ask_for_Monitor return boolean;

  -- DESCRIPTION :
  --   Asks if the user wants to monitor the progress of the path tracking
  --   on screen.  On return is true if the user does, false otherwise.

  procedure Read_Degree_of_Witness_Set
               ( file : in file_type; d,n : out natural32 );

  -- DESCRIPTION :
  --   Scans the input file for the solutions banner and then reads
  --   the number of solutions and the length of the solution vectors,
  --   returned respectively in the values of the parameters d and n.

  procedure Read_Witness_Set 
               ( file : out file_type; lp : out Link_to_Poly_Sys;
                 n,k,d : out natural32 );

  -- DESCRIPTION :
  --   Prompts the user for an input file which is then opened for input.
  --   On return is the embedded system and the relevant numerical info:
  --   ambient dimension n, dimension k and degree d of the witness set.
  --   The embedded slack variables are permuted to the end.

  procedure Read_Two_Witness_Sets 
               ( wf1,wf2 : out file_type; lp1,lp2 : out Link_to_Poly_Sys;
                 n1,n2,k1,k2,d1,d2 : out natural32;
                 sym1,sym2 : out Symbol_Table.Link_to_Array_of_Symbols );

  -- DESCRIPTION :
  --   Prompts the user to provide two input files for two witness sets.
  --   Returns the symbols used in both witness sets,
  --   needed to compute the matching permutation.

  -- ON RETURN :
  --   wf1,wf2   input files for the two witness sets;
  --   lp1,lp2   embedded systems;
  --   n1,n2     total number of variables and equations in the systems;
  --   k1,k2     dimension of the two solution sets;
  --   d1,d2     degrees of the two solution sets;
  --   sym1      symbols used in the first witness set;
  --   sym2      symbols used in the second witness set.

  -- REQUIRED :
  --   The original ambient dimension n equals n1 - k1 = n2 - k2.

  procedure Match_Symbols
               ( n1,n2,k1,k2 : in natural32; p2 : in out Poly_Sys;
                 ls1,ls2 : in Symbol_Table.Link_to_Array_of_Symbols;
                 solsym : out Symbol_Table.Link_to_Array_of_Symbols );

  -- DESCRIPTION :
  --   Computes the matching permutation between the two symbol arrays,
  --   used to read in the two polynomial systems.
  --   Permutes the second polynomial systems so that the symbols match
  --   and reassigns the symbol table.

  -- ON ENTRY :
  --   n1        number of variables and equations in the first system;
  --   n2        number of variables and equations in the second system;
  --   k1        dimension of the first witness set;
  --   k2        dimension of the second witness set;
  --   p2        embedded polynomial system for the second witness set;
  --   ls1       symbol table for the first witness set;
  --   ls2       symbol table for the second witness set.

  -- ON RETURN :
  --   p2        permuted 2nd system with symbols matching the 1st system;
  --   solsym    symbol table to read all witness points from file.

  -- REQUIRED :
  --   The original ambient dimension n equals n1 - k1 = n2 = k2.

  procedure Reset_to_Reread_Solutions
               ( file : in out file_type; d,n : out natural32 );

  -- DESCRIPTION :
  --   Resets the input file to reread the solutions,
  --   i.e.: scans the banner and the solution dimensions,
  --   so that a Read_Next invoked after this reset gives
  --   the first solution.  The length of the solution list and
  --   the dimension of each solution vector are respectively
  --   in the values of d and n on return.

  procedure Start_Cascade
               ( file,a_f : in file_type; b_f : in out file_type;
                 report : in boolean; start,target : in Poly_Sys;
                 a_n,b_n,a_dim,b_dim,a_deg,b_deg : in natural32;
                 sbt : in Symbol_Table.Link_to_Array_of_Symbols );

  -- DESCRIPTION :
  --   Uses jumpstarting in the homotopy to start the cascade.

  -- ON ENTRY :
  --   file      output file to write the results;
  --   a_f       file for the witness points of the 1st set;
  --   b_f       file for the witness points of the 2nd set;
  --   report    true if the user will monitor progress on screen;
  --   start     start system in the diagonal cascade;
  --   target    target system in the diagonal cascade;
  --   a_n       ambient dimension of the first system;
  --   b_n       ambient dimension of the second system;
  --   a_dim     dimension of the first algebraic set;
  --   b_dim     dimension of the second algebraic set;
  --   a_deg     degree of the first algebraic set;
  --   b_deg     degree of the second algebraic set;
  --   sbt       symbol table to use for reading witness points.

  -- ON RETURN :
  --   b_f       may have been reset to reread solutions.

  procedure Intersect_Witness_Sets
               ( file,a_f : in file_type;
                 b_f : in out file_type; a_p,b_p : in Poly_Sys;
                 a_n,b_n,a_dim,b_dim,a_deg,b_deg : in natural32;
                 s1e,s2e,sbt : in Symbol_Table.Link_to_Array_of_Symbols );

  -- DESCRIPTION :
  --   Uses an extrinsic diagonal homotopy to intersect two witness sets.

  -- REQUIRED : a_dim >= b_dim.

  -- ON ENTRY :
  --   file      output file;
  --   a_f       input file for the first witness set;
  --   b_f       input file for the second witness set;
  --   a_p       embedded system for the first witness set;
  --   b_p       embedded system for the second witness set;
  --   a_n       #variables and equations in the first system;
  --   b_n       #variables and equations in the second system;
  --   a_dim     dimension of the first witness set;
  --   b_dim     dimension of the second witness set;
  --   a_deg     degree of the first witness set;
  --   b_deg     degree of the second witness set;
  --   s1e       symbols used to read in first witness set;
  --   s2e       symbols used to read in second witness set;
  --   sbt       symbol table to use to read all witness points.

  -- ON RETURN :
  --   b_f       is "in out" for resetting to reread solutions.

  procedure Jumpstart_Diagonal_Homotopy;

  -- DESCRIPTION :
  --   Prompts the user for two witness sets and runs the first homotopy
  --   to start the cascade.

  procedure Jumpstart_Diagonal_Homotopy
               ( infile : in out file_type;
	         outfile : in file_type; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   Given the input file and target system for the first witness set,
  --   this procedure will prompt the user for another file name
  --   for the second witness set and then run the first homotopy
  --   to start the cascade.

  -- ON ENTRY :
  --   infile    input file which contains the first witness set,
  --             positioned after reading the polynomial system p;
  --   outfile   output file for the results;
  --   p         the polynomial system read from infile.

  -- ON RETURN :
  --   infile    may be reset for input to reread solution list;
  --   p         slack variables have been permuted to the end.

  procedure Track_Paths_Down
               ( infile,outfile : in file_type;
                 start : in Poly_Sys; n,k,d : in natural32 );

  -- DESCRIPTION :
  --   This procedure is called by Jumpstart_Cascade_Homotopy.

  procedure Jumpstart_Cascade_Homotopy;

  -- DESCRIPTION :
  --   Prompts the user to enter an embedded system
  --   and removes the last slice.

  procedure Jumpstart_Cascade_Homotopy
               ( infile,outfile : in file_type; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   Uses the embedded system p and the solutions from infile
  --   to go one level down in a homotopy cascade.

  -- ON ENTRY :
  --   infile    input file which contains the first witness set,
  --             positioned after reading the polynomial system p;
  --   outfile   output file for the results;
  --   p         the polynomial system read from infile.

  -- ON RETURN :
  --   p         slack variables have been permuted to the end.

  procedure Remove_Last_Slack_Variable
               ( infile,outfile : in file_type; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   Removes the last slack variable from the embedded system p
  --   and from the solutions from infile, as needed after going
  --   one level down in the cascade.

  -- ON ENTRY :
  --   infile    input file which contains one witness set,
  --             where the last slack variable has been set to zero;
  --   outfile   output file to write the new witness set,
  --             after removal of the last slack variable;
  --   p         the polynomial system read from input file infile.

  -- ON RETURN :
  --   p         last slack variable has been removed.

end Jumpstart_Diagonal_Homotopies;
