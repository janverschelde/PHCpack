with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;
with QuadDobl_Point_Lists;               use QuadDobl_Point_Lists;
with QuadDobl_Quad_Trees;                use QuadDobl_Quad_Trees;

package QuadDobl_Condition_Report is

-- DESCRIPTION :
--   The condition report on a computed list of isolated solutions,
--   computed with double double floating-point arithmetic consists
--   of the (err,rco,res) fields as reported by running Newton's method
--   with as well a comparison to other solutions in the list to see if
--   all solutions are distinct or if there are clusters of solutions.

  procedure Write_Diagnostics ( sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes all diagnostics of the solution list to the screen.

  procedure Compute_Condition_Tables ( sols : in Solution_List );

  -- DESCRIPTION :
  --   Computes the condition tables for the list of solutions
  --   and writes the condition tables directly to screen.

  procedure Write_Condition_Results
               ( file : file_type;
                 e,c,r : in Standard_Natural_Vectors.Vector;
                 cnt_real : in natural32; tol_real : in double_float );

  -- DESCRIPTION :
  --   Writes the results of the frequency analysis to file.

  -- ON ENTRY :
  --   e         frequency table for norm of last correction vector;
  --   c         frequency table for condition numbers;
  --   r         frequency table for norm of residual vector;
  --   cnt_real  number of real solutions found;
  --   tol_real  tolerance to decide whether a solution is real.

  procedure Write_Condition_Results
               ( file : file_type; i,f : in natural32;
                 e,c,r : in Standard_Natural_Vectors.Vector;
                 cnt_real : in natural32; tol_real : in double_float );

  -- DESCRIPTION :
  --   Writes the results of the frequency analysis to file.

  -- ON ENTRY :
  --   i         index to the current solution;
  --   f         frequency updater of solution numbers;
  --   e         frequency table for norm of last correction vector;
  --   c         frequency table for condition numbers;
  --   r         frequency table for norm of residual vector;
  --   cnt_real  number of real solutions found;
  --   tol_real  tolerance to decide whether a solution is real.

  procedure Count_Clusters
               ( infile : in out file_type; outfile : in file_type;
                 bannered,to_file : in boolean;
                 root : in Link_to_Quad_Node; tol : in double_float );

  -- DESCRIPTION :
  --   Counts the clusters and presents the user with a menu
  --   to handle the candidate clustered pairs.

  -- ON ENTRY :
  --   infile    file where the solutions are, may have to be reset;
  --   outfile   file for output;
  --   bannered  indicates if system preceeds the solutions;
  --   to_file   if the cluster report must be written to file;
  --   root      root of a quad tree;
  --   tol       tolerance to decide clustering.

  -- ON RETURN :
  --   infile    may have been reset to write clustered solutions.

  procedure Write_Cluster_Report
               ( infile : in out file_type; outfile : in file_type;
                 bannered,to_file : in boolean;
                 pl : in Point_List; size : in natural32;
                 tol : in double_float );

  -- DESCRIPTION :
  --   Creates a quad tree for the list of points in pl and uses this
  --   to search for pairs of clustered solutions.

  -- ON ENTRY :
  --   infile    input file where the solutions are;
  --   outfile   file to write the results to;
  --   bannered  true if the solutions on the input file are preceeded
  --             by a banner, false otherwise,
  --   to_file   true if the results need to go to file;
  --   pl        projected list of points;
  --   size      number of points in pl;
  --   tol       tolerance for clustering.

  procedure Write_Cluster_Report
               ( outfile : in file_type; to_file : in boolean;
                 sols : in Solution_List;
                 pl : in Point_List; tol : in double_float );

  -- DESCRIPTION :
  --   Creates a quad tree for the list of points in pl and uses this
  --   to search for pairs of clustered solutions.

  -- ON ENTRY :
  --   outfile   file to write the results to;
  --   to_file   true if the results need to go to file;
  --   sols      list of solutions;
  --   pl        projected list of points;
  --   tol       tolerance for clustering.

  procedure Is_Clustered
               ( s : in Solution; nb : in natural32;
                 sols : in Solution_List; tol : in double_float; 
                 h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 );
  procedure Is_Clustered
               ( s : in Solution; nb : in natural32;
                 sols : in Solution_Array; tol : in double_float; 
                 h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 );

  -- DESCRIPTION :
  --   Compared to the function Standard_Solution_Diagnostics.Is_Clustered,
  --   this procedure scales better for larger solution lists.

  -- ON ENTRY :
  --   s         a solution which occurs in sols at position nb;
  --   nb        position of the solution s in the list or array sols;
  --   sols      list or array of solutions;
  --   tol       tolerance to decide whether two solutions are clustered.
  --   h1        first hash key for the point list;
  --   h2        second has key for the point list;
  --   pl        list of points updated up to position nb-1.
  
  -- ON RETURN :
  --   pl        the hashed version of the solution has been inserted
  --             into the point list pl;
  --   val       if equal to nb, then the solution is not clustered,
  --             else, val is the index of the first other occurrence
  --             of the solution s in sols.

  procedure Multiplicity
               ( s : in out Solution; nb : in natural32;
                 sols : in Solution_List; tol : in double_float; 
                 h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 );
  procedure Multiplicity
               ( s : in out Solution; nb : in natural32;
                 sols : in Solution_Array; tol : in double_float; 
                 h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 );

  -- DESCRIPTION :
  --   Compared to the function Standard_Solution_Diagnostics.Multiplicity,
  --   this procedure scales better for larger solution lists.

  -- ON ENTRY :
  --   s         a solution which occurs in sols at position nb;
  --   nb        position of the solution s in the list or array sols;
  --   sols      list or array of solutions;
  --   tol       tolerance to decide whether two solutions are clustered.
  --   h1        first hash key for the point list;
  --   h2        second has key for the point list;
  --   pl        list of points updated up to position nb-1.
  
  -- ON RETURN :
  --   s         multiplicity field of s may have been increased to val;
  --   sols      adjusted multiplicity field for clusters;
  --   pl        the hashed version of the solution has been inserted
  --             into the point list pl;
  --   val       if equal to 1, if the solution s is not clustered,
  --             else, val counts the number of occurrences of the
  --             solution s in sols.

  procedure Scan_for_Condition_Tables 
               ( infile : in out file_type; outfile : in file_type;
                 bannered,to_file : in boolean;
                 len,dim : in natural32;
                 tol_real,tol_clus : in double_float;
                 i_end,f_val : out natural32;
                 e,c,r : out Standard_Natural_Vectors.Vector;
                 nb_real : out natural32; pl : out Point_List );

  -- DESCRIPTION :
  --   Scans the input file for solutions to compute a condition table.
  --   For huge solution lists, this procedure will not store
  --   the solutions in main memory.

  -- REQUIRED : e'range = c'range = r'range = 0..15.

  -- ON ENTRY :
  --   infile    input file for the solutions, where the two numbers on
  --             top of the list (len and dim) have already been read, 
  --             the type is "in out" because it may have to be reset
  --             for reading again;
  --   outfile   output file to write the condition tables on;
  --   bannered  true if the solutions were preceeded by a banner;
  --   to_file   true if output file is not standard_output;
  --   len       length of the solution list on the input file;
  --   dim       dimension of the solution vectors;
  --   tol_real  tolerance to decide whether a solution is real or not;
  --   tol_clus  tolerance to decide whether two solutions are clustered.

  -- ON RETURN :
  --   i_end     value of index to current solution at the end;
  --   f_val     value of the frequency updater;
  --   e         frequency counts of the err values;
  --   c         frequency counts of the rco values;
  --   r         frequency counts of the res values;
  --   nb_real   equals the number of real solutions,
  --             with respect to the tolerance in tol_real;
  --   pl        projected list of points for cluster detection.

  procedure Scan_for_Condition_Tables 
               ( file : in file_type; sols : in Solution_List;
                 tol_real,tol_clus : in double_float;
                 e,c,r : out Standard_Natural_Vectors.Vector;
                 nb_real : out natural32; pl : out Point_List );

  -- DESCRIPTION :
  --   Scans the solutions to compute a condition table.
  --   Unlike the other routine, all solutions in sols are
  --   in main memory.

  -- REQUIRED : e'range = c'range = r'range = 0..15.

  -- ON ENTRY :
  --   file      output file to write the condition tables on;
  --   tol_real  tolerance to decide whether a solution is real or not;
  --   tol_clus  tolerance to decide whether two solutions are clustered.

  -- ON RETURN :
  --   e         frequency counts of the err values;
  --   c         frequency counts of the rco values;
  --   r         frequency counts of the res values;
  --   nb_real   equals the number of real solutions,
  --             with respect to the tolerance in tol_real;
  --   pl        projected list of points for cluster detection.

end QuadDobl_Condition_Report;
