with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with QuadDobl_Complex_Numbers;         use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Solutions;       use QuadDobl_Complex_Solutions;

package QuadDobl_Solution_Filters is

-- DESCRIPTION :
--   This package offers tools to create new solution lists by
--   filtering lists subject to certain criteria.

-- CRITERIA :

  function Vanishing ( sol : Solution; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if either the error on the last correction term
  --   or the residual are smaller than the given tolerance.

  function Zero_Component ( sol : Solution; k : natural32;
                            tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   The k-th component of a solution is considered as zero if its
  --   absolute value is less than the given tolerance.

  function Regular ( sol : Solution; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   A solution is regular if the estimate for the inverse condition
  --   number is larger than the given tolerance.

  function On_Target ( sol : Solution; target : Complex_Number;
                       tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if | sol.t - target | < tol.

-- PROTOTYPE FILTERS :

  generic

    with function Criterium ( i : natural32; s : Solution ) return boolean;

    -- DESCRIPTION :
    --   Returns true (or false) if the i-th solution s fits (or not)
    --   any particular criterium implemented by this function.

  function List_Filter ( sols : Solution_List ) return Solution_List;

  -- DESCRIPTION :
  --   Returns a solution list, filtered by the given criterium.
  --   The list on return shares no data with the given list.

  generic
  
    with function Criterium ( i : natural32; s : Solution ) return boolean;

    -- DESCRIPTION :
    --   Returns true (or false) if the i-th solution s fits (or not)
    --   any particular criterium implemented by this function.

  procedure Scan_Filter
              ( infile,outfile : in file_type; len,dim : in natural32;
                cnt : out natural32 );

  -- DESCRIPTION :
  --   Filters a solution lists of the given length, reading from infile
  --   and writing to outfile, subject to the given criterium.

  -- REQUIRED :
  --   The position at infile is ready for reading the next solution.

  -- ON ENTRY :
  --   infile   input file where the solution list is;
  --   outfile  output file, must be opened for output;
  --   len      length of the solution list on infile;
  --   dim      dimension of the solution vectors on infile.

  -- ON RETURN :
  --   cnt      number of solutions written to outfile.

  generic
  
    with function Criterium ( i : natural32; s : Solution ) return boolean;

    -- DESCRIPTION :
    --   Returns true (or false) if the i-th solution s fits (or not)
    --   any particular criterium implemented by this function.

  procedure Scan_Filter1
              ( infile,outfile : in file_type; len,dim : in natural32;
                ls : in Link_to_Solution; cnt : out natural32 );

  -- DESCRIPTION :
  --   Filters a solution lists of the given length, reading from infile
  --   and writing to outfile, subject to the given criterium, with the
  --   first solution already given.

  -- REQUIRED :
  --   The position at infile is ready for reading the next solution.

  -- ON ENTRY :
  --   infile   input file where the solution list is;
  --   outfile  output file, must be opened for output;
  --   len      length of the solution list on infile;
  --   dim      dimension of the solution vectors on infile;
  --   ls       pointer to the first solution already read from file.

  -- ON RETURN :
  --   cnt      number of solutions written to outfile.

-- FILTERS :
--   functions follow the "List_Filter" prototype, while
--   procedures are along the "Scan_Filter" prototype.

  function Vanishing_Filter 
             ( sols : Solution_List; tol : double_float )
             return Solution_List;
  procedure Vanishing_Filter 
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that are considered as vanishing
  --   with respect to the given tolerance on errors and residuals.

  function Spurious_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List;
  procedure Spurious_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that are considered as spurious
  --   that is nonvanishing for the given tolerance.

  function Zero_Component_Filter
             ( sols : Solution_List; k : natural32; tol : double_float )
             return Solution_List;
  procedure Zero_Component_Filter
             ( infile,outfile : in file_type; len,dim,k : in natural32;
               tol : in double_float; ls : in Link_to_Solution;
               cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that are considered as having their
  --   k-th component zero with respect to the given tolerance.
  --   The procedure follows the Scan_Filter1 prototype.

  function Free_Component_Filter
             ( sols : Solution_List; k : natural32; tol : double_float )
             return Solution_List;
  procedure Free_Component_Filter
             ( infile,outfile : in file_type; len,dim,k : natural32;
               tol : in double_float; ls : in Link_to_Solution;
               cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that are considered as having their
  --   k-th component free (nonzero) with respect to the given tolerance.
  --   The procedure follows the Scan_Filter1 prototype.

  function Regular_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List;
  procedure Regular_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that are considered as regular with
  --   respect to the given tolerance for the inverse condition number.

  function Singular_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List;
  procedure Singular_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that are considered as singular with
  --   respect to the given tolerance for the inverse condition number.

  function On_Target_Filter
             ( sols : Solution_List; target : Complex_Number;
               tol : double_float ) return Solution_List;
  procedure On_Target_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               target : in Complex_Number; tol : in double_float;
               cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that have reached the given target.

  function Off_Target_Filter 
             ( sols : Solution_List; target : Complex_Number;
               tol : double_float ) return Solution_List;
  procedure Off_Target_Filter 
             ( infile,outfile : in file_type; len,dim : in natural32;
               target : in Complex_Number; tol : in double_float;
               cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions that have not reached the given target.

  function Real_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List;
  procedure Real_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 );

  -- DESCRIPTION :
  --   Returns the list of solutions whose components all have an
  --   imaginary part smaller than the given tolerance.

  function Select_Solutions
             ( sols : Solution_List;
               nb : Standard_Natural_Vectors.Vector ) return Solution_List;
  procedure Select_Solutions
             ( infile,outfile : in file_type; len,dim : in natural32;
               nb : in Standard_Natural_Vectors.Vector; cnt : out natural32 );

  -- DESCRIPTION :
  --   Selects the solutions according to the numbers in the vector nb.

  -- REQUIRED : the entries in nb are sorted in increasing order.

  function Select_Failed_Solutions
              ( psols,qsols : Solution_List; tol : double_float;
                verbose : boolean := false ) return Solution_List;

  -- DESCRIPTION :
  --   Given a list of solutions at the end of paths in psols,
  --   returns the corresponding start solutions in qsols,
  --   relative to the tolerance tol.
  --   If verbose, then at every failed path, one line is written to screen.

  -- REQUIRED : Length_Of(psols) = Length_Of(qsols).

  -- ON ENTRY :
  --   psols    solutions at the end of the paths;
  --   qsols    start solutions corresponding to psols;
  --   tol      tolerance using in functions On_Target and Vanishing;
  --   verbose  if true, then extra output is written for every failed path.

end QuadDobl_Solution_Filters;
