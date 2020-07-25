with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Monomial_Maps_Interface is

-- DESCRIPTION :
--   Solutions to polynomial systems defined by binomial equations
--   are monomial maps and can be computed efficiently.

  function Monomial_Maps_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the blackbox solver for binomial systems
  --   on the Laurent system stored with double precision coefficients.

  -- ON ENTRY :
  --   a       in a[0] is the option puretopdim,
  --           to indicate that the ideal is pure dimensional;
  --   vrblvl  is the verbose level.

  function Monomial_Maps_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the maps in the container to screen.
  --   The input parameter is the verbose level.

  function Monomial_Maps_Top_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the top dimension of the monomial maps.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the top dimension of the maps.

  function Monomial_Maps_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of maps of the given dimension.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of interest;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of maps of dimension a[0].

  function Monomial_Maps_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of a map.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the map,
  --           in a[1] is the index of the map;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the degree of the map.

  function Monomial_Maps_Coefficients
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the coefficients of a map.

  function Monomial_Maps_Exponents
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
   --   Returns the exponents of a map.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the map,
  --           in a[1] is the index of the map;
  --           in a[2] is the number of variables;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b are the exponents of the map,
  --           as one flattened vector of dimension many numbers
  --           times the number of variables.

  function Monomial_Maps_Data
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Return coefficients and exponents of a map.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the map,
  --           in a[1] is the index of the map;
  --           in a[2] is the number of variables;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b are the exponents of the map,
  --           as one flattened vector of dimension many numbers
  --           times the number of variables;
  --   c       in c are the coefficient of the map,
  --           as one flattened vector of real and imaginary parts
  --           as many doubles as 2 times the number of variables.

  function Monomial_Maps_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the maps.
  --   The verbose level in in vrblvl.

end Monomial_Maps_Interface;
