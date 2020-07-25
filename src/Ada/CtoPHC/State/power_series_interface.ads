with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Power_Series_Interface is

-- DESCRIPTION :
--   Functions in this package define the interface to the power
--   series methods to compute Pade approximants for solution curves.

  function Series_Standard_Newton_at_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Newton's method in double precision starting at a point.
  --   The start point is a solution vector in the container.
  --   The series are stored in the systems pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter,
  --           in a[1] is the maximal degree of the series,
  --           in a[2] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_DoblDobl_Newton_at_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Newton's method in double double precision starting at a point.
  --   The start point is a solution vector in the container.
  --   The series are stored in the systems pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter,
  --           in a[1] is the maximal degree of the series,
  --           in a[2] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_QuadDobl_Newton_at_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Newton's method in quad double precision starting at a point.
  --   The start point is a solution vector in the container.
  --   The series are stored in the systems pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter,
  --           in a[1] is the maximal degree of the series,
  --           in a[2] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_Standard_Newton_at_Series
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Newton's method in double precision starting at a series.
  --   The start series are assumed to be in the systems pool container.
  --   The series are stored in the systems pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter,
  --           in a[1] is the maximal degree of the series,
  --           in a[2] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_DoblDobl_Newton_at_Series
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Newton's method in double double precision starting at a series.
  --   The start series are assumed to be in the systems pool container.
  --   The series are stored in the systems pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter,
  --           in a[1] is the maximal degree of the series,
  --           in a[2] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_QuadDobl_Newton_at_Series
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Newton's method in quad double precision starting at a series.
  --   The start series are assumed to be in the systems pool container.
  --   The series are stored in the systems pool.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter,
  --           in a[1] is the maximal degree of the series,
  --           in a[2] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_Standard_Pade
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Constructs a Pade approximant in double precision.
  --   Results are in the pool for double precision systems.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter;
  --           in a[1] is the degree of the numerator;
  --           in a[2] is the degree of the denominator;
  --           in a[3] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_DoblDobl_Pade
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Constructs a Pade approximant in double double precision.
  --   Results are in the pool for double double precision systems.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter;
  --           in a[1] is the degree of the numerator;
  --           in a[2] is the degree of the denominator;
  --           in a[3] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

  function Series_QuadDobl_Pade
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Constructs a Pade approximant in quad double precision.
  --   Results are in the pool for quad double precision systems.

  -- ON ENTRY :
  --   a       in a[0] is the index of the series parameter;
  --           in a[1] is the degree of the numerator;
  --           in a[2] is the degree of the denominator;
  --           in a[3] is the number of steps in Newton's method;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then intermediate results are written to screen;
  --   vrblvl  is the verbose level.

end Power_Series_Interface; 
