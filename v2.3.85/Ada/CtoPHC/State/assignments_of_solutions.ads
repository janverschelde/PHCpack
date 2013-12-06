with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Assignments_of_Solutions is

-- DESCRIPTION :
--   The package provides routines to convert the (b,c)-format for
--   a solution into its PHCpack format, and vice versa.

  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return Standard_Complex_Solutions.Solution;
  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return DoblDobl_Complex_Solutions.Solution;
  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return QuadDobl_Complex_Solutions.Solution;
  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return Standard_Complex_Solutions.Link_to_Solution;
  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return DoblDobl_Complex_Solutions.Link_to_Solution;
  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return QuadDobl_Complex_Solutions.Link_to_Solution;

  -- DESCRIPTION :
  --   Converts the input parameters in b and in c into a solution type.

  -- REQUIRED :
  --   for standard solutions of dimension n, c holds 2*n + 5 doubles,
  --   for double double solutions of dimension n, c holds 4*n + 10 doubles,
  --   for quad double solutions of dimension n, c holds 8*n + 20 doubles.

  -- ON ENTRY :
  --   b       pointer to the dimension and multiplicity:
  --           b[0] has the dimension n of the solution vector,
  --           b[1] contains the multiplicity flag of a solution;
  --   c       real and imaginary part of the continuation parameter,
  --           the real and imaginary parts of the solution coordinates,
  --           and the diagnostics (err,rco,res) as last 3 doubles.

  procedure Assign_Solution
             ( s : in Standard_Complex_Solutions.Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer );
  procedure Assign_Solution
             ( s : in DoblDobl_Complex_Solutions.Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer );
  procedure Assign_Solution
             ( s : in QuadDobl_Complex_Solutions.Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer );
  procedure Assign_Solution
             ( ls : in Standard_Complex_Solutions.Link_to_Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer );
  procedure Assign_Solution
             ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer );
  procedure Assign_Solution
             ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer );

  -- DESCRIPTION :
  --   Assign the data in s and ls to the output parameters b and c.

  -- REQUIRED : c has space for 2*n + 5 doubles for standard solutions,
  --   or 4*n + 10 for double double solutions,
  --   or 8*n + 20 for quad double solutions.

  -- ON ENTRY :
  --   ls      pointer to a solution in PHCpack format;
  --   b       address to multiplicity field;
  --   c       address to the coordinates of t, solution vector,
  --           and the diagnostics (for more details, see above).

  -- ON RETURN :
  --   b       contains multiplicity of the solution;
  --   c       real and imaginary parts of the continuation parameter,
  --           the real and imaginary parts of the solution coordinates,
  --           and the diagnostics (err,rco,res) as last 3 doubles.

end Assignments_of_Solutions;
