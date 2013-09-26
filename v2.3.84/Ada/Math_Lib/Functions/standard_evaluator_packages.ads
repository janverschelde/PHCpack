with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Standard_Evaluator_Packages is

-- DESCRIPTION :
--   Provides creators of a package to evaluate a system and its Jacobian.

-- PRIMITIVE OPERATIONS :

  procedure Replace_Symbols;

  -- DESCRIPTION :
  --   Replaces all symbols in the symbol table with vector entries:
  --   x(1), x(2), up to x(n).

  procedure Create_Inline_System_Evaluator
               ( file : in file_type; funname : in String; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Writes the body of a function for an evaluator for p on file.
  --   The name of the function is parametrized by "funname".

  procedure Create_Inline_Jacobian_Evaluator
               ( file : in file_type; funname : in String; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Writes the body of a function to evaluate the Jacobian matrix of
  --   p on file.  The name of the function is parametrized by "funname".

  function Read_Package_Name return String;

  -- DESCRIPTION :
  --   Reads the package name from standard input and returns the string.

-- TARGET ROUTINES :

  procedure Create ( packname : in String; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Creates a package with name packname to evaluate p and its
  --   Jacobian matrix.

  procedure Create ( p : in Poly_Sys );

  -- DESCRIPTION :
  --   Creates a package to evaluate the system p and its Jacobian matrix.
  --   The package name will be read and the file will be created.

end Standard_Evaluator_Packages;
