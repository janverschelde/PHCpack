with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Witness_Interface is

-- DESCRIPTION :
--   Provides functions interface with witness sets for polynomial and
--   Laurent systems in double, double double, and quad double precision.

  function Witness_Standard_Polynomial_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system witness set
  --   in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the total number of equations;
  --   b       in b[0] is the dimension of the witness set,
  --           in b[1] is the degree of the witness set.

  function Witness_DoblDobl_Polynomial_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system witness set
  --   in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the total number of equations;
  --   b       in b[0] is the dimension of the witness set,
  --           in b[1] is the degree of the witness set.

  function Witness_QuadDobl_Polynomial_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system witness set
  --   in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the total number of equations;
  --   b       in b[0] is the dimension of the witness set,
  --           in b[1] is the degree of the witness set.

  function Witness_Standard_Laurent_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system witness set
  --   in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the total number of equations;
  --   b       in b[0] is the dimension of the witness set,
  --           in b[1] is the degree of the witness set.

  function Witness_DoblDobl_Laurent_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system witness set
  --   in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the total number of equations;
  --   b       in b[0] is the dimension of the witness set,
  --           in b[1] is the degree of the witness set.

  function Witness_QuadDobl_Laurent_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system witness set
  --   in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the total number of equations;
  --   b       in b[0] is the dimension of the witness set,
  --           in b[1] is the degree of the witness set.

  function Witness_Standard_Polynomial_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a polynomial system witness set in double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       the file name for the witness set;
  --   vrblvl  is the verbose level.

  function Witness_DoblDobl_Polynomial_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a polynomial system witness set 
  --   in double double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       the file name for the witness set;
  --   vrblvl  is the verbose level.

  function Witness_QuadDobl_Polynomial_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a polynomial system witness set 
  --   in quad double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       the file name for the witness set;
  --   vrblvl  is the verbose level.

  function Witness_Standard_Laurent_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a Laurent system witness set 
  --   in double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       the file name for the witness set;
  --   vrblvl  is the verbose level.

  function Witness_DoblDobl_Laurent_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a Laurent system witness set 
  --   in double double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       the file name for the witness set;
  --   vrblvl  is the verbose level.

  function Witness_QuadDobl_Laurent_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a Laurent system witness set 
  --   in quad double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       the file name for the witness set;
  --   vrblvl  is the verbose level.

  function Witness_Standard_Polynomial_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps the variables in the polynomial system and solutions stored
  --   in double precision so the slack variables are last.

  -- ON ENTRY :
  --   a       in a[0] are the number of variables;
  --   b       in b[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Witness_DoblDobl_Polynomial_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps the variables in the polynomial system and solutions stored
  --   in double double precision so the slack variables are last.

  -- ON ENTRY :
  --   a       in a[0] are the number of variables;
  --   b       in b[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Witness_QuadDobl_Polynomial_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps the variables in the polynomial system and solutions stored
  --   in quad double precision so the slack variables are last.

  -- ON ENTRY :
  --   a       in a[0] are the number of variables;
  --   b       in b[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Witness_Standard_Laurent_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps the variables in the Laurent system and solutions stored
  --   in double precision so the slack variables are last.

  -- ON ENTRY :
  --   a       in a[0] are the number of variables;
  --   b       in b[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Witness_DoblDobl_Laurent_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps the variables in the Laurent system and solutions stored
  --   in double double precision so the slack variables are last.

  -- ON ENTRY :
  --   a       in a[0] are the number of variables;
  --   b       in b[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Witness_QuadDobl_Laurent_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps the variables in the Laurent system and solutions stored
  --   in quad double precision so the slack variables are last.

  -- ON ENTRY :
  --   a       in a[0] are the number of variables;
  --   b       in b[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Witness_Standard_Polynomial_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Embeds the polynomial system in double precision 
  --   adding as many slack variables as the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables to add;
  --   vrblvl  is the verbose level.

  function Witness_DoblDobl_Polynomial_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Embeds the polynomial system in double precision 
  --   adding as many slack variables as the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables to add;
  --   vrblvl  is the verbose level.

  function Witness_QuadDobl_Polynomial_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Embeds the polynomial system in quad double precision 
  --   adding as many slack variables as the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables to add;
  --   vrblvl  is the verbose level.

  function Witness_Standard_Laurent_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Embeds the Laurent system in double precision 
  --   adding as many slack variables as the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables to add;
  --   vrblvl  is the verbose level.

  function Witness_DoblDobl_Laurent_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Embeds the Laurent system in double precision 
  --   adding as many slack variables as the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables to add;
  --   vrblvl  is the verbose level.

  function Witness_QuadDobl_Laurent_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Embeds the Laurent system in double precision 
  --   adding as many slack variables as the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables to add;
  --   vrblvl  is the verbose level.

  function Witness_Standard_Polynomial_Write
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes a witness set for a polynomial system in double precision
  --   to the file with the given name.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the file name;
  --   b       contains the characters of the file name;
  --   vrblvl  is the verbose level.

  function Witness_DoblDobl_Polynomial_Write
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes a witness set for a polynomial system in double double
  --   precision to the file with the given name.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the file name;
  --   b       contains the characters of the file name;
  --   vrblvl  is the verbose level.

  function Witness_QuadDobl_Polynomial_Write
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes a witness set for a polynomial system in quad double
  --   precision to the file with the given name.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the file name;
  --   b       contains the characters of the file name;
  --   vrblvl  is the verbose level.

end Witness_Interface;
