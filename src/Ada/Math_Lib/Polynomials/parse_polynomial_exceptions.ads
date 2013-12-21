package Parse_Polynomial_Exceptions is

-- DESCRIPTION :
--   This package defined the exceptional situations which may
--   occur when reading a polynomial in the wrong format.

  ILLEGAL_CHARACTER : exception;
      -- occurs when an unexpected character is found

  ILLEGAL_OPERATION : exception;
      -- occurs when illegal operations with polynomials is attempted

  OVERFLOW_OF_UNKNOWNS : exception;
      -- occurs when the number of unknowns exceeds number in symbol table

  BAD_BRACKET : exception;
      -- occurs when a bracket, like '(' or ')' is misplaced

end Parse_Polynomial_Exceptions;
