with text_io;                            use text_io;
with Multprec_Complex_Term_Lists;        use Multprec_Complex_Term_Lists;

package Multprec_Complex_Term_Lists_io is

  procedure put ( p : in Term_List );
  procedure put ( file : in file_type; p : in Term_List );

  -- DESCRIPTION :
  --   Writes the terms in the list p in tableau format:
  --   every term in p occupies exactly one line,
  --   the real and imaginary parts of the complex coefficient 
  --   come first, written plainly as two multiprecision floating point
  --   numbers, followed by the integer exponents.

end Multprec_Complex_Term_Lists_io;
