with text_io;                            use text_io;
with DoblDobl_Complex_Term_Lists;        use DoblDobl_Complex_Term_Lists;

package DoblDobl_Complex_Term_Lists_io is

  procedure put ( p : in Term_List );
  procedure put ( file : in file_type; p : in Term_List );

  -- DESCRIPTION :
  --   Writes the terms in the list p in tableau format:
  --   every term in p occupies exactly one line,
  --   the real and imaginary parts of the complex coefficient 
  --   come first, written plainly as two double doubles,
  --   followed by the integer exponents.

  procedure put ( p : in Array_of_Term_Lists );
  procedure put ( file : in file_type; p : in Array_of_Term_Lists );

  -- DESCRIPTION :
  --   Writes the array of term lists stored in p to standard output
  --   or to file, using tableau format, with one blank line to
  --   separate the polynomials in p.

end DoblDobl_Complex_Term_Lists_io;
