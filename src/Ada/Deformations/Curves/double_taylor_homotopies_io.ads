with text_io;                            use text_io;
with Double_Taylor_Homotopies;           use Double_Taylor_Homotopies;

package Double_Taylor_Homotopies_io is

-- DESCRIPTION :
--   Offers output procedures for double Taylor homotopies,
--   mainly for testing purposes.

  procedure put ( tm : in Taylor_Monomial );
  procedure put ( file : in file_type; tm : in Taylor_Monomial );
  procedure put ( tm : in Link_to_Taylor_Monomial );
  procedure put ( file : in file_type; tm : in Link_to_Taylor_Monomial );

  -- DESCRIPTION :
  --   Writes the Taylor monomial tm to file or to standard output.

  procedure put ( tmv : in Taylor_Monomial_Vector );
  procedure put ( file : in file_type; tmv : in Taylor_Monomial_Vector );
  procedure put ( tmv : in Link_to_Taylor_Monomial_Vector );
  procedure put ( file : in file_type;
                  tmv : in Link_to_Taylor_Monomial_Vector );

  -- DESCRIPTION :
  --   Writes the Taylor monomial vector tmv to file or to standard output.

  procedure put ( th : in Taylor_Homotopy );
  procedure put ( file : in file_type; th : in Taylor_Homotopy );
  procedure put ( th : in Link_to_Taylor_Homotopy );
  procedure put ( file : in file_type; th : in Link_to_Taylor_Homotopy );

  -- DESCRIPTION :
  --   Writes the Taylor homotopy th to file or to standard output.

end Double_Taylor_Homotopies_io;
