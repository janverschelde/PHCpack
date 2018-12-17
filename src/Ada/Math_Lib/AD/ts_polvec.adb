with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Polynomial_Vectors;
with Standard_Polynomial_Vectors_io;    use Standard_Polynomial_Vectors_io;
with DoblDobl_Polynomial_Vectors;
with DoblDobl_Polynomial_Vectors_io;    use DoblDobl_Polynomial_Vectors_io;
with QuadDobl_Polynomial_Vectors;
with QuadDobl_Polynomial_Vectors_io;    use QuadDobl_Polynomial_Vectors_io;
with Random_Polynomial_Vectors;         use Random_Polynomial_Vectors;

procedure ts_polvec is

-- DESCRIPTION :
--   Tests the operations on monomials.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random polynomial vector of monomials with coefficients
  --   in standard double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    p : Standard_Polynomial_Vectors.Link_to_Polynomial_Vector;

  begin
    put_line("Testing monomial operations in double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant Standard_Polynomial_Vectors.Polynomial_Vector
             := Standard_Random_Polynomial_Vector(size,dim,expmax,true);
    begin
      p := new Standard_Polynomial_Vectors.Polynomial_Vector'(ranvec);
      put_line("a random polynomial vector : "); put(p);
    end;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random polynomial vector of monomials with coefficients
  --   in double double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    p : DoblDobl_Polynomial_Vectors.Link_to_Polynomial_Vector;

  begin
    put_line("Testing monomial operations in double double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant DoblDobl_Polynomial_Vectors.Polynomial_Vector
             := DoblDobl_Random_Polynomial_Vector(size,dim,expmax,true);
    begin
      p := new DoblDobl_Polynomial_Vectors.Polynomial_Vector'(ranvec);
      put_line("a random polynomial vector : "); put(p);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random polynomial vector of monomials with coefficients
  --   in quad double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    p : QuadDobl_Polynomial_Vectors.Link_to_Polynomial_Vector;

  begin
    put_line("Testing monomial operations in quad double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant QuadDobl_Polynomial_Vectors.Polynomial_Vector
             := QuadDobl_Random_Polynomial_Vector(size,dim,expmax,true);
    begin
      p := new QuadDobl_Polynomial_Vectors.Polynomial_Vector'(ranvec);
      put_line("a random monomial vector :"); put(p);
    end;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then lauches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_polvec;
