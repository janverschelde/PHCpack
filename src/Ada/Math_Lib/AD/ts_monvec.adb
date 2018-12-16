with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Monomial_Vectors;
with Standard_Monomial_Vectors_io;      use Standard_Monomial_Vectors_io;
with DoblDobl_Monomial_Vectors;
with DoblDobl_Monomial_Vectors_io;      use DoblDobl_Monomial_Vectors_io;
with QuadDobl_Monomial_Vectors;
with QuadDobl_Monomial_Vectors_io;      use QuadDobl_Monomial_Vectors_io;
with Random_Monomial_Vectors;           use Random_Monomial_Vectors;

procedure ts_monvec is

-- DESCRIPTION :
--   Tests the operations on monomials.

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial vector of monomials with coefficients
  --   in standard double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    v : Standard_Monomial_Vectors.Link_to_Monomial_Vector;

  begin
    put_line("Testing monomial operations in double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant Standard_Monomial_Vectors.Monomial_Vector
             := Standard_Random_Monomial_Vector(size,dim,expmax,true);
    begin
      v := new Standard_Monomial_Vectors.Monomial_Vector'(ranvec);
      put_line("a random monomial vector : "); put(v);
    end;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial vector of monomials with coefficients
  --   in double double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    v : DoblDobl_Monomial_Vectors.Link_to_Monomial_Vector;

  begin
    put_line("Testing monomial operations in double double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant DoblDobl_Monomial_Vectors.Monomial_Vector
             := DoblDobl_Random_Monomial_Vector(size,dim,expmax,true);
    begin
      v := new DoblDobl_Monomial_Vectors.Monomial_Vector'(ranvec);
      put_line("a random monomial vector : "); put(v);
    end;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial vector of monomials with coefficients
  --   in quad double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    v : QuadDobl_Monomial_Vectors.Link_to_Monomial_Vector;

  begin
    put_line("Testing monomial operations in quad double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    declare
      ranvec : constant QuadDobl_Monomial_Vectors.Monomial_Vector
             := QuadDobl_Random_Monomial_Vector(size,dim,expmax,true);
    begin
      v := new QuadDobl_Monomial_Vectors.Monomial_Vector'(ranvec);
      put_line("a random monomial vector :"); put(v);
    end;
  end QuadDobl_Main;

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
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_monvec;
