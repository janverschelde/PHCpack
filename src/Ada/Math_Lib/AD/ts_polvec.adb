with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with Standard_Monomial_Vectors;
with Standard_Polynomial_Vectors;
with Standard_Polynomial_Vectors_io;    use Standard_Polynomial_Vectors_io;
with DoblDobl_Monomial_Vectors;
with DoblDobl_Polynomial_Vectors;
with DoblDobl_Polynomial_Vectors_io;    use DoblDobl_Polynomial_Vectors_io;
with QuadDobl_Monomial_Vectors;
with QuadDobl_Polynomial_Vectors;
with QuadDobl_Polynomial_Vectors_io;    use QuadDobl_Polynomial_Vectors_io;
with Random_Polynomial_Vectors;         use Random_Polynomial_Vectors;

procedure ts_polvec is

-- DESCRIPTION :
--   Tests the operations on polynomial vectors.

  procedure Standard_Eval
              ( p : in Standard_Polynomial_Vectors.Link_to_Polynomial_Vector )
  is

  -- DESCRIPTION :
  --   Evaluates p at a random vector, in double precision.

    dim : constant integer32 := integer32(p(p'first).dim);
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    y : Standard_Complex_Vectors.Vector(p'range);

  begin
    y := Standard_Polynomial_Vectors.Eval(p,x);
    put_line("y : "); put_line(y); new_line;
  end Standard_Eval;

  procedure DoblDobl_Eval
              ( p : in DoblDobl_Polynomial_Vectors.Link_to_Polynomial_Vector )
  is

  -- DESCRIPTION :
  --   Evaluates p at a random vector, in double double precision.

    dim : constant integer32 := integer32(p(p'first).dim);
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y : DoblDobl_Complex_Vectors.Vector(p'range);

  begin
    y := DoblDobl_Polynomial_Vectors.Eval(p,x);
    put_line("y : "); put_line(y); new_line;
  end DoblDobl_Eval;

  procedure QuadDobl_Eval
              ( p : in QuadDobl_Polynomial_Vectors.Link_to_Polynomial_Vector )
  is

  -- DESCRIPTION :
  --   Evaluates m at a random vector, in quad double precision.

    dim : constant integer32 := integer32(p(p'first).dim);
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : QuadDobl_Complex_Vectors.Vector(p'range);

  begin
    y := QuadDobl_Polynomial_Vectors.Eval(p,x);
    put_line("y : "); put_line(y); new_line;
  end QuadDobl_Eval;

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
      Standard_Eval(p);
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
      DoblDobl_Eval(p);
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
      QuadDobl_Eval(p);
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
