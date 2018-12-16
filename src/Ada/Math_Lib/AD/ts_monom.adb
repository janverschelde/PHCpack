with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Monomials;
with DoblDobl_Complex_Monomials;
with QuadDobl_Complex_Monomials;

procedure ts_monom is

-- DESCRIPTION :
--   Tests the operations on monomials.

  function Random_Exponents
             ( dim : integer32; expmax : natural32 )  
             return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of random exponents in the range 0..expmax,
  --   of dimension dim.

    res : Standard_Natural_Vectors.Vector(1..dim);
    rnd : integer32;

  begin
    for i in 1..dim loop
      rnd := Standard_Random_Numbers.Random(0,integer32(expmax));
      res(i) := natural32(rnd);
    end loop;
    return res;
  end Random_Exponents;

  function Standard_Random_Monomial
             ( dim : integer32; expmax : natural32 )  
             return Standard_Complex_Monomials.Monomial is

  -- DESCRIPTION :
  --   Generates a random complex coefficients and an exponent vector
  --   of dimension dim, with entries in the range 0..expmax.

    res : Standard_Complex_Monomials.Monomial;
    exp : constant Standard_Natural_Vectors.Vector
        := Random_Exponents(dim,expmax);
    cff : constant Standard_Complex_Numbers.Complex_Number
        := Standard_Random_Numbers.Random1;

  begin
    put("Coefficient : "); put(cff); new_line;
    put("Exponents : "); put(exp); new_line;
    res := Standard_Complex_Monomials.Create(cff,exp);
    return res;
  end Standard_Random_Monomial;

  function DoblDobl_Random_Monomial
             ( dim : integer32; expmax : natural32 )  
             return DoblDobl_Complex_Monomials.Monomial is

  -- DESCRIPTION :
  --   Generates a random complex coefficients and an exponent vector
  --   of dimension dim, with entries in the range 0..expmax.

    res : DoblDobl_Complex_Monomials.Monomial;
    exp : constant Standard_Natural_Vectors.Vector
        := Random_Exponents(dim,expmax);
    cff : constant DoblDobl_Complex_Numbers.Complex_Number
        := DoblDobl_Random_Numbers.Random1;

  begin
    put("Coefficient : "); put(cff); new_line;
    put("Exponents : "); put(exp); new_line;
    res := DoblDobl_Complex_Monomials.Create(cff,exp);
    return res;
  end DoblDobl_Random_Monomial;

  function QuadDobl_Random_Monomial
             ( dim : integer32; expmax : natural32 )  
             return QuadDobl_Complex_Monomials.Monomial is

  -- DESCRIPTION :
  --   Generates a random complex coefficients and an exponent vector
  --   of dimension dim, with entries in the range 0..expmax.

    res : QuadDobl_Complex_Monomials.Monomial;
    exp : constant Standard_Natural_Vectors.Vector
        := Random_Exponents(dim,expmax);
    cff : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Random_Numbers.Random1;

  begin
    put("Coefficient : "); put(cff); new_line;
    put("Exponents : "); put(exp); new_line;
    res := QuadDobl_Complex_Monomials.Create(cff,exp);
    return res;
  end QuadDobl_Random_Monomial;

  procedure put ( m : in Standard_Complex_Monomials.Monomial ) is

  -- DESCRIPTION :
  --   Writes the contents of the monomial.

  begin
    put("coefficient : "); put(m.cff); new_line;
    put("dimension : "); put(m.dim,1); new_line;
    put("number of variables : "); put(m.nvr,1); new_line;
    put("number of variables with exponent > 1 : "); put(m.n_base,1); new_line;
    if m.nvr > 0 then
      put("positions : "); put(m.pos); new_line;
      put("exponents : "); put(m.exp); new_line;
      if m.n_base > 0 then
        put_line("for variables with exponents > 1 :");
        put("positions : "); put(m.pos_base); new_line;
        put("exponents : "); put(m.exp_base); new_line;
        put("exponents minus 2 : "); put(m.exp_tbl_base); new_line;
      end if;
    end if;
  end put;

  procedure put ( m : in DoblDobl_Complex_Monomials.Monomial ) is

  -- DESCRIPTION :
  --   Writes the contents of the monomial.

  begin
    put("coefficient : "); put(m.cff); new_line;
    put("dimension : "); put(m.dim,1); new_line;
    put("number of variables : "); put(m.nvr,1); new_line;
    put("number of variables with exponent > 1 : "); put(m.n_base,1); new_line;
    if m.nvr > 0 then
      put("positions : "); put(m.pos); new_line;
      put("exponents : "); put(m.exp); new_line;
      if m.n_base > 0 then
        put_line("for variables with exponents > 1 :");
        put("positions : "); put(m.pos_base); new_line;
        put("exponents : "); put(m.exp_base); new_line;
        put("exponents minus 2 : "); put(m.exp_tbl_base); new_line;
      end if;
    end if;
  end put;

  procedure put ( m : in QuadDobl_Complex_Monomials.Monomial ) is

  -- DESCRIPTION :
  --   Writes the contents of the monomial.

  begin
    put("coefficient : "); put(m.cff); new_line;
    put("dimension : "); put(m.dim,1); new_line;
    put("number of variables : "); put(m.nvr,1); new_line;
    put("number of variables with exponent > 1 : "); put(m.n_base,1); new_line;
    if m.nvr > 0 then
      put("positions : "); put(m.pos); new_line;
      put("exponents : "); put(m.exp); new_line;
      if m.n_base > 0 then
        put_line("for variables with exponents > 1 :");
        put("positions : "); put(m.pos_base); new_line;
        put("exponents : "); put(m.exp_base); new_line;
        put("exponents minus 2 : "); put(m.exp_tbl_base); new_line;
      end if;
    end if;
  end put;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial with coefficient in standard double precision.

    dim : integer32 := 0;
    expmax : natural32 := 0;
    m : Standard_Complex_Monomials.Monomial;

  begin
    put_line("Testing monomial operations in double precision ...");
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    m := Standard_Random_Monomial(dim,expmax);
    put_line("A random monomial :");
    put(m);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial with coefficient in double double precision.

    dim : integer32 := 0;
    expmax : natural32 := 0;
    m : DoblDobl_Complex_Monomials.Monomial;

  begin
    put_line("Testing monomial operations in double double precision ...");
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    m := DoblDobl_Random_Monomial(dim,expmax);
    put_line("A random monomial :"); put(m);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial with coefficient in quad double precision.

    dim : integer32 := 0;
    expmax : natural32 := 0;
    m : QuadDobl_Complex_Monomials.Monomial;

  begin
    put_line("Testing monomial operations in quad double precision ...");
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    m := QuadDobl_Random_Monomial(dim,expmax);
    put_line("A random monomial :"); put(m);
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
end ts_monom;
