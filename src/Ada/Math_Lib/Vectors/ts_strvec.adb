with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Complex_Vector_Strings;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Vector_Strings;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Strings;

procedure ts_strvec is

-- DESCRIPTION :
--   Writing vectors to strings and parsing strings from vectors.

  procedure Standard_Random_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Test on a random vector of dimension n.

    v : Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    s : constant string
      := Standard_Complex_Vector_Strings.Write(v);

  begin
    put_line("The random vector : "); put_line(v);
    put_line("The vector written as string :"); put_line(s);
    declare
      w : constant Standard_Complex_Vectors.Vector
        := Standard_Complex_Vector_Strings.Parse(s);
    begin
      put_line("After parsing from string, the vector :"); put_line(w);
    end;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Test on a random vector of dimension n.

    v : DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    s : constant string
      := DoblDobl_Complex_Vector_Strings.Write(v);

  begin
    put_line("The random vector : "); put_line(v);
    put_line("The vector written as string :"); put_line(s);
    declare
      w : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Complex_Vector_Strings.Parse(s);
    begin
      put_line("After parsing from string, the vector :"); put_line(w);
    end;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Test on a random vector of dimension n.

    v : QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    s : constant string
      := QuadDobl_Complex_Vector_Strings.Write(v);

  begin
    put_line("The random vector : "); put_line(v);
    put_line("The vector written as string :"); put_line(s);
    declare
      w : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Complex_Vector_Strings.Parse(s);
    begin
      put_line("After parsing from string, the vector :"); put_line(w);
    end;
  end QuadDobl_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension n and then launches a test.

    n : integer32 := 0;
    ans : character;

  begin
    put("Give the dimension of the vector : "); get(n);
    new_line;
    put_line("MENU for the precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Random_Test(n);
      when '1' => DoblDobl_Random_Test(n);
      when '2' => QuadDobl_Random_Test(n);
      when others => null;
    end case;
  end Main;
 
begin
  Main;
end ts_strvec;
