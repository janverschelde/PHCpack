with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Complex_Vector_Strings;

procedure ts_strvec is

-- DESCRIPTION :
--   Writing vectors to strings and parsing strings from vectors.

  procedure Standard_Random_Test( n : in integer32 ) is

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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension n and then launches a test.

    n : integer32 := 0;

  begin
    put("Give the dimension of the vector : "); get(n);
    Standard_Random_Test(n);
  end Main;
 
begin
  Main;
end ts_strvec;
