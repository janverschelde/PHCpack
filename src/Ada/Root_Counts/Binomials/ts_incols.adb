with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;      use Standard_Floating_VecVecs_io;
with Standard_Floating_Column_Span;     use Standard_Floating_Column_Span;

procedure ts_incols is

-- DESCRIPTION :
--   Interactive testing of in column span package.

  procedure Read_Vectors
              ( n,m : in integer32;
                v : in out Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Prompts the user to enter m vectors of length n.

  begin
    for i in 1..m loop
      put("Give vector "); put(i,1); put(" : ");
      v(i) := new Standard_Floating_Vectors.Vector(1..n);
      get(v(i).all);
    end loop;
  end Read_Vectors;

  procedure Test_Span ( n,m : in integer32 ) is

    v : Standard_Floating_VecVecs.VecVec(1..m);
    x : Standard_Floating_Vectors.Vector(1..n);

  begin
    Read_Vectors(n,m,v);
    put_line("The vectors : "); put(v);
    put("Give "); put(n,1); put("-vector : "); get(x);
    put_line("Test vector : "); put(x); new_line;
    if In_Span(standard_output,v,x,1.0e-8)
     then put_line("Test vector lies in span.");
     else put_line("Test vector does not lie in span.");
    end if;
  end Test_Span;


  procedure Main is

    n,m : integer32 := 0;

  begin
    put("Give ambient dimension : "); get(n);
    put("Give number of columns : "); get(m);
    Test_Span(n,m);
  end Main;

begin
  Main;
end ts_incols;
