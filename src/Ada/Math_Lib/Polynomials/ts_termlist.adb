with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Random_Polynomials;
with Standard_Complex_Term_Lists;

procedure ts_termlist is

-- DESCRIPTION :
--   Test on the package list of terms.

  procedure put ( p : in Standard_Complex_Term_Lists.Term_List ) is

  -- DESCRIPTION :
  --   Writes the terms in the list p in tableau format.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Term_Lists;

    tmp : Term_List := p;
    t : Term;

  begin
    while not Is_Null(tmp) loop
      t := Head_Of(tmp);
      put(t.cf);
      put(t.dg.all);
      new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure Main is

    n,d,m : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    t : Standard_Complex_Term_Lists.Term_List;

  begin
    new_line;
    put_line("Generating a random polynomial with complex coefficients ...");
    put("Give the number of variables : "); get(n);
    put("Give the largest degree : "); get(d);
    put("Give the number of terms : "); get(m);
    p := Standard_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
    put("The random polynomial : "); put_line(p);
    t := Standard_Complex_Term_Lists.Create(p);
    put_line("The list of terms :"); put(t);
    Standard_Complex_Term_Lists.Clear(t);
    put_line("The list of terms after clear :"); put(t);
  end Main;

begin
  Main;
end ts_termlist;
