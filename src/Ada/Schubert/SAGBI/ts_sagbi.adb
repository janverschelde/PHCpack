with text_io;                            use text_io;
--with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Matrix_Indeterminates;              use Matrix_Indeterminates;
with SAGBI_Homotopies;                   use SAGBI_Homotopies;
--with Drivers_for_SAGBI_Homotopies;       use Drivers_for_SAGBI_Homotopies;

procedure ts_sagbi is

  procedure Write ( p : in Poly ) is

    procedure Write_Term ( t : in Term; cont : out boolean ) is
    begin
      put("t.cf : "); put(t.cf); new_line;
      put("t.dg : ");
      if t.dg /= null
       then put(t.dg.all); new_line;
      end if;
      cont := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
  end Write;

  procedure Write_Coeff ( p : in Poly ) is

    cff : constant Vector := Coeff(p);

  begin
    put_line("The coefficient vector : ");
    put_line(cff);
  end Write_Coeff;

  procedure Test_SAGBI_Homotopy ( n,d : in natural32 ) is

    p,l,p1 : Poly;
    mat : Matrix(1..integer32(n),1..integer32(d));
    ans : character;
    sb : Symbol;

  begin
    sb := (sb'range => ' ');
    sb(sb'first) := 't';
    loop
      Matrix_Indeterminates.Initialize_Symbols(n,d);
      p := Lifted_Localized_Laplace_Expansion(n,d);
      Symbol_Table.Replace((n-d)*d+1,sb);
      new_line;
      put_line("The generic polynomial : "); put(p); new_line;
      Write_Coeff(p);
      p1 := Lifted_Localized_Laplace_Expansion(n,d);
      put_line("The new generic polynomial : "); put(p1); new_line;
      Write(p1);
      put("Give a floating-point "); put(n,1); put("x"); put(d,1);
      put_line("-matrix : "); get(mat);
      l := Intersection_Condition(mat,p);
      put_line("The specific polynomial : "); put(l); new_line;
      Write_Coeff(l);
      declare
        cff : constant Vector := Intersection_Coefficients(mat,Coeff(p));
      begin
        put_line("The computed coefficient vector : ");
        put_line(cff);
      end;
      Matrix_Indeterminates.Clear_Symbols;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_SAGBI_Homotopy;

  procedure Main is

    n,d : natural32 := 0;
   -- file : file_type;

  begin
    new_line;
    put_line("SAGBI Homotopies to intersect planes in projective space.");
   -- new_line;
   -- put_line("Reading the name of the output file.");
   -- Read_Name_and_Create_File(file);
    new_line;
    put("Give number of elements to choose from : "); get(n);
    put("Give the number of entries in bracket : "); get(d);
   -- put(file,"SAGBI Homotopies for n = "); put(file,n,1);
   -- put(file," and d = "); put(file,d,1); new_line(file);
    Test_SAGBI_Homotopy(n,d);
   -- Driver_for_SAGBI_Homotopies(file,n,d);
  end Main;

begin
  Main;
end ts_sagbi;
