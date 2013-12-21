with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Bracket_Expansions;                 use Bracket_Expansions;
with Matrix_Indeterminates;

procedure ts_local is

  procedure Main is

    n,d : natural32 := 0;
    ans : character;

  begin
    put("Give n, the dimension of the space : "); get(n);
    put("Give d, the dimension of the plane : "); get(d);
    declare
      locmap : constant Matrix(1..integer32(n),1..integer32(d))
             := Localization_Map(n,d);
      b : Bracket(1..integer32(d));
      e : Poly;
    begin
      put_line("The localization map : "); put(locmap);
      loop
        Matrix_Indeterminates.Initialize_Symbols(n,d);
        put("Give "); put(d,1); put(" numbers for a bracket : "); get(b);
        e := Expand(locmap,b);
        put("The expansion of minor "); put(b);
        put_line(" following the localization map : ");
        put(e); new_line;
        Matrix_Indeterminates.Reduce_Symbols(locmap);
        Reduce_Variables(locmap,e);
        put_line("The expansion polynomial in reduced format : ");
        put(e); new_line;
        put("Do you wish to give other brackets ? (y/n) "); get(ans);
        exit when (ans /= 'y');
        Matrix_Indeterminates.Clear_Symbols;
      end loop;
    end;
  end Main;

begin
  new_line;
  put_line("Test of setting up a spanning localization.");
  new_line;
  Main;
end ts_local;
