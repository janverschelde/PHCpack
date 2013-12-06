with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_Lists;        use Standard_Complex_Laur_Lists;

procedure ts_laurlist is

-- DESCRIPTION :
--   Basic test on storing standard complex Laurent polynomials in a list.

  procedure Write ( pl : in Poly_List ) is

    tmp : Poly_List := pl;

  begin
    put("Number of polynomials : ");
    put(Length_Of(pl),1); new_line;
    while not Is_Null(tmp) loop
      put(Head_Of(tmp)); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end Write;

  procedure Main is

    lp : Link_to_Laur_Sys;
    pl : Poly_List;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    put_line("your system : "); put(lp.all);
    pl := Create(lp.all);
    Write(pl);
  end Main;

begin
  Main;
end ts_laurlist;
