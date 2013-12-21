with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Communications_with_User;           use Communications_with_User;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;

procedure ts_square is

-- DESCRIPTION :
--   A polynomial system with any number of equations and unknowns
--   can be turned into a square system,
--    1) by adding random slices in case it is underdetermined;
--    2) by adding slack variables in case it is overdetermined. 
-- NOTE : adding zz1 = 0 as last equation is not very healty since
--    the last slice is used for slicing...

  lp : Link_to_Poly_Sys;

  procedure Write_Symbol_Table is
  begin
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(Symbol_Table.Get(i)); put(" ");
    end loop;
    new_line;
  end Write_Symbol_Table;

  function Maximum ( a,b : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximum of a and b.

  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

  procedure Make_Square ( p : in Poly_Sys ) is

    nbequ : constant natural32 := natural32(p'length);
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Poly_Sys(1..integer32(max)) := Square(p);
    ans : character;
    outfile : file_type;

  begin
    put("Number of unknowns  : "); put(nbunk,1); new_line; 
    put_line("The unknowns :"); Write_Symbol_Table;
    put("Number of equations : "); put(nbequ,1); new_line; 
    if nbequ > nbunk
     then Add_Slack_Symbols(nbequ-nbunk);
    end if;
    put_line("The system after squaring : "); put_line(sp);
    put("Do you want to save the square system to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(outfile);
      put_line(outfile,sp);
    end if;
  end Make_Square;

  procedure Main is

  begin
    new_line;
    put_line("Turning a polynomial system into a square system,");
    put_line("by adding random slices or slack variables.");
    new_line;
    get(lp);
    put_line("your polynomial system : "); put(lp.all);
    Make_Square(lp.all);
  end Main;

begin
  Main;
end ts_square;
