with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with String_Splitters;                   use String_Splitters;
with Symbol_Table;
with Standard_Complex_Laur_Strings;      use Standard_Complex_Laur_Strings;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;

procedure ts_nbrvar is

-- DESCRIPTION :
--   Interactive development of the determination of
--   the number of variables in a polynomial, given as a string.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a string and the maximum number of
  --   of variables, used to initialize the symbol table.

    strpol : constant String := Read_String;
    maxvar : natural32 := 0;
    p : Poly;
    size : natural32;
  
  begin
    new_line;
    put_line("Your string : " & strpol);
    new_line;
    put("Give the maximum number of variables : ");
    get(maxvar);
    Symbol_Table.Init(maxvar);
    p := Parse(maxvar,strpol);
    new_line;
    put("Your polynomial : "); put(p); new_line;
    size := Size_of_Support(p);
    put("The number of variables in p : "); put(size,1); new_line;
  end Main;

begin
  new_line;
  put_line("Reading a string with a polynomial ...");
  Main;
end ts_nbrvar;
