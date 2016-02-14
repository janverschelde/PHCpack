with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with String_Splitters;                   use String_Splitters;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Parse_Dimensions;

procedure ts_nbrvar is

-- DESCRIPTION :
--   Interactive development of the determination of
--   the number of variables in a polynomial, given as a string.

  procedure Test_Polynomial is

  -- DESCRIPTION :
  --   Prompts the user for a string and the maximum number of
  --   of variables, used to initialize the symbol table.

    strpol : constant string := Read_String;
    maxvar : natural32 := 0;
    size : natural32;
    p : Poly;
  
  begin
    new_line;
    put_line("Your string : " & strpol);
    new_line;
    put("Give the maximum number of variables : ");
    get(maxvar);
    size := Parse_Dimensions.Dim(maxvar,strpol);
    put("The number of symbols : "); put(size,1); new_line;
    p := Parse_Dimensions.Get;
    put("The parsed polynomial : "); put(p); new_line;
    Parse_Dimensions.Clear;
  end Test_Polynomial;

  procedure Test_System is

  -- DESCRIPTION :
  --   Prompts the user for the number of polynomials
  --   and maximum number of variables.
  --   Reads as many strings as the number of polynomials.

    nq : integer32 := 0;
    ls : Link_to_Array_of_Strings;
    maxvar,size : natural32 := 0;
    s : Link_to_Laur_Sys;

  begin
    put("Give the number of polynomials : "); get(nq);
    skip_line;
    ls := new Array_of_Strings(1..integer(nq));
    for k in ls'range loop
      put("  reading polynomial "); put(integer32(k),1); 
      put_line(" ...");
      declare
        lsk : constant string := Read_String;
      begin
        ls(k) := new string'(lsk);
      end;
    end loop;
    new_line;
    put_line("Your system :");
    for k in ls'range loop
      put_line(ls(k).all);
    end loop;
    new_line;
    put("Give the maximum number of variables : ");
    get(maxvar);
    size := Parse_Dimensions.Dim(maxvar,ls.all);
    put("The number of symbols : "); put(size,1); new_line;
    s := Parse_Dimensions.Get;
    put_line("The parsed system :"); put(s.all); new_line;
    Parse_Dimensions.Clear;
  end Test_System;

  procedure Main is
  begin
    new_line;
    put_line("Reading a string with a polynomial ...");
    Test_Polynomial;
    new_line;
    put_line("Testing the parsing of a system ...");
    new_line;
    Test_System;
  end Main;

begin
  Main;
end ts_nbrvar;
