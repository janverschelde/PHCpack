with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Tableau_Formats;

procedure convert_formats is

-- DESCRIPTION :
--   This program allows to convert between formats of polynomial systems.
--
-- SYMBOLIC FORMAT :
--
--   2
--    x*y**2 - x**2 + 3;
--    x + 2;
--
-- TABLEAU FORMAT :
--
--   2
--    x y
--   3
--    1 2
--    2 0
--    0 0
--  2 
--    1 0
--    0 0
--    1.0 0.0
--   -1.0 0.0
--    3.0 0.0
--    1.0 0.0
--    2.0 0.0

  procedure Tableau_to_Symbolic ( infile,outfile : in file_type;
                                  flt : in boolean ) is

    lp : Link_to_Poly_Sys;

  begin
    Tableau_Formats.get(infile,flt,lp);
    put(outfile,lp'last,lp.all);
  end Tableau_to_Symbolic;

  procedure Symbolic_to_Tableau ( infile,outfile : in file_type;
                                  flt : in boolean ) is

    lp : Link_to_Poly_Sys;

  begin
    get(infile,lp);
    Tableau_Formats.put(outfile,flt,lp.all);
  end Symbolic_to_Tableau;

  procedure Main is

    infile,outfile : file_type;
    ans : character;
    flt,sym2tab : boolean;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Conversion between symbolic and tableau formats.");
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. Convert from symbolic to tableau format");
    put_line("  2. Convert from tableau to symbolic format");
    put("Type 1 or 2 to select : "); Ask_Alternative(ans,"12");
    sym2tab := (ans = '1');
    new_line;
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(infile);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Are the coefficients real or complex ? (r/c) ");
    Ask_Alternative(ans,"rc");
    flt := (ans = 'r');
    if sym2tab
     then Symbolic_to_Tableau(infile,outfile,flt);
     else Tableau_to_Symbolic(infile,outfile,flt);
    end if;
  end Main;

begin
  Main;
end convert_formats;
