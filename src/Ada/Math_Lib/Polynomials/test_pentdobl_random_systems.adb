with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Symbol_Table;
with PentDobl_Complex_Poly_Systems;
with Random_Polynomial_Systems;         use Random_Polynomial_Systems;

package body Test_PentDobl_Random_Systems is

  procedure Main is

    n,d,m,c : natural32 := 0;
    e : integer32 := 0;
    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Generation of random dense and sparse polynomial systems.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    put("Give the number of polynomials : "); get(e);
    new_line;
    put_line("MENU to generate coefficients : ");
    put_line("  0 = on complex unit circle;");
    put_line("  1 = equal to the constant one;");
    put_line("  2 = random real numbers between -1 and +1;");
    put("Give natural number to make your choice : "); get(c);
    new_line;
    PentDobl_Generate_and_Show(n,d,m,c,e,lp);
  end Main;

end Test_PentDobl_Random_Systems;
