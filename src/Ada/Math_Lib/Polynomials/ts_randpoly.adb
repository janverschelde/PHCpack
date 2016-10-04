with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Random_Polynomial_Systems;         use Random_Polynomial_Systems;

procedure ts_randpoly is

-- DESCRIPTION :
--   Lets the user generate random sparse and dense polynomials.

  procedure Main is

    n,d,m,c : natural32 := 0;
    e : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Generation of random dense and sparse polynomial systems.");
   -- loop
      new_line;
      put("Give the number of variables : "); get(n);
   --   put("Give the number of variables (0 to exit): "); get(n);
   --   exit when (n = 0);
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
      put_line("MENU for type of coefficients : ");
      put_line("  0. standard double floats;");
      put_line("  1. double double coefficients;");
      put_line("  2. quad double coefficients.");
      put_line("  3. multiprecision coefficients.");
      put("Type 0, 1, 2, or 3 to make a choice : ");
      Ask_Alternative(ans,"0123");
      case ans is
        when '0' =>
          declare
            lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
          begin
            Standard_Generate_and_Show(n,d,m,c,e,lp);
          end;
        when '1' => 
          declare
            lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
          begin
            DoblDobl_Generate_and_Show(n,d,m,c,e,lp);
          end;
        when '2' =>
          declare
            lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
          begin
            QuadDobl_Generate_and_Show(n,d,m,c,e,lp);
          end;
        when others =>
          declare
            lp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
          begin
            Multprec_Generate_and_Show(n,d,m,c,e,lp);
          end;
      end case;
   -- end loop;
  end Main;

begin
  Main;
end ts_randpoly;
