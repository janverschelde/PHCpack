with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Lists_of_Floating_Vectors;         use Lists_of_Floating_Vectors;
with Lists_of_Floating_Vectors_io;      use Lists_of_Floating_Vectors_io;
with Arrays_of_Floating_Vector_Lists;   use Arrays_of_Floating_Vector_Lists;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Supported_Subsystems;              use Supported_Subsystems;

procedure ts_supsub is

-- DESCRIPTION :
--   Test on the operations in supported subsystems.

  procedure Select_Supported_Subpolynomial ( p : in Poly ) is

    n : constant natural32 := Number_of_Unknowns(p);
    k : natural32 := 0;
    s : List;
    sub_p : Poly;

  begin
    new_line;
    put("Give #selected monomials from polynomial : ");
    get(k);
    put("Give "); put(k,1); put(" vectors of length ");
    put(n,1); put_line(" :"); get(n,k,s);
    put_line("The selected support : "); put(s);
    sub_p := Select_Terms(p,s);
    put_line("The selected polynomial : "); put(sub_p); new_line;
  end Select_Supported_Subpolynomial;

  procedure Select_Supported_Subsystem
              ( p : in Poly_Sys; r : in integer32 ) is

    n : constant natural32 := natural32(p'last);
    k : natural32;
    m : Standard_Integer_Vectors.Vector(1..r);
    s : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    sub_p : Poly_Sys(p'range);

  begin
    put("Give the type of mixture : "); get(m);
    for i in 1..r loop
      new_line;
      put("Give #selected monomials from support "); put(i,1);
      put(" : "); get(k);
      put("Give "); put(k,1); put(" vectors of length ");
      put(n,1); put_line(" :"); get(n,k,s(i));
      put_line("The selected support : "); put(s(i));
    end loop;
    sub_p := Select_Terms(p,m,s);
    put_line("The selected subsystem is "); put(sub_p); 
  end Select_Supported_Subsystem;

  procedure Select_Supported_Subsystem ( p : in Poly_Sys ) is

    n : constant natural32 := natural32(p'last);
    k : natural32;
    s : Arrays_of_Floating_Vector_Lists.Array_of_Lists(p'range);
    sub_p : Poly_Sys(p'range);

  begin
    for i in p'range loop
      new_line;
      put("Give #selected monomials from polynomial "); put(i,1);
      put(" : "); get(k);
      put("Give "); put(k,1); put(" vectors of length ");
      put(n,1); put_line(" :"); get(n,k,s(i));
      put_line("The selected support : "); put(s(i));
    end loop;
    sub_p := Select_Terms(p,s);
    put_line("The selected subsystem is "); put(sub_p); 
  end Select_Supported_Subsystem;

  procedure Main is

    ans : character;
    lp : Link_to_Poly_Sys;
    n : natural32 := 0;
    r : integer32 := 0;
    p : Poly;

  begin
    new_line;
    put_line("Testing the selection of supported subsystems...");
    new_line;
    put_line("Choose one of the following operations : ");
    put_line("  1. Select terms from a given polynomial; or");
    put_line("  2. Select a supported subsystem of a fully-mixed system.");
    put_line("  3. Select a supported subsystem of a semi-mixed system.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    if ans = '1' then
      put("Give the number of unknowns : "); get(n);
      Symbol_Table.Init(n);
      put("Give a polynomial in "); put(n,1); put(" variables,");
      put_line(" terminate with semi-colon:");
      get(p);
      Select_Supported_Subpolynomial(p);
    elsif ans = '2' then
      get(lp);
      Select_Supported_Subsystem(lp.all);
    else
      get(lp);
      new_line;
      put("Give the number of different supports : "); get(r);
      Select_Supported_Subsystem(lp.all,r);
    end if;
  end Main;

begin
  Main;
end ts_supsub;
