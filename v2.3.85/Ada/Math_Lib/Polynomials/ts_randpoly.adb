with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Random_Polynomials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Random_Polynomials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Random_Polynomials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Multprec_Random_Polynomials;

procedure ts_randpoly is

-- DESCRIPTION :
--   Lets the user generate random sparse and dense polynomials.

  procedure Save_System_to_File
              ( s : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Save_System_to_File
              ( s : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Save_System_to_File
              ( s : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Save_System_to_File
              ( s : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Standard_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32 ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients and show it.
  --   The number of polynomials equals e.

    s : Standard_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : Standard_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := Standard_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := Standard_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
  end Standard_Generate_and_Show;

  procedure DoblDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32 ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients and show it.

    s : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : DoblDobl_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := DoblDobl_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := DoblDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
  end DoblDobl_Generate_and_Show;

  procedure QuadDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32 ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients and show it.

    s : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : QuadDobl_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := QuadDobl_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := QuadDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
  end QuadDobl_Generate_and_Show;

  procedure Multprec_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32 ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients and show it.

    s : Multprec_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : Multprec_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := Multprec_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := Multprec_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
  end Multprec_Generate_and_Show;

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
      if ans = '0' then
        Standard_Generate_and_Show(n,d,m,c,e);
      elsif ans = '1' then
        DoblDobl_Generate_and_Show(n,d,m,c,e);
      elsif ans = '2' then
        QuadDobl_Generate_and_Show(n,d,m,c,e);
      else
        Multprec_Generate_and_Show(n,d,m,c,e);
      end if;
   -- end loop;
  end Main;

begin
  Main;
end ts_randpoly;
