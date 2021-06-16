with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Systems_Pool;
with DoblDobl_Systems_Pool;
with QuadDobl_Systems_Pool;

procedure ts_syspool is

-- DESCRIPTION :
--   Tests some basic operations on the systems pools.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Runs a basic test to store systems of polynomials with
  --   complex coefficients in standard double precision.

    use Standard_Complex_Poly_Systems;

    n : integer32 := 0;
    p : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the operations in the systems pool...");
    new_line;
    put("-> The size of the systems pool : ");
    put(Standard_Systems_Pool.Size,1); new_line;
    put("Give the number of systems : "); get(n);
    Standard_Systems_Pool.Initialize(n);
    put("-> The size of the systems pool : ");
    put(Standard_Systems_Pool.Size,1); new_line;
    new_line;
    put_line("Reading a polynomial system..."); get(p);
    for k in 1..n loop
      Standard_Systems_Pool.Create(k,p.all);
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Runs a basic test to store systems of polynomials with
  --   complex coefficients in double double precision.

    use DoblDobl_Complex_Poly_Systems;

    n : integer32 := 0;
    p : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the operations in the systems pool...");
    new_line;
    put("-> The size of the systems pool : ");
    put(DoblDobl_Systems_Pool.Size,1); new_line;
    put("Give the number of systems : "); get(n);
    DoblDobl_Systems_Pool.Initialize(n);
    put("-> The size of the systems pool : ");
    put(DoblDobl_Systems_Pool.Size,1); new_line;
    new_line;
    put_line("Reading a polynomial system..."); get(p);
    for k in 1..n loop
      DoblDobl_Systems_Pool.Create(k,p.all);
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Runs a basic test to store systems of polynomials with
  --   complex coefficients in quad double precision.

    use QuadDobl_Complex_Poly_Systems;

    n : integer32 := 0;
    p : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the operations in the systems pool...");
    new_line;
    put("-> The size of the systems pool : ");
    put(QuadDobl_Systems_Pool.Size,1); new_line;
    put("Give the number of systems : "); get(n);
    QuadDobl_Systems_Pool.Initialize(n);
    put("-> The size of the systems pool : ");
    put(QuadDobl_Systems_Pool.Size,1); new_line;
    new_line;
    put_line("Reading a polynomial system..."); get(p);
    for k in 1..n loop
      QuadDobl_Systems_Pool.Create(k,p.all);
    end loop;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the corresponding test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_syspool;
