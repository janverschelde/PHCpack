with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;    use Standard_Integer64_Matrices_io;
with Multprec_Integer_Matrices;
with Multprec_Integer_Matrices_io;      use Multprec_Integer_Matrices_io;
with Convex_Hull_Methods;               use Convex_Hull_Methods;
with Point_Lists_and_Strings;

procedure ts_ptlstr is

-- DESCRIPTION :
--   Tests the conversions from integer matrices to Python lists of tuples.

  procedure Standard_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Tests with standard 64-bit integers.

    A : Standard_Integer64_Matrices.Matrix(1..n,1..m);
   
  begin
    A := Random_Data(n,m);
    put_line("The point configuration : "); put(A);
    put_line("The string representation : ");
    declare
      s : constant string := Point_Lists_and_Strings.convert(A);
    begin
      put_line(s);
    end;
  end Standard_Test;

  procedure Multprec_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Test with arbitrary multiprecision integers.

    A : Multprec_Integer_Matrices.Matrix(1..n,1..m);
   
  begin
    A := Random_Data(n,m);
    put_line("The point configuration : "); put(A);
    put_line("The string representation : ");
    declare
      s : constant string := Point_Lists_and_Strings.convert(A);
    begin
      put_line(s);
    end;
  end Multprec_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for dimensions and whether multiprecision or not.

    n,m : integer32 := 0;
    ans : character;

  begin
    put("Give the dimension : "); get(n);
    put("Give the number of points : "); get(m);
    put("Multiprecision arithmetic ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multprec_Test(n,m);
     else Standard_Test(n,m);
    end if;
  end Main;

begin
  Main;
end ts_ptlstr;
