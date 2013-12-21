with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with Standard_Complex_Row_Reduction;
with DoblDobl_Complex_Row_Reduction;
with QuadDobl_Complex_Row_Reduction;

procedure ts_rowred is

-- DESCRIPTION :
--   Interactive development of row reduction on matrices of
--   standard complex numbers.  Pivoting is done explicitly.

  procedure Standard_Read_Row
              ( A : in out Standard_Complex_Matrices.Matrix;
                i : in integer32;
                piv : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Prompts the user for m complex numbers to fill up row i of A.
  --   The j-th entry is placed according to its position in piv.

  -- REQUIRED : 1 <= i <= A'last(1).

    m : constant integer32 := A'last(2);

  begin
    put("Reading "); put(m,1); put(" complex numbers for row ");
    put(i,1); put_line(" ...");
    for j in 1..m loop
      put("A("); put(i,1); put(","); put(j,1); put(") : ");
      get(A(i,piv(j)));
    end loop;
  end Standard_Read_Row;

  procedure DoblDobl_Read_Row
              ( A : in out DoblDobl_Complex_Matrices.Matrix;
                i : in integer32;
                piv : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Prompts the user for m complex numbers to fill up row i of A.
  --   The j-th entry is placed according to its position in piv.

  -- REQUIRED : 1 <= i <= A'last(1).

    m : constant integer32 := A'last(2);

  begin
    put("Reading "); put(m,1); put(" complex numbers for row ");
    put(i,1); put_line(" ...");
    for j in 1..m loop
      put("A("); put(i,1); put(","); put(j,1); put(") : ");
      get(A(i,piv(j)));
    end loop;
  end DoblDobl_Read_Row;

  procedure QuadDobl_Read_Row
              ( A : in out QuadDobl_Complex_Matrices.Matrix;
                i : in integer32;
                piv : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Prompts the user for m complex numbers to fill up row i of A.
  --   The j-th entry is placed according to its position in piv.

  -- REQUIRED : 1 <= i <= A'last(1).

    m : constant integer32 := A'last(2);

  begin
    put("Reading "); put(m,1); put(" complex numbers for row ");
    put(i,1); put_line(" ...");
    for j in 1..m loop
      put("A("); put(i,1); put(","); put(j,1); put(") : ");
      get(A(i,piv(j)));
    end loop;
  end QuadDobl_Read_Row;

  procedure Standard_Interactive_Row_Reduction ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Asks the user to give n row of m complex numbers
  --   and performs row reduction on each row.

    use Standard_Complex_Row_Reduction;

    A : Standard_Complex_Matrices.Matrix(1..n,1..m);
    piv : Standard_Integer_Vectors.Vector(1..n) := Start_Pivots(n);
    tol : constant double_float := 1.0E-8;
    singular : boolean;
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    for i in 1..n loop
      Standard_Read_Row(A,i,piv);
      if ans = 'y'
       then Reduce_Row(Standard_Output,A,i,piv,tol,singular);
       else Reduce_Row(A,i,piv,tol,singular);
      end if;
      exit when singular;
    end loop;
    put_line("The matrix A after reduction :"); put(A,3);
    if singular
     then put_line("is singular.");
     else put_line("is not singular.");
    end if;
  end Standard_Interactive_Row_Reduction;

  procedure DoblDobl_Interactive_Row_Reduction ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Asks the user to give n row of m complex numbers
  --   and performs row reduction on each row.

    use DoblDobl_Complex_Row_Reduction;

    A : DoblDobl_Complex_Matrices.Matrix(1..n,1..m);
    piv : Standard_Integer_Vectors.Vector(1..n) := Start_Pivots(n);
    tol : constant double_float := 1.0E-8;
    singular : boolean;
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    for i in 1..n loop
      DoblDobl_Read_Row(A,i,piv);
      if ans = 'y'
       then Reduce_Row(Standard_Output,A,i,piv,tol,singular);
       else Reduce_Row(A,i,piv,tol,singular);
      end if;
      exit when singular;
    end loop;
    put_line("The matrix A after reduction :"); put(A,3);
    if singular
     then put_line("is singular.");
     else put_line("is not singular.");
    end if;
  end DoblDobl_Interactive_Row_Reduction;

  procedure QuadDobl_Interactive_Row_Reduction ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Asks the user to give n row of m complex numbers
  --   and performs row reduction on each row.

    use QuadDobl_Complex_Row_Reduction;

    A : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);
    piv : Standard_Integer_Vectors.Vector(1..n) := Start_Pivots(n);
    tol : constant double_float := 1.0E-8;
    singular : boolean;
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    for i in 1..n loop
      QuadDobl_Read_Row(A,i,piv);
      if ans = 'y'
       then Reduce_Row(Standard_Output,A,i,piv,tol,singular);
       else Reduce_Row(A,i,piv,tol,singular);
      end if;
      exit when singular;
    end loop;
    put_line("The matrix A after reduction :"); put(A,3);
    if singular
     then put_line("is singular.");
     else put_line("is not singular.");
    end if;
  end QuadDobl_Interactive_Row_Reduction;

  procedure Main is

    n,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Incremental row reduction of a standard complex matrix ...");
    new_line;
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    put_line("MENU to select the precision : ");
    put_line("  0. standard double floating point arithmetic;");
    put_line("  1. double double floating point arithmetic;");
    put_line("  2. quad double floating point arithmetic.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is 
      when '0' => Standard_Interactive_Row_Reduction(n,m);
      when '1' => DoblDobl_Interactive_Row_Reduction(n,m);
      when '2' => QuadDobl_Interactive_Row_Reduction(n,m);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_rowred;
