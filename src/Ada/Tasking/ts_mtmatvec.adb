with text_io;                            use text_io;
with Multitasking;                       use Multitasking;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Random_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Random_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Random_Matrices;
with Multitasking_Matrix_x_Vector;       use Multitasking_Matrix_x_Vector;

procedure ts_mtmatvec is

-- DESCRIPTION :
--   Interactive development of multitasking and matrix-vector operations.

  procedure Run_Standard_Multiply
              ( A : in Standard_Complex_Matrices.Matrix;
                v : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Matrices;

    y : constant Standard_Complex_Vectors.Vector := A*v;
    nt,m : integer32 := 0;
    ans : character;
    z : Standard_Complex_Vectors.Vector(A'range(1));

  begin
    put("Give number of tasks : "); get(nt);
    put("Give number of times : "); get(m);
    put("Do you want output ? (y/n) "); get(ans);
    if ans = 'y' then
      for i in 1..m loop
        z := Reporting_Multiply(nt,A,v);
      end loop;
    else
      for i in 1..m loop
        z := Silent_Multiply(nt,A,v);
      end loop;
    end if;
    put_line("y : "); put_line(y);
    put_line("z : "); put_line(z);
  end Run_Standard_Multiply;

  procedure Run_DoblDobl_Multiply
              ( A : in DoblDobl_Complex_Matrices.Matrix;
                v : in DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Matrices;

    y : constant DoblDobl_Complex_Vectors.Vector := A*v;
    nt,m : integer32 := 0;
    ans : character;
    z : DoblDobl_Complex_Vectors.Vector(A'range(1));

  begin
    put("Give number of tasks : "); get(nt);
    put("Give number of times : "); get(m);
    put("Do you want output ? (y/n) "); get(ans);
    if ans = 'y' then
      for i in 1..m loop
        z := Reporting_Multiply(nt,A,v);
      end loop;
    else
      for i in 1..m loop
        z := Silent_Multiply(nt,A,v);
      end loop;
    end if;
    put_line("y : "); put_line(y);
    put_line("z : "); put_line(z);
  end Run_DoblDobl_Multiply;

  procedure Run_QuadDobl_Multiply
              ( A : in QuadDobl_Complex_Matrices.Matrix;
                v : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Matrices;

    y : constant QuadDobl_Complex_Vectors.Vector := A*v;
    nt,m : integer32 := 0;
    ans : character;
    z : QuadDobl_Complex_Vectors.Vector(A'range(1));

  begin
    put("Give number of tasks : "); get(nt);
    put("Give number of times : "); get(m);
    put("Do you want output ? (y/n) "); get(ans);
    if ans = 'y' then
      for i in 1..m loop
        z := Reporting_Multiply(nt,A,v);
      end loop;
    else
      for i in 1..m loop
        z := Silent_Multiply(nt,A,v);
      end loop;
    end if;
    put_line("y : "); put_line(y);
    put_line("z : "); put_line(z);
  end Run_QuadDobl_Multiply;

  procedure Standard_Test ( r,c : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a test using standard complex arithmetic
  --   on a r-by-c matrix of complex numbers.

    A : constant Standard_Complex_Matrices.Matrix
      := Standard_Random_Matrices.Random_Matrix(r,c); 
    v : constant Standard_Complex_Vectors.Vector
      := Standard_Random_Vectors.Random_Vector(1,integer32(c));

  begin
    Run_Standard_Multiply(A,v);
  end Standard_Test;

  procedure DoblDobl_Test ( r,c : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a test using complex double double arithmetic,
  --   on a r-by-c matrix of complex double doubles.

    A : constant DoblDobl_Complex_Matrices.Matrix
      := DoblDobl_Random_Matrices.Random_Matrix(r,c); 
    v : constant DoblDobl_Complex_Vectors.Vector
      := DoblDobl_Random_Vectors.Random_Vector(1,integer32(c));

  begin
    Run_DoblDobl_Multiply(A,v);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( r,c : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a test using complex quad double arithmetic,
  --   on a r-by-c matrix of complex quad doubles.

    A : constant QuadDobl_Complex_Matrices.Matrix
      := QuadDobl_Random_Matrices.Random_Matrix(r,c); 
    v : constant QuadDobl_Complex_Vectors.Vector
      := QuadDobl_Random_Vectors.Random_Vector(1,integer32(c));

  begin
    Run_QuadDobl_Multiply(A,v);
  end QuadDobl_Test;

  procedure Main is

    ans : character;
    r,c : natural32 := 0;

  begin
    new_line;
    put_line("MENU to test multitasked matrix-vector multiplications :");
    put_line("  1. use standard complex aritmetic;");
    put_line("  2. use complex double double aritmetic;");
    put_line("  3. use complex quad double aritmetic.");
    put("Type 1, 2, or 3 to choose : "); get(ans);
    new_line;
    put("Give number of rows : "); get(r);
    put("Give number of columns : "); get(c);
    if ans = '1' then
      Standard_Test(r,c);
    elsif ans = '2' then
      DoblDobl_Test(r,c);
    else
      QuadDobl_Test(r,c);
    end if;
  end Main;

begin
  Main;
end ts_mtmatvec;
