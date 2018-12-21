with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Matrix_Norms;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrix_Norms;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrix_Norms;
with Standard_Monomial_Vectors;
with Standard_Polynomial_Vectors;
with Standard_Polynomial_Vectors_io;    use Standard_Polynomial_Vectors_io;
with DoblDobl_Monomial_Vectors;
with DoblDobl_Polynomial_Vectors;
with DoblDobl_Polynomial_Vectors_io;    use DoblDobl_Polynomial_Vectors_io;
with QuadDobl_Monomial_Vectors;
with QuadDobl_Polynomial_Vectors;
with QuadDobl_Polynomial_Vectors_io;    use QuadDobl_Polynomial_Vectors_io;
with Random_Polynomial_Vectors;         use Random_Polynomial_Vectors;

procedure ts_polvec is

-- DESCRIPTION :
--   Tests the operations on polynomial vectors.

  procedure Write ( A,B : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes all entries of the two matrices A and B consecutively,
  --   useful for comparing differences.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2).

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put(") : ");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put(") : ");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( A,B : in DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes all entries of the two matrices A and B consecutively,
  --   useful for comparing differences.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2).

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put(") : ");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put(") : ");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( A,B : in QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes all entries of the two matrices A and B consecutively,
  --   useful for comparing differences.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2).

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put(") : ");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put(") : ");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Standard_Eval ( p : in Standard_Polynomial_Vectors.System ) is

  -- DESCRIPTION :
  --   Evaluates p at a random vector, in double precision.

    x : constant Standard_Complex_Vectors.Vector(1..p.dim)
      := Standard_Random_Vectors.Random_Vector(1,p.dim);
    y,z : Standard_Complex_Vectors.Vector(1..p.nbr);
    dm,ym,zm : Standard_Complex_Matrices.Matrix(y'range,x'range);
    nrm : double_float;
    ans : character;

  begin
    y := Standard_Polynomial_Vectors.Eval(p,x);
    put_line("y : "); put_line(y); new_line;
    Standard_Polynomial_Vectors.Diff(p,x,ym);
    Standard_Polynomial_Vectors.Speel(p,x,z,zm);
    put_line("z : "); put_line(z); new_line;
    Standard_Complex_Vectors.Sub(y,z);
    nrm := Standard_Complex_Vector_Norms.Max_Norm(y);
    put("Max norm of the difference :"); put(nrm,3); new_line;
    Standard_Complex_Matrices.Copy(ym,dm);
    Standard_Complex_Matrices.Sub(dm,zm);
    nrm := Standard_Complex_Matrix_Norms.Max_Norm(dm);
    put("Max norm of the difference in partial derivatives :");
    put(nrm,3); new_line;
    put("Compare the matrix of all partial derivatives ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Comparing the partial derivatives : "); Write(ym,zm);
    end if;
  end Standard_Eval;

  procedure DoblDobl_Eval ( p : in DoblDobl_Polynomial_Vectors.System ) is

  -- DESCRIPTION :
  --   Evaluates p at a random vector, in double double precision.

    x : constant DoblDobl_Complex_Vectors.Vector(1..p.dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,p.dim);
    y,z : DoblDobl_Complex_Vectors.Vector(1..p.nbr);
    dm,ym,zm : DoblDobl_Complex_Matrices.Matrix(y'range,x'range);
    nrm : double_float;
    ans : character;

  begin
    y := DoblDobl_Polynomial_Vectors.Eval(p,x);
    put_line("y : "); put_line(y); new_line;
    DoblDobl_Polynomial_Vectors.Diff(p,x,ym);
    DoblDobl_Polynomial_Vectors.Speel(p,x,z,zm);
    put_line("z : "); put_line(z); new_line;
    DoblDobl_Complex_Vectors.Sub(y,z);
    nrm := hi_part(DoblDobl_Complex_Vector_Norms.Max_Norm(y));
    put("Max norm of the difference :"); put(nrm,3); new_line;
    DoblDobl_Complex_Matrices.Copy(ym,dm);
    DoblDobl_Complex_Matrices.Sub(dm,zm);
    nrm := hi_part(DoblDobl_Complex_Matrix_Norms.Max_Norm(dm));
    put("Max norm of the difference in partial derivatives :");
    put(nrm,3); new_line;
    put("Compare the matrix of all partial derivatives ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Comparing the partial derivatives : "); Write(ym,zm);
    end if;
  end DoblDobl_Eval;

  procedure QuadDobl_Eval ( p : in QuadDobl_Polynomial_Vectors.System ) is

  -- DESCRIPTION :
  --   Evaluates p at a random vector, in quad double precision.

    x : constant QuadDobl_Complex_Vectors.Vector(1..p.dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,p.dim);
    y,z : QuadDobl_Complex_Vectors.Vector(1..p.nbr);
    dm,ym,zm : QuadDobl_Complex_Matrices.Matrix(y'range,x'range);
    nrm : double_float;
    ans : character;

  begin
    y := QuadDobl_Polynomial_Vectors.Eval(p,x);
    put_line("y : "); put_line(y); new_line;
    QuadDobl_Polynomial_Vectors.Diff(p,x,ym);
    QuadDobl_Polynomial_Vectors.Speel(p,x,z,zm);
    put_line("z : "); put_line(z); new_line;
    QuadDobl_Complex_Vectors.Sub(y,z);
    nrm := hihi_part(QuadDobl_Complex_Vector_Norms.Max_Norm(y));
    put("Max norm of the difference :"); put(nrm,3); new_line;
    QuadDobl_Complex_Matrices.Copy(ym,dm);
    QuadDobl_Complex_Matrices.Sub(dm,zm);
    nrm := hihi_part(QuadDobl_Complex_Matrix_Norms.Max_Norm(dm));
    put("Max norm of the difference in partial derivatives :");
    put(nrm,3); new_line;
    put("Compare the matrix of all partial derivatives ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Comparing the partial derivatives : "); Write(ym,zm);
    end if;
  end QuadDobl_Eval;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random polynomial vector of monomials with coefficients
  --   in standard double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    s : Standard_Polynomial_Vectors.Link_to_System;

  begin
    put_line("Testing monomial operations in double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    s := Standard_Random_System(size,dim,expmax,true);
    put_line("a random polynomial system : "); put(s);
    Standard_Eval(s.all);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random polynomial vector of monomials with coefficients
  --   in double double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    s : DoblDobl_Polynomial_Vectors.Link_to_System;

  begin
    put_line("Testing monomial operations in double double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    s := DoblDobl_Random_System(size,dim,expmax,true);
    put_line("a random polynomial vector : "); put(s);
    DoblDobl_Eval(s.all);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random polynomial vector of monomials with coefficients
  --   in quad double precision.

    size,dim : integer32 := 0;
    expmax : natural32 := 0;
    s : QuadDobl_Polynomial_Vectors.Link_to_System;

  begin
    put_line("Testing monomial operations in quad double precision ...");
    put("Give the size of the vector : "); get(size);
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    s := QuadDobl_Random_System(size,dim,expmax,true);
    put_line("a random monomial vector :"); put(s);
    QuadDobl_Eval(s.all);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then lauches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_polvec;
