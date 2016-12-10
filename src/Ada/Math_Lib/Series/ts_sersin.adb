with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Matrices;
with Standard_Dense_Vector_Series;
with Standard_Dense_Vector_Series_io;    use Standard_Dense_Vector_Series_io;
with Standard_Random_Series;
with DoblDobl_Random_Series;
with QuadDobl_Random_Series;
with Standard_Dense_Matrix_Series;
with Standard_Dense_Matrix_Series_io;    use Standard_Dense_Matrix_Series_io;
with DoblDobl_Dense_Vector_Series;
with DoblDobl_Dense_Vector_Series_io;    use DoblDobl_Dense_Vector_Series_io;
with DoblDobl_Dense_Matrix_Series;
with DoblDobl_Dense_Matrix_Series_io;    use DoblDobl_Dense_Matrix_Series_io;
with QuadDobl_Dense_Vector_Series;
with QuadDobl_Dense_Vector_Series_io;    use QuadDobl_Dense_Vector_Series_io;
with QuadDobl_Dense_Matrix_Series;
with QuadDobl_Dense_Matrix_Series_io;    use QuadDobl_Dense_Matrix_Series_io;
with Random_Matrix_Series;               use Random_Matrix_Series;
with Standard_Interpolating_Series;
with DoblDobl_Interpolating_Series;
with QuadDobl_Interpolating_Series;
with Standard_Echelon_Forms;             use Standard_Echelon_Forms;
with DoblDobl_Echelon_Forms;             use DoblDobl_Echelon_Forms;
with QuadDobl_Echelon_Forms;             use QuadDobl_Echelon_Forms;

procedure ts_sersin is

-- DESCRIPTION :
--   Development of solving linear systems of series where the leading
--   coefficient matrices are likely to be singular.

  function Standard_Integer2Complex
             ( A : Standard_Integer_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Converts the integer matrix to a complex matrix,
  --   in standard double precision.

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Standard_Complex_Numbers.Create(A(i,j));
      end loop;
    end loop;
    return res;
  end Standard_Integer2Complex;

  function DoblDobl_Integer2Complex
             ( A : Standard_Integer_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Converts the integer matrix to a complex matrix,
  --   in double double precision.

    res : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := DoblDobl_Complex_Numbers.Create(A(i,j));
      end loop;
    end loop;
    return res;
  end DoblDobl_Integer2Complex;

  function QuadDobl_Integer2Complex
             ( A : Standard_Integer_Matrices.Matrix )
             return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Converts the integer matrix to a complex matrix,
  --   in quad double precision.

    res : QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := QuadDobl_Complex_Numbers.Create(A(i,j));
      end loop;
    end loop;
    return res;
  end QuadDobl_Integer2Complex;

  procedure Standard_Check
              ( A,L,U : in Standard_Complex_Matrices.Matrix;
                rowp,colp : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given in L the lower triangular echelon form of A,
  --   with multipliers in U, row and column pivots in rowp, colp,
  --   the L is reconstructed from the data in A, U, rowp and colp.

    P : constant Standard_Integer_Matrices.Matrix(rowp'range,rowp'range)
      := Row_Permutation_Matrix(rowp);
    Q : constant Standard_Integer_Matrices.Matrix(colp'range,colp'range)
      := Column_Permutation_Matrix(colp);
    cP : constant Standard_Complex_Matrices.Matrix(P'range(1),P'range(2))
       := Standard_Integer2Complex(P);
    cQ : constant Standard_Complex_Matrices.Matrix(Q'range(1),Q'range(2))
       := Standard_Integer2Complex(Q);
    W : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));
    checksum : double_float := 0.0;
   
    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;

  begin
    put_line("Checking the echelon form ...");
    put_line("The row permutation matrix : "); put(P);
    put_line("The column permutation matrix : "); put(Q);
    put_line("The structure of the multiplier matrix : ");
    Write_Integer_Matrix(U);
    W := cP*A;
    W := W*cQ;
    W := W*U;
    for i in L'range(1) loop
      for j in L'range(2) loop
        put("L("); put(i,1); put(","); put(j,1); put(") :");
        put(L(i,j)); new_line;
        put("W("); put(i,1); put(","); put(j,1); put(") :");
        put(W(i,j)); new_line;
        checksum := checksum + AbsVal(L(i,j) - W(i,j));
      end loop;
    end loop;
    put("Check sum of differences : "); put(checksum); new_line;
  end Standard_Check;

  procedure DoblDobl_Check
              ( A,L,U : in DoblDobl_Complex_Matrices.Matrix;
                rowp,colp : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given in L the lower triangular echelon form of A,
  --   with multipliers in U, row and column pivots in rowp, colp,
  --   the L is reconstructed from the data in A, U, rowp and colp.

    P : constant Standard_Integer_Matrices.Matrix(rowp'range,rowp'range)
      := Row_Permutation_Matrix(rowp);
    Q : constant Standard_Integer_Matrices.Matrix(colp'range,colp'range)
      := Column_Permutation_Matrix(colp);
    cP : constant DoblDobl_Complex_Matrices.Matrix(P'range(1),P'range(2))
       := DoblDobl_Integer2Complex(P);
    cQ : constant DoblDobl_Complex_Matrices.Matrix(Q'range(1),Q'range(2))
       := DoblDobl_Integer2Complex(Q);
    W : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));
    checksum : double_double := create(0.0);

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Matrices;

  begin
    put_line("Checking the echelon form ...");
    put_line("The row permutation matrix : "); put(P);
    put_line("The column permutation matrix : "); put(Q);
    put_line("The structure of the multiplier matrix : ");
    Write_Integer_Matrix(U);
    W := cP*A;
    W := W*cQ;
    W := W*U;
    for i in L'range(1) loop
      for j in L'range(2) loop
        put("L("); put(i,1); put(","); put(j,1); put(") : ");
        put(L(i,j)); new_line;
        put("W("); put(i,1); put(","); put(j,1); put(") : ");
        put(W(i,j)); new_line;
        checksum := checksum + AbsVal(L(i,j) - W(i,j));
      end loop;
    end loop;
    put("Check sum of differences : "); put(checksum); new_line;
  end DoblDobl_Check;

  procedure QuadDobl_Check
              ( A,L,U : in QuadDobl_Complex_Matrices.Matrix;
                rowp,colp : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given in L the lower triangular echelon form of A,
  --   with multipliers in U, row and column pivots in rowp, colp,
  --   the L is reconstructed from the data in A, U, rowp and colp.

    P : constant Standard_Integer_Matrices.Matrix(rowp'range,rowp'range)
      := Row_Permutation_Matrix(rowp);
    Q : constant Standard_Integer_Matrices.Matrix(colp'range,colp'range)
      := Column_Permutation_Matrix(colp);
    cP : constant QuadDobl_Complex_Matrices.Matrix(P'range(1),P'range(2))
       := QuadDobl_Integer2Complex(P);
    cQ : constant QuadDobl_Complex_Matrices.Matrix(Q'range(1),Q'range(2))
       := QuadDobl_Integer2Complex(Q);
    W : QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));
    checksum : quad_double := create(0.0);

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Matrices;

  begin
    put_line("Checking the echelon form ...");
    put_line("The row permutation matrix : "); put(P);
    put_line("The column permutation matrix : "); put(Q);
    put_line("The structure of the multiplier matrix : ");
    Write_Integer_Matrix(U);
    W := cP*A;
    W := W*cQ;
    W := W*U;
    for i in L'range(1) loop
      for j in L'range(2) loop
        put("L("); put(i,1); put(","); put(j,1); put(") : ");
        put(L(i,j)); new_line;
        put("W("); put(i,1); put(","); put(j,1); put(") : ");
        put(W(i,j)); new_line;
        checksum := checksum + AbsVal(L(i,j) - W(i,j));
      end loop;
    end loop;
    put("Check sum of differences : "); put(checksum); new_line;
  end QuadDobl_Check;

  procedure Standard_Hermite_Laurent
              ( mat : in Standard_Dense_Matrix_Series.Matrix;
                rhs : in Standard_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system defined by the Hermite-Laurent
  --   interpolation conditions, in standard double precision.

    use Standard_Interpolating_Series;

    deg : constant integer32 := mat.deg;
    nbr : constant integer32 := mat.cff(0)'last(1);
    nbc : constant integer32 := mat.cff(0)'last(2);
    nrows : constant integer32 := nbr*(2*deg+1);
    ncols : constant integer32 := nbc*(2*deg+1);
    A : Standard_Complex_Matrices.Matrix(1..nrows,1..ncols)
      := Hermite_Laurent_Matrix(mat.cff(0..deg));
    b : Standard_Complex_Vectors.Vector(1..nrows)
      := Hermite_Laurent_Vector(rhs.cff(0..deg));
    L,U : Standard_Complex_Matrices.Matrix(1..nrows,1..ncols);
    row_ipvt : Standard_Integer_Vectors.Vector(1..nrows);
    col_ipvt : Standard_Integer_Vectors.Vector(1..ncols);

  begin
    put_line("The Hermite-Laurent matrix :"); Write_Integer_Matrix(A);
    put_line("The Hermite-Laurent right hand side vector :");
    put_line(b);
    L := A;
    Lower_Triangular_Echelon_Form(L,U,row_ipvt,col_ipvt);
    put_line("The matrix in echelon form :"); Write_Integer_Matrix(L);
    Standard_Check(A,L,U,row_ipvt,col_ipvt);
  end Standard_Hermite_Laurent;

  procedure DoblDobl_Hermite_Laurent
              ( mat : in DoblDobl_Dense_Matrix_Series.Matrix;
                rhs : in DoblDobl_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system defined by the Hermite-Laurent
  --   interpolation conditions, in double double precision.

    use DoblDobl_Interpolating_Series;

    deg : constant integer32 := mat.deg;
    nbr : constant integer32 := mat.cff(0)'last(1);
    nbc : constant integer32 := mat.cff(0)'last(2);
    nrows : constant integer32 := nbr*(2*deg+1);
    ncols : constant integer32 := nbc*(2*deg+1);
    A : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols)
      := Hermite_Laurent_Matrix(mat.cff(0..deg));
    b : DoblDobl_Complex_Vectors.Vector(1..nrows)
      := Hermite_Laurent_Vector(rhs.cff(0..deg));
    L,U : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    row_ipvt : Standard_Integer_Vectors.Vector(1..nrows);
    col_ipvt : Standard_Integer_Vectors.Vector(1..ncols);

  begin
    put_line("The Hermite-Laurent matrix :");
    Write_Integer_Matrix(A);
    put_line("The Hermite-Laurent right hand side vector :");
    put_line(b);
    L := A;
    Lower_Triangular_Echelon_Form(L,U,row_ipvt,col_ipvt);
    put_line("The matrix in echelon form :"); Write_Integer_Matrix(L);
    DoblDobl_Check(A,L,U,row_ipvt,col_ipvt);
  end DoblDobl_Hermite_Laurent;

  procedure QuadDobl_Hermite_Laurent
              ( mat : in QuadDobl_Dense_Matrix_Series.Matrix;
                rhs : in QuadDobl_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system defined by the Hermite-Laurent
  --   interpolation conditions, in quad double precision.

    use QuadDobl_Interpolating_Series;

    deg : constant integer32 := mat.deg;
    nbr : constant integer32 := mat.cff(0)'last(1);
    nbc : constant integer32 := mat.cff(0)'last(2);
    nrows : constant integer32 := nbr*(2*deg+1);
    ncols : constant integer32 := nbc*(2*deg+1);
    A : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..ncols)
      := Hermite_Laurent_Matrix(mat.cff(0..deg));
    b : QuadDobl_Complex_Vectors.Vector(1..nrows)
      := Hermite_Laurent_Vector(rhs.cff(0..deg));
    L,U : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    row_ipvt : Standard_Integer_Vectors.Vector(1..nrows);
    col_ipvt : Standard_Integer_Vectors.Vector(1..ncols);

  begin
    put_line("The Hermite-Laurent matrix :");
    Write_Integer_Matrix(A);
    put_line("The Hermite-Laurent right hand side vector :");
    put_line(b);
    L := A;
    Lower_Triangular_Echelon_Form(L,U,row_ipvt,col_ipvt);
    put_line("The matrix in echelon form :"); Write_Integer_Matrix(L);
    QuadDobl_Check(A,L,U,row_ipvt,col_ipvt);
  end QuadDobl_Hermite_Laurent;

  procedure Standard_Integer_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with randomly generated 0/1 matrices as coefficients.
  --   Solves the problem in standard double precision.

    mat : constant Standard_Dense_Matrix_Series.Matrix
        := Standard_Random_Matrix_Series(deg,dim,0,1);
    sol : constant Standard_Dense_Vector_Series.Vector
        := Standard_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series.Vector
        := Standard_Dense_Matrix_Series.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    Standard_Hermite_Laurent(mat,rhs);
  end Standard_Integer_Test;

  procedure DoblDobl_Integer_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with randomly generated 0/1 matrices as coefficients.
  --   Solves the problem in standard double precision.

    mat : constant DoblDobl_Dense_Matrix_Series.Matrix
        := DoblDobl_Random_Matrix_Series(deg,dim,0,1);
    sol : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Dense_Matrix_Series.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    DoblDobl_Hermite_Laurent(mat,rhs);
  end DoblDobl_Integer_Test;

  procedure QuadDobl_Integer_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with randomly generated 0/1 matrices as coefficients.
  --   Solves the problem in standard double precision.

    mat : constant QuadDobl_Dense_Matrix_Series.Matrix
        := QuadDobl_Random_Matrix_Series(deg,dim,0,1);
    sol : constant QuadDobl_Dense_Vector_Series.Vector
        := QuadDobl_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant QuadDobl_Dense_Vector_Series.Vector
        := QuadDobl_Dense_Matrix_Series.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    QuadDobl_Hermite_Laurent(mat,rhs);
  end QuadDobl_Integer_Test;

  procedure Standard_Random_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with random coefficients on the complex unit circle.
  --   Solves the problem in standard double precision.

    mat : constant Standard_Dense_Matrix_Series.Matrix
        := Standard_Random_Matrix_Series(deg,dim);
    sol : constant Standard_Dense_Vector_Series.Vector
        := Standard_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series.Vector
        := Standard_Dense_Matrix_Series.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    Standard_Hermite_Laurent(mat,rhs);
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with random coefficients on the complex unit circle.
  --   Solves the problem in standard double precision.

    mat : constant DoblDobl_Dense_Matrix_Series.Matrix
        := DoblDobl_Random_Matrix_Series(deg,dim);
    sol : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Dense_Matrix_Series.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    DoblDobl_Hermite_Laurent(mat,rhs);
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with random coefficients on the complex unit circle.
  --   Solves the problem in standard double precision.

    mat : constant QuadDobl_Dense_Matrix_Series.Matrix
        := QuadDobl_Random_Matrix_Series(deg,dim);
    sol : constant QuadDobl_Dense_Vector_Series.Vector
        := QuadDobl_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant QuadDobl_Dense_Vector_Series.Vector
        := QuadDobl_Dense_Matrix_Series.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    QuadDobl_Hermite_Laurent(mat,rhs);
  end QuadDobl_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the series.
  --   The degree is the highest exponent in the power series.
  --   The dimension is the number of variables in the series.

    deg,dim : integer32 := 0;
    ans,i01 : character;

  begin
    new_line;
    put_line("Testing solving of singular linear series systems.");
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    put("Test on zero/one matrices ? (y/n) ");
    Ask_Yes_or_No(i01);
    new_line;
    if i01 = 'y' then
      case ans is
        when '0' => Standard_Integer_Test(deg,dim);
        when '1' => DoblDobl_Integer_Test(deg,dim);
        when '2' => QuadDobl_Integer_Test(deg,dim);
        when others => null;
      end case;
    else
      case ans is
        when '0' => Standard_Random_Test(deg,dim);
        when '1' => DoblDobl_Random_Test(deg,dim);
        when '2' => QuadDobl_Random_Test(deg,dim);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_sersin;
