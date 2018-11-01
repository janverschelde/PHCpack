with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_VecMats;
with Standard_Dense_Vector_Series2;
with Standard_Dense_Vector_Series2_io;   use Standard_Dense_Vector_Series2_io;
with Standard_Dense_Matrix_Series2;
with Standard_Dense_Matrix_Series2_io;   use Standard_Dense_Matrix_Series2_io;
with Random_Series_Vectors;
with Random_Series_Matrices;
with Standard_Interpolating_Series2;
with Standard_Echelon_Forms;             use Standard_Echelon_Forms;

procedure ts_sersin2 is

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

  function Column_Permutations
              ( pivots : Standard_Integer_Vectors.Vector )
              return Standard_Complex_VecMats.VecMat is

  -- DESCRIPTION :
  --   Returns the sequence of column permutation matrices,
  --   corresponding to the sequence of pivot selections in pivots.

    res : Standard_Complex_VecMats.VecMat(pivots'range);

  begin
    for i in pivots'range loop
      declare
        Q : Standard_Integer_Matrices.Matrix(pivots'range,pivots'range);
        cQ : Standard_Complex_Matrices.Matrix(Q'range(1),Q'range(2));
      begin
        for i in Q'range(1) loop
          for j in Q'range(2) loop
            Q(i,j) := 0;
          end loop;
          Q(i,i) := 1;
        end loop;
        if pivots(i) /= i then  -- swap columns i with pivots(i)
          Q(i,i) := 0;
          Q(i,pivots(i)) := 1;
          Q(pivots(i),pivots(i)) := 0;
          Q(pivots(i),i) := 1;
        end if;
        cQ := Standard_Integer2Complex(Q);
        res(i) := new Standard_Complex_Matrices.Matrix'(cQ);
      end;
    end loop;
    return res;
  end Column_Permutations;

  function Multiplier_Matrices
             ( U : Standard_Complex_Matrices.Matrix ) 
             return Standard_Complex_VecMats.VecMat is

  -- DESCRIPTION :
  --   Returns the sequence of multiplier matrices,
  --   using the accumulated multiplier data in U.

    use Standard_Complex_Numbers;

    res : Standard_Complex_VecMats.VecMat(U'range(2));

  begin
    for k in res'range loop
      declare
        M : Standard_Complex_Matrices.Matrix(U'range(1),U'range(2));
      begin
        for i in M'range(1) loop
          for j in M'range(2) loop
            M(i,j) := Create(0.0);
          end loop;
          M(i,i) := Create(1.0);
        end loop;
        for j in U'range(2) loop
          M(k,j) := U(k,j);
        end loop;
        res(k) := new Standard_Complex_Matrices.Matrix'(M);
      end;
    end loop;
    return res;
  end Multiplier_Matrices;

  procedure Standard_Check
              ( A,L,U : in Standard_Complex_Matrices.Matrix;
                rowp,colp,pivs : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given in L the lower triangular echelon form of A,
  --   with multipliers in U, row and column pivots in rowp, colp, pivs,
  --   the L is reconstructed from the data in A, U, rowp, colp and pivs.

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
    R : Standard_Complex_VecMats.VecMat(pivs'range)
      := Column_Permutations(pivs);
    M : Standard_Complex_VecMats.VecMat(U'range(2))
      := Multiplier_Matrices(U);
   
    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;

  begin
    put_line("Checking the echelon form ...");
    put_line("The row permutation matrix : "); put(P);
    put_line("The column permutation matrix : "); put(Q);
    put_line("The structure of the multiplier matrix : ");
    Write_Integer_Matrix(U);
    W := cP*A;
    put_line("The matrix at the start :"); Write_Integer_Matrix(W);
    for k in R'range loop
      put("Permutation matrix "); put(k,1); put_line(" :");
      Write_Integer_Matrix(R(k).all);
      W := W*R(k).all; -- find the pivot and swap if needed
      put_line("After the permutation :"); Write_Integer_Matrix(W);
      W := W*M(k).all; -- eliminate at the right of the pivot
      put("Elimination matrix "); put(k,1); put_line(" :");
      Write_Integer_matrix(M(k).all);
      put_line("After the elimination :"); Write_Integer_Matrix(W);
    end loop;
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
    if checksum + 1.0 /= 1.0 then
      put_line("The computed echelon form : "); Write_Integer_Matrix(L);
      put_line("The reconstructed echelon form : "); Write_Integer_Matrix(W);
    end if;
  end Standard_Check;

  procedure Standard_Solve
              ( dim : in integer32;
                A,L,U : in Standard_Complex_Matrices.Matrix;
                b : in Standard_Complex_Vectors.Vector;
                rowp,colp,pivs : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given in L the lower triangular echelon form of A,
  --   with multipliers in U, row and column pivots in rowp, colp, pivs,
  --   and right hand side vector in b, the system A*x = b is solved,
  --   in standard double precision.

    P : constant Standard_Integer_Matrices.Matrix(rowp'range,rowp'range)
      := Row_Permutation_Matrix(rowp);
    Q : constant Standard_Integer_Matrices.Matrix(colp'range,colp'range)
      := Column_Permutation_Matrix(colp);
    cP : constant Standard_Complex_Matrices.Matrix(P'range(1),P'range(2))
       := Standard_Integer2Complex(P);
    cQ : constant Standard_Complex_Matrices.Matrix(Q'range(1),Q'range(2))
       := Standard_Integer2Complex(Q);
    R : Standard_Complex_VecMats.VecMat(pivs'range)
      := Column_Permutations(pivs);
    M : Standard_Complex_VecMats.VecMat(U'range(2))
      := Multiplier_Matrices(U);
    Pb : Standard_Complex_Vectors.Vector(b'range);
    x,rv,y,z : Standard_Complex_Vectors.Vector(b'range);
    res,checksum : double_float;
    ans : character;
   
    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

  begin
    put_line("Checking the echelon form ...");
    put_line("The row permutation matrix : "); put(P);
    Pb := cP*b;
    Solve_with_Echelon_Form(L,Pb,x);
    put_line("The echelon form : "); put(L);
    put_line("The permuted right hand side vector :"); put_line(Pb);
    put_line("The solution : "); put_line(x);
    rv := Pb - L*x;
    put_line("The residual vector : "); put_line(rv);
    res := Standard_Complex_Vector_Norms.Max_Norm(rv);
    put("The residual : "); put(res,3); new_line;
    put("Continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      z := x;
      for k in reverse R'range loop
        x := M(k).all*x;
        x := R(k).all*x;
      end loop;
      put_line("The solution vector after the transformations :");
      put_line(x);
      rv := b - A*x;
      put_line("The residual vector : "); put_line(rv);
      res := Standard_Complex_Vector_Norms.Max_Norm(rv);
      put("The residual : "); put(res,3); new_line;
      put("Continue ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Multiply_and_Permute(z,U,pivs);
        put_line("After Multiply_and_Permute :");
        put_line(z);
        Solve(dim,A,b,y,false);
        put_line("The outcome of solve :");
        put_line(y);
        put_line("Component by component comparison :");
        checksum := 0.0;
        for i in x'range loop
          put("x("); put(i,1); put(") : "); put(x(i)); new_line;
          put("y("); put(i,1); put(") : "); put(y(i)); new_line;
          put("z("); put(i,1); put(") : "); put(z(i)); new_line;
          checksum := checksum + AbsVal(x(i) - y(i)) + AbsVal(x(i) - z(i));
        end loop;
        put("The check sum of differences : "); put(checksum,3); new_line;
      end if;
    end if;
  end Standard_Solve;

  procedure Standard_Hermite_Laurent
              ( dim : in integer32;
                mat : in Standard_Dense_Matrix_Series2.Matrix;
                rhs : in Standard_Dense_Vector_Series2.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system defined by the Hermite-Laurent
  --   interpolation conditions, in standard double precision.

    use Standard_Interpolating_Series2;

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
    col_ipvt,pivots : Standard_Integer_Vectors.Vector(1..ncols);
    ans : character;

  begin
    put_line("The Hermite-Laurent matrix :"); Write_Integer_Matrix(A);
    put_line("The Hermite-Laurent right hand side vector :");
    put_line(b);
    L := A;
    Lower_Triangular_Echelon_Form(dim,L,U,row_ipvt,col_ipvt,pivots);
    put_line("The matrix in echelon form :"); Write_Integer_Matrix(L);
    put("Continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Standard_Check(A,L,U,row_ipvt,col_ipvt,pivots);
      put("Continue ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Standard_Solve(dim,A,L,U,b,row_ipvt,col_ipvt,pivots);
      end if;
    end if;
  end Standard_Hermite_Laurent;

  procedure Standard_Integer_Test ( deg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with randomly generated 0/1 matrices as coefficients.
  --   Solves the problem in standard double precision.

    mat : constant Standard_Dense_Matrix_Series2.Matrix
        := Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant Standard_Dense_Vector_Series2.Vector
        := Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series2.Vector
        := Standard_Dense_Matrix_Series2.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    Standard_Hermite_Laurent(dim,mat,rhs);
  end Standard_Integer_Test;

  procedure Standard_Random_Test ( deg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension,
  --   with random coefficients on the complex unit circle.
  --   Solves the problem in standard double precision.

    mat : constant Standard_Dense_Matrix_Series2.Matrix
        := Random_Series_Matrices.Random_Matrix_Series(deg,dim);
    sol : constant Standard_Dense_Vector_Series2.Vector
        := Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series2.Vector
        := Standard_Dense_Matrix_Series2.Multiply(mat,sol);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    Standard_Hermite_Laurent(dim,mat,rhs);
  end Standard_Random_Test;

  procedure Read_Standard_Vector_Series
              ( deg,dim : in integer32;
                vec : out Standard_Dense_Vector_Series2.Vector ) is

  -- DESCRIPTION :
  --   Prompts the user for deg+1 vectors of dimension dim,
  --   for the coefficients of a vector series of degree deg,
  --   in standard double precision.

    wrk : Standard_Complex_Vectors.Vector(1..dim);

  begin
    for i in 0..deg loop
      put("Reading a complex "); put(dim,1);
      put("-vector for coefficient "); put(i,1); put_line(" :");
      get(wrk);
      vec.cff(i) := new Standard_Complex_Vectors.Vector'(wrk);
    end loop;
  end Read_Standard_Vector_Series;

  procedure Read_Standard_Matrix_Series
              ( deg,dim : in integer32;
                mat : out Standard_Dense_Matrix_Series2.Matrix ) is

  -- DESCRIPTION :
  --   Prompts the user for deg+1 matrices of dimension dim,
  --   for the coefficients of a matrix series of degree deg,
  --   in standard double precision.

    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in 0..deg loop
      put("Reading a complex "); put(dim,1); put("-by-"); put(dim,1);
      put(" matrix for coefficient "); put(i,1); put_line(" :");
      get(wrk);
      mat.cff(i) := new Standard_Complex_Matrices.Matrix'(wrk);
    end loop;
  end Read_Standard_Matrix_Series;

  procedure Standard_User_Test ( deg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for a matrix series of given degree deg and dimension dim.
  --   Prompts for a corresponding right hand side vector series.
  --   Solves the problem in standard double precision.

    mat : Standard_Dense_Matrix_Series2.Matrix(deg);
    rhs : Standard_Dense_Vector_Series2.Vector(deg);

  begin
    Read_Standard_Matrix_Series(deg,dim,mat);
    put_line("Your matrix series :"); put(mat);
    Read_Standard_Vector_Series(deg,dim,rhs);
    put_line("Your right hand side vector : "); put(rhs);
    Standard_Hermite_Laurent(dim,mat,rhs);
  end Standard_User_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the series.
  --   The degree is the highest exponent in the power series.
  --   The dimension is the number of variables in the series.

    deg,dim : integer32 := 0;
    i01,usr : character;

  begin
    new_line;
    put_line("Testing solving of singular linear series systems.");
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    new_line;
    put("Test on random zero/one matrices ? (y/n) ");
    Ask_Yes_or_No(i01);
    new_line;
    if i01 = 'y' then
      Standard_Integer_Test(deg,dim);
    else
      put("Test on user given matrices ? (y/n) ");
      Ask_Yes_or_No(usr);
      new_line;
      if usr = 'y'
       then Standard_User_Test(deg,dim);
       else Standard_Random_Test(deg,dim);
      end if;
    end if;
  end Main;

begin
  Main;
end ts_sersin2;
