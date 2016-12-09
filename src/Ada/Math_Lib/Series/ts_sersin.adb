with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
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

procedure ts_sersin is

-- DESCRIPTION :
--   Development of solving linear systems of series where the leading
--   coefficient matrices are likely to be singular.

  procedure Write_Integer_Matrix
              ( A : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the integer matrix to screen.

    use Standard_Complex_Numbers;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put(" "); put(integer32(REAL_PART(A(i,j))),1);       
      end loop;
      new_line;
    end loop;
  end Write_Integer_Matrix;

  procedure Write_Integer_Matrix
              ( A : in DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the integer matrix to screen.

    use DoblDobl_Complex_Numbers;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put(" "); put(integer32(hi_part(REAL_PART(A(i,j)))),1);       
      end loop;
      new_line;
    end loop;
  end Write_Integer_Matrix;

  procedure Write_Integer_Matrix
              ( A : in QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the integer matrix to screen.

    use QuadDobl_Complex_Numbers;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put(" "); put(integer32(hihi_part(REAL_PART(A(i,j)))),1);       
      end loop;
      new_line;
    end loop;
  end Write_Integer_Matrix;

  function Is_Zero_Row 
              ( A : Standard_Complex_Matrices.Matrix;
                i : integer32; tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all elements on the i-th row of A
  --   are less than tol in magnitude.  Returns false otherwise.

  begin
    for j in A'range(2) loop
      if Standard_Complex_Numbers.AbsVal(A(i,j)) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero_Row;

  procedure Swap_Rows
              ( A : in out Standard_Complex_Matrices.Matrix;
                i,j : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps row i with row j in A.

    tmp : Standard_Complex_Numbers.Complex_Number;   

  begin
    for k in A'range(2) loop
      tmp := A(i,k);
      A(i,k) := A(j,k);
      A(j,k) := tmp;
    end loop;
  end Swap_Rows;

  procedure Swap_Elements
              ( v : in out Standard_Complex_Vectors.Vector;
                i,j : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps element i with j in v.

    tmp : Standard_Complex_Numbers.Complex_Number;   

  begin
    tmp := v(i);
    v(i) := v(j);
    v(j) := tmp;
  end Swap_Elements;

  procedure Swap_Zero_Rows
              ( A : in out Standard_Complex_Matrices.Matrix;
                b : in out Standard_Complex_Vectors.Vector;
                tol : in double_float; pivrow : out integer32 ) is

  -- DESCRIPTION :
  --   Moves zero rows in A to the top of the matrix,
  --   swapping also the corresponding entries in the right hand side b.

  -- REQUIRED : A'range(1) = b'range.

  -- ON ENTRY :
  --   A        block upper triangular coefficient matrix;
  --   b        corresponding right hand side of a linear system;
  --   tol      tolerance to decide whether a number is zero or not.

  -- ON RETURN :
  --   A        matrix with all zero rows on top;
  --   b        correspondingly swapped right hand side;
  --   pivrow   index of the first nonzero row in A.

    idx : integer32 := A'first(1); -- first nonzero row

  begin
    for i in A'range(1) loop
      if Is_Zero_Row(A,i,tol) then
        if i /= idx then
          Swap_Rows(A,i,idx);
          Swap_Elements(b,i,idx);
        end if;
        idx := idx + 1;
      end if;
    end loop;
    pivrow := idx;
  end Swap_Zero_Rows;

  function Max_on_Row
             ( A : Standard_Complex_Matrices.Matrix;
               i,j : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index k >= j on row i for which A(i,k) is largest.

    use Standard_Complex_Numbers;

    res : integer32 := j;
    maxval : double_float := AbsVal(A(i,j));
    val : double_float;

  begin
    for k in j+1..A'last(2) loop
      val := AbsVal(A(i,k));
      if val > maxval
       then maxval := val; res := k;
      end if;
    end loop;
    return res;
  end Max_on_Row;

  procedure Swap_Columns
              ( A : in out Standard_Complex_Matrices.Matrix;
                ipvt : in out Standard_Integer_Vectors.Vector;
                j,k : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns j and k in A, and the corresponding pivoting
  --   information in ipvt as well.

    Atmp : Standard_Complex_Numbers.Complex_Number;
    itmp : integer32;

  begin
    for i in A'range(1) loop
      Atmp := A(i,j);
      A(i,j) := A(i,k);
      A(i,k) := Atmp;
    end loop;
    itmp := ipvt(j);
    ipvt(j) := ipvt(k);
    ipvt(k) := itmp;
  end Swap_Columns;

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
    tol : constant double_float := 1.0E-12;
    pivrow,pivcol,idx : integer32;
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));

  begin
    put_line("The Hermite-Laurent matrix :");
    Write_Integer_Matrix(A);
    put_line("The Hermite-Laurent right hand side vector :");
    put_line(b);
    for k in ipvt'range loop
      ipvt(k) := k;
    end loop;
    Swap_Zero_Rows(A,b,tol,pivrow);
    put_line("After swapping zero rows :"); 
    Write_Integer_Matrix(A);
    put("The pivot row : "); put(pivrow,1); new_line;
    idx := A'first(2);
    pivcol := Max_on_Row(A,pivrow,idx);
    put("The pivot column : "); put(pivcol,1); new_line; 
    if pivcol /= idx then
      Swap_Columns(A,ipvt,idx,pivcol);
      put_line("After swapping columns : ");
      Write_Integer_Matrix(A);
      put("The pivoting information : "); put(ipvt); new_line;
    end if;
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

  begin
    put_line("The Hermite-Laurent matrix :");
    Write_Integer_Matrix(A);
    put_line("The Hermite-Laurent right hand side vector :");
    put_line(b);
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

  begin
    put_line("The Hermite-Laurent matrix :");
    Write_Integer_Matrix(A);
    put_line("The Hermite-Laurent right hand side vector :");
    put_line(b);
  end QuadDobl_Hermite_Laurent;

  procedure Standard_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
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
  end Standard_Test;

  procedure DoblDobl_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
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
  end DoblDobl_Test;

  procedure QuadDobl_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
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
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the series.
  --   The degree is the highest exponent in the power series.
  --   The dimension is the number of variables in the series.

    deg,dim : integer32 := 0;
    ans : character;

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
    case ans is
      when '0' => Standard_Test(deg,dim);
      when '1' => DoblDobl_Test(deg,dim);
      when '2' => QuadDobl_Test(deg,dim);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sersin;
