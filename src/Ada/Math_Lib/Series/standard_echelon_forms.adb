with text_io;                          use text_io;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;      use Standard_Integer_Vectors_io;
with Standard_Complex_Numbers;

package body Standard_Echelon_Forms is

  procedure Write_Integer_Matrix
              ( A : in Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put(" "); put(integer32(REAL_PART(A(i,j))),2);
      end loop;
      new_line;
    end loop;
  end Write_Integer_Matrix;

  function Is_Zero_Row 
              ( A : Standard_Complex_Matrices.Matrix;
                i : integer32; tol : double_float ) return boolean is
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
               i,j : integer32; tol : double_float ) return integer32 is

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
    if maxval > tol
     then return res;
     else return -1;
    end if;
  end Max_on_Row;

  procedure Swap_Columns
              ( A : in out Standard_Complex_Matrices.Matrix;
                ipvt : in out Standard_Integer_Vectors.Vector;
                j,k : in integer32 ) is

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

  procedure Eliminate_on_Row
              ( A : in out Standard_Complex_Matrices.Matrix;
                i,j : in integer32; tol : in double_float ) is

     use Standard_Complex_Numbers;

     fac : Complex_Number;

  begin
    for k in j+1..A'last(2) loop
      if AbsVal(A(i,k)) > tol then
        fac := A(i,k)/A(i,j);
        for row in i..A'last(1) loop
          A(row,k) := A(row,k) - fac*A(row,j);
        end loop;
      end if;
    end loop;
  end Eliminate_on_Row;

  procedure Lower_Triangular_Echelon_Form
              ( A : in out Standard_Complex_Matrices.Matrix;
                b : in out Standard_Complex_Vectors.Vector;
                verbose : in boolean := true ) is

    tol : constant double_float := 1.0E-12;
    pivrow,pivcol,colidx : integer32;
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));

  begin
    for k in ipvt'range loop
      ipvt(k) := k;
    end loop;
    Swap_Zero_Rows(A,b,tol,pivrow);
    if verbose then
      put_line("After swapping zero rows :"); Write_Integer_Matrix(A);
    end if;
    colidx := A'first(2);
    loop
      pivcol := Max_on_Row(A,pivrow,colidx,tol);
      if verbose then
        put("The pivot row : "); put(pivrow,1); 
        put("  pivot column : "); put(pivcol,1); 
        put("  column index : "); put(colidx,1); new_line;
      end if;
      if pivcol /= -1 then -- if no pivot, then skip row
        if pivcol /= colidx then
          Swap_Columns(A,ipvt,colidx,pivcol);
          if verbose then
            put_line("After swapping columns : "); Write_Integer_Matrix(A);
            put("The pivoting information : "); put(ipvt); new_line;
          end if;
        end if;
        Eliminate_on_Row(A,pivrow,colidx,tol);
        if verbose then
          put_line("After elimination on the pivot row :");
          Write_Integer_Matrix(A);
        end if;
        colidx := colidx + 1;
      end if;
      pivrow := pivrow + 1;
      exit when ((pivrow > A'last(1)) or (colidx > A'last(2)));
    end loop;
  end Lower_Triangular_Echelon_Form;

end Standard_Echelon_Forms;
