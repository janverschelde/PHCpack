with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Common_Divisors;
with Standard64_Common_Divisors;

package body Standard_Smith_Normal_Form is

-- NOTE :
--   The implementation consists in three parts:
--     1. The triangulation of one row and column by left and right
--        multiplication with unimodular matrices.
--     2. Finding pivots and permuting rows and columns.
--     3. Repeatedly invoking steps 1 and 2 yields the Smith normal form.

-- DYNAMIC TRIANGULATORS :

  procedure Upper_Triangulate
              ( U,A : in out Standard_Integer_Matrices.Matrix;
                row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Makes all elements under the current row in the indicated column zero,
  --   by application of unimodular transformations.

  -- REQUIRED : A(row,col) /= 0.

    aa,bb,ka,lb,d,a_rowk,a_ik : integer32;
    use Standard_Common_Divisors;

  begin
    for i in (row+1)..A'last(1) loop
      if A(i,col) /= 0 then                           -- make A(i,col) zero
        gcd(A(row,col),A(i,col),ka,lb,d);             -- compute multipliers
        aa := A(row,col)/d;     bb := A(i,col)/d;
        if (aa = bb) and then (ka = 0)
         then ka := lb; lb := 0;
        end if;
        if (aa = -bb) and then (ka = 0)
         then ka := -lb; lb := 0;
        end if;
        for k in A'range(2) loop                      -- perform combinations
          a_rowk := A(row,k);
          a_ik := A(i,k);
          A(row,k) :=    ka*a_rowk + lb*a_ik;
          A(i,k)   := (-bb)*a_rowk + aa*a_ik;
        end loop;
        for k in U'range(2) loop
          a_rowk := U(row,k);
          a_ik := U(i,k);
          U(row,k) :=    ka*a_rowk + lb*a_ik;
          U(i,k)   := (-bb)*a_rowk + aa*a_ik;
        end loop;
      end if;
    end loop;
  end Upper_Triangulate;

  procedure Upper_Triangulate
              ( U,A : in out Standard_Integer64_Matrices.Matrix;
                row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Makes all elements under the current row in the indicated column zero,
  --   by application of unimodular transformations.

  -- REQUIRED : A(row,col) /= 0.

    use Standard64_Common_Divisors;
    aa,bb,ka,lb,d,a_rowk,a_ik : integer64;

  begin
    for i in (row+1)..A'last(1) loop
      if A(i,col) /= 0 then                           -- make A(i,col) zero
        gcd(A(row,col),A(i,col),ka,lb,d);             -- compute multipliers
        aa := A(row,col)/d;     bb := A(i,col)/d;
        if (aa = bb) and then (ka = 0)
         then ka := lb; lb := 0;
        end if;
        if (aa = -bb) and then (ka = 0)
         then ka := -lb; lb := 0;
        end if;
        for k in A'range(2) loop                      -- perform combinations
          a_rowk := A(row,k);
          a_ik := A(i,k);
          A(row,k) :=    ka*a_rowk + lb*a_ik;
          A(i,k)   := (-bb)*a_rowk + aa*a_ik;
        end loop;
        for k in U'range(2) loop
          a_rowk := U(row,k);
          a_ik := U(i,k);
          U(row,k) :=    ka*a_rowk + lb*a_ik;
          U(i,k)   := (-bb)*a_rowk + aa*a_ik;
        end loop;
      end if;
    end loop;
  end Upper_Triangulate;

  procedure Lower_Triangulate
              ( A,V : in out Standard_Integer_Matrices.Matrix;
                row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Makes all elements at the right of the current column in the indicated
  --   row zero, by application of unimodular transformations.

  -- REQUIRED : A(row,col) /= 0.

    aa,bb,ka,lb,d,a_kcol,a_kj : integer32;
    use Standard_Common_Divisors;

  begin
    for j in (col+1)..A'last(2) loop
      if A(row,j) /= 0 then                            -- make A(row,j) zero
        gcd(A(row,col),A(row,j),ka,lb,d);              -- compute multipliers
        aa := A(row,col)/d;        bb := A(row,j)/d;
        if (aa = bb) and then (ka = 0)
         then ka := lb; lb := 0;
        end if;
        if (aa = -bb) and then (ka = 0)
         then ka := -lb; lb := 0;
        end if;
        for k in A'range(1) loop                       -- perform combinations
          a_kcol := A(k,col);
          a_kj := A(k,j);
          A(k,col) := a_kcol*ka    + a_kj*lb;
          A(k,j)   := a_kcol*(-bb) + a_kj*aa;
        end loop;
        for k in V'range(1) loop
          a_kcol := V(k,col);
          a_kj := V(k,j);
          V(k,col) := a_kcol*ka    + a_kj*lb;
          V(k,j)   := a_kcol*(-bb) + a_kj*aa;
        end loop;
      end if;
    end loop;
  end Lower_Triangulate;

  procedure Lower_Triangulate
              ( A,V : in out Standard_Integer64_Matrices.Matrix;
                row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Makes all elements at the right of the current column in the indicated
  --   row zero, by application of unimodular transformations.

  -- REQUIRED : A(row,col) /= 0.

    use Standard64_Common_Divisors;
    aa,bb,ka,lb,d,a_kcol,a_kj : integer64;

  begin
    for j in (col+1)..A'last(2) loop
      if A(row,j) /= 0 then                            -- make A(row,j) zero
        gcd(A(row,col),A(row,j),ka,lb,d);              -- compute multipliers
        aa := A(row,col)/d;        bb := A(row,j)/d;
        if (aa = bb) and then (ka = 0)
         then ka := lb; lb := 0;
        end if;
        if (aa = -bb) and then (ka = 0)
         then ka := -lb; lb := 0;
        end if;
        for k in A'range(1) loop                       -- perform combinations
          a_kcol := A(k,col);
          a_kj := A(k,j);
          A(k,col) := a_kcol*ka    + a_kj*lb;
          A(k,j)   := a_kcol*(-bb) + a_kj*aa;
        end loop;
        for k in V'range(1) loop
          a_kcol := V(k,col);
          a_kj := V(k,j);
          V(k,col) := a_kcol*ka    + a_kj*lb;
          V(k,j)   := a_kcol*(-bb) + a_kj*aa;
        end loop;
      end if;
    end loop;
  end Lower_Triangulate;

-- PIVOTING ROUTINES :

  function Absolute_Value ( i : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the absolute value of the integer.

  begin
    if i < 0
     then return -i;
     else return i;
    end if;
  end Absolute_Value;

  function Absolute_Value ( i : integer64 ) return integer64 is

  -- DESCRIPTION :
  --   Returns the absolute value of the integer.

  begin
    if i < 0
     then return -i;
     else return i;
    end if;
  end Absolute_Value;

  procedure Find_Pivots
              ( A : in Standard_Integer_Matrices.Matrix;
                pivrow,pivcol : in out integer32 ) is

  -- DESCRIPTION :
  --   Finds the smallest nonzero entry in the matrix A,
  --   starting at the current pivrow and pivcol.

    row : integer32 := pivrow;
    col : integer32 := pivcol;
    first_time : boolean := true;
    min : integer32;

  begin
    for i in pivrow..A'last(1) loop
      for j in pivcol..A'last(2) loop
        if (A(i,j) /= 0) then
          if first_time then
            min := Absolute_Value(A(i,j));
            row := i; col := j;
            first_time := false;
          elsif (Absolute_Value(A(i,j)) < min) then
            min := Absolute_Value(A(i,j));
            row := i; col := j;
          end if;
        end if;
      end loop;
    end loop;
    pivrow := row;
    pivcol := col;
  end Find_Pivots;

  procedure Find_Pivots
              ( A : in Standard_Integer64_Matrices.Matrix;
                pivrow,pivcol : in out integer32 ) is

  -- DESCRIPTION :
  --   Finds the smallest nonzero entry in the matrix A,
  --   starting at the current pivrow and pivcol.

    row : integer32 := pivrow;
    col : integer32 := pivcol;
    first_time : boolean := true;
    min : integer64;

  begin
    for i in pivrow..A'last(1) loop
      for j in pivcol..A'last(2) loop
        if (A(i,j) /= 0) then
          if first_time then
            min := Absolute_Value(A(i,j));
            row := i; col := j;
            first_time := false;
          elsif (Absolute_Value(A(i,j)) < min) then
            min := Absolute_Value(A(i,j));
            row := i; col := j;
          end if;
        end if;
      end loop;
    end loop;
    pivrow := row;
    pivcol := col;
  end Find_Pivots;

  procedure Switch_Rows_and_Columns
              ( A,U,V : in out Standard_Integer_Matrices.Matrix;
                pivrow,row,pivcol,col : in integer32 ) is

  -- DESCRIPTION :
  --   Switches the pivot rows and columns in a to the current row and column.
  --   The unimodular matrices u and v are pivoted along.

  -- ON ENTRY :
  --   A        matrix diagonal up to (row,col), (row,col) not included;
  --   U        unimodular matrix for left multiplication;
  --   V        unimodular matrix for right multiplication;
  --   pivrow   pivot row;
  --   row      current row;
  --   pivcol   pivot column;
  --   col      current column.

  -- ON RETURN :
  --   A        matrix with (pivrow,pivcol) at place of (row,col);
  --   U        updated unimodular matrix for left multiplication;
  --   V        updated unimodular matrix for right multiplication.

    temp : integer32;

  begin
    if pivrow /= row then                       -- interchange rows
      for k in A'range(2) loop
        temp := A(row,k);
        A(row,k) := A(pivrow,k);
        A(pivrow,k) := temp;
      end loop;
      for k in U'range(2) loop
        temp := U(row,k);
        U(row,k) := U(pivrow,k);
        U(pivrow,k) := temp;
      end loop;
    end if;
    if pivcol /= col then                       -- interchange columns
      for k in A'range(1) loop
        temp := A(k,col);
        A(k,col) := A(k,pivcol);
        A(k,pivcol) := temp;
      end loop;
      for k in V'range(1) loop
        temp := V(k,col);
        V(k,col) := V(k,pivcol);
        V(k,pivcol) := temp;
      end loop;
    end if;
  end Switch_Rows_and_Columns;

  procedure Switch_Rows_and_Columns
              ( A,U,V : in out Standard_Integer64_Matrices.Matrix;
                pivrow,row,pivcol,col : in integer32 ) is

  -- DESCRIPTION :
  --   Switches the pivot rows and columns in a to the current row and column.
  --   The unimodular matrices u and v are pivoted along.

  -- ON ENTRY :
  --   A        matrix diagonal up to (row,col), (row,col) not included;
  --   U        unimodular matrix for left multiplication;
  --   V        unimodular matrix for right multiplication;
  --   pivrow   pivot row;
  --   row      current row;
  --   pivcol   pivot column;
  --   col      current column.

  -- ON RETURN :
  --   A        matrix with (pivrow,pivcol) at place of (row,col);
  --   U        updated unimodular matrix for left multiplication;
  --   V        updated unimodular matrix for right multiplication.

    temp : integer64;

  begin
    if pivrow /= row then                       -- interchange rows
      for k in A'range(2) loop
        temp := A(row,k);
        A(row,k) := A(pivrow,k);
        A(pivrow,k) := temp;
      end loop;
      for k in U'range(2) loop
        temp := U(row,k);
        U(row,k) := U(pivrow,k);
        U(pivrow,k) := temp;
      end loop;
    end if;
    if pivcol /= col then                       -- interchange columns
      for k in A'range(1) loop
        temp := A(k,col);
        A(k,col) := A(k,pivcol);
        A(k,pivcol) := temp;
      end loop;
      for k in V'range(1) loop
        temp := V(k,col);
        V(k,col) := V(k,pivcol);
        V(k,pivcol) := temp;
      end loop;
    end if;
  end Switch_Rows_and_Columns;

-- TARGET ROUTINES :

  function Identity
             ( n : natural32 ) return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(1..integer32(n),1..integer32(n));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := 1;
         else res(i,j) := 0;
        end if;
      end loop;
    end loop;
    return res;
  end Identity;

  function Identity
             ( n : natural32 ) return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(1..integer32(n),1..integer32(n));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := 1;
         else res(i,j) := 0;
        end if;
      end loop;
    end loop;
    return res;
  end Identity;

  function Diagonal ( A : Standard_Integer_Matrices.Matrix ) return boolean is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if i /= j then
          if A(i,j) /= 0
           then return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Diagonal;

  function Diagonal
             ( A : Standard_Integer64_Matrices.Matrix ) return boolean is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if i /= j then
          if A(i,j) /= 0
           then return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Diagonal;

  function Rank_of_Diagonal_Matrix
             ( D : Standard_Integer_Matrices.Matrix ) return natural32 is

    res : natural32 := 0;

  begin
    for i in D'range(1) loop
      exit when i > D'last(2);
      if D(i,i) = 0
       then return res;
       else res := res + 1;
      end if;
    end loop;
    return res;
  end Rank_of_Diagonal_Matrix;

  function Rank_of_Diagonal_Matrix
             ( D : Standard_Integer64_Matrices.Matrix ) return natural32 is

    res : natural32 := 0;

  begin
    for i in D'range(1) loop
      exit when i > D'last(2);
      if D(i,i) = 0
       then return res;
       else res := res + 1;
      end if;
    end loop;
    return res;
  end Rank_of_Diagonal_Matrix;

  procedure Diagonalize ( U,A,V : in out Standard_Integer_Matrices.Matrix ) is

    row : integer32 := A'first(1);
    col : integer32 := A'first(2);
    pivrow,pivcol : integer32;

  begin
    while not Diagonal(A) loop
      row := a'first(1);
      col := a'first(2);
      loop
        pivrow := row; pivcol := col;
        Find_Pivots(A,pivrow,pivcol);
        exit when (A(pivrow,pivcol) = 0);
        Switch_Rows_and_Columns(A,U,V,pivrow,row,pivcol,col);
        Upper_Triangulate(U,A,row,col);
        Lower_Triangulate(A,V,row,col);
        row := row+1;
        col := col+1;
        exit when ((row > a'last(1)) or (col > a'last(2)));
      end loop;
    end loop;
  end Diagonalize;

  procedure Diagonalize
              ( U,A,V : in out Standard_Integer64_Matrices.Matrix ) is

    row : integer32 := A'first(1);
    col : integer32 := A'first(2);
    pivrow,pivcol : integer32;

  begin
    while not Diagonal(A) loop
      row := a'first(1);
      col := a'first(2);
      loop
        pivrow := row; pivcol := col;
        Find_Pivots(A,pivrow,pivcol);
        exit when (A(pivrow,pivcol) = 0);
        Switch_Rows_and_Columns(A,U,V,pivrow,row,pivcol,col);
        Upper_Triangulate(U,A,row,col);
        Lower_Triangulate(A,V,row,col);
        row := row+1;
        col := col+1;
        exit when ((row > a'last(1)) or (col > a'last(2)));
      end loop;
    end loop;
  end Diagonalize;

end Standard_Smith_Normal_Form;
