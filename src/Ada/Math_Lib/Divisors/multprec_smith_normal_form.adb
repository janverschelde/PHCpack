with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Common_Divisors;           use Multprec_Common_Divisors;

package body Multprec_Smith_Normal_Form is

-- NOTE :
--   The implementation consists in three parts:
--     1. The triangulation of one row and column by left and right
--        multiplication with unimodular matrices.
--     2. Finding pivots and permuting rows and columns.
--     3. Repeatedly invoking steps 1 and 2 yields the Smith normal form.

-- DYNAMIC TRIANGULATORS :

  procedure Upper_Triangulate ( u,a : in out Matrix; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Makes all elements under the current row in the indicated column zero,
  --   by application of unimodular transformations.

  -- REQUIRED : a(row,col) /= 0.

    aa,bb,ka,lb,d,a_rowk,a_ik,tmp : Integer_Number;

  begin
    for i in (row+1)..a'last(1) loop
      if not Multprec_Integer_Numbers.Equal(a(i,col),0) then
                                 -- if a(i,col) /= 0 then make a(i,col) zero
        gcd(a(row,col),a(i,col),ka,lb,d);             -- compute multipliers
        aa := a(row,col)/d;
        bb := a(i,col)/d;
        if Multprec_Integer_Numbers.Equal(aa,bb) then
          if Multprec_Integer_Numbers.Equal(ka,0) then
            Multprec_Integer_Numbers.Copy(lb,ka);
            Multprec_Integer_Numbers.Clear(lb);
            lb := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end if;
        tmp := -bb;
        if Multprec_Integer_Numbers.Equal(aa,tmp) then
          if Multprec_Integer_Numbers.Equal(ka,0) then
            Multprec_Integer_Numbers.Copy(lb,ka);
            Multprec_Integer_Numbers.Min(ka);
            Multprec_Integer_Numbers.Clear(lb);
            lb := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end if;
        for k in a'range(2) loop                      -- perform combinations
          Multprec_Integer_Numbers.Copy(a(row,k),a_rowk);
          Multprec_Integer_Numbers.Copy(a(i,k),a_ik);
         -- a(row,k) :=    ka*a_rowk + lb*a_ik;
          tmp := ka*a_rowk;
          Multprec_Integer_Numbers.Copy(tmp,a(row,k));
          Multprec_Integer_Numbers.Clear(tmp);
          tmp := lb*a_ik;
          Multprec_Integer_Numbers.Add(a(row,k),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          if Multprec_Integer_Numbers.Empty(a(row,k))
           then a(row,k) := Multprec_Integer_Numbers.Create(integer(0));
          end if;
         -- a(i,k)   := (-bb)*a_rowk + aa*a_ik;
          tmp := bb*a_rowk;
          Multprec_Integer_Numbers.Min(tmp);
          Multprec_Integer_Numbers.Copy(tmp,a(i,k));
          Multprec_Integer_Numbers.Clear(tmp);
          tmp := aa*a_ik;
          Multprec_Integer_Numbers.Add(a(i,k),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          Multprec_Integer_Numbers.Clear(a_rowk);
          Multprec_Integer_Numbers.Clear(a_ik);
          if Multprec_Integer_Numbers.Empty(a(i,k))
           then a(i,k) := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end loop;
        for k in u'range(2) loop
          Multprec_Integer_Numbers.Copy(u(row,k),a_rowk);
          Multprec_Integer_Numbers.Copy(u(i,k),a_ik);
         -- u(row,k) :=    ka*a_rowk + lb*a_ik;
          tmp := ka*a_rowk;
          Multprec_Integer_Numbers.Copy(tmp,u(row,k));
          Multprec_Integer_Numbers.Clear(tmp);
          tmp := lb*a_ik;
          Multprec_Integer_Numbers.Add(u(row,k),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          if Multprec_Integer_Numbers.Empty(u(row,k))
           then u(row,k) := Multprec_Integer_Numbers.Create(integer(0));
          end if;
         -- u(i,k)   := (-bb)*a_rowk + aa*a_ik;
          tmp := bb*a_rowk;
          Multprec_Integer_Numbers.Min(tmp);
          Multprec_Integer_Numbers.Copy(tmp,u(i,k));
          Multprec_Integer_Numbers.Clear(tmp);
          tmp := aa*a_ik;
          Multprec_Integer_Numbers.Add(u(i,k),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          Multprec_Integer_Numbers.Clear(a_rowk);
          Multprec_Integer_Numbers.Clear(a_ik);
          if Multprec_Integer_Numbers.Empty(u(i,k))
           then u(i,k) := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end loop;
        Multprec_Integer_Numbers.Clear(aa);
        Multprec_Integer_Numbers.Clear(bb);
        Multprec_Integer_Numbers.Clear(ka);
        Multprec_Integer_Numbers.Clear(lb);
        Multprec_Integer_Numbers.Clear(d);
      end if;
    end loop;
  end Upper_Triangulate;

  procedure Lower_Triangulate ( a,v : in out Matrix; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Makes all elements at the right of the current column in the indicated
  --   row zero, by application of unimodular transformations.

  -- REQUIRED : a(row,col) /= 0.

    aa,bb,ka,lb,d,a_kcol,a_kj,tmp : Integer_Number;

  begin
    for j in (col+1)..a'last(2) loop
      if not Multprec_Integer_Numbers.Equal(a(row,j),0) then
                                  -- if a(row,j) /= 0 then make a(row,j) zero
        gcd(a(row,col),a(row,j),ka,lb,d);              -- compute multipliers
        aa := a(row,col)/d;        
        bb := a(row,j)/d;
        if Multprec_Integer_Numbers.Equal(aa,bb) then
          if Multprec_Integer_Numbers.Equal(ka,0) then
            Multprec_Integer_Numbers.Copy(lb,ka);
            Multprec_Integer_Numbers.Clear(lb);
            lb := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end if;
        tmp := -bb;
        if Multprec_Integer_Numbers.Equal(aa,tmp) then
          if Multprec_Integer_Numbers.Equal(ka,0) then
            Multprec_Integer_Numbers.Copy(lb,ka);
            Multprec_Integer_Numbers.Min(ka);
            Multprec_Integer_Numbers.Clear(lb);
            lb := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end if;
        Multprec_Integer_Numbers.Clear(tmp);
        for k in a'range(1) loop                       -- perform combinations
          Multprec_Integer_Numbers.Copy(a(k,col),a_kcol);
          Multprec_Integer_Numbers.Copy(a(k,j),a_kj);
         -- a(k,col) := a_kcol*ka    + a_kj*lb;
          tmp := a_kcol*ka;
          Multprec_Integer_Numbers.Copy(tmp,a(k,col));
          Multprec_Integer_Numbers.Clear(tmp);
          tmp := a_kj*lb;
          Multprec_Integer_Numbers.Add(a(k,col),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          if Multprec_Integer_Numbers.Empty(a(k,col))
           then a(k,col) := Create(integer(0));
          end if;
         -- a(k,j)   := a_kcol*(-bb) + a_kj*aa;
          tmp := a_kcol*bb;
          Multprec_Integer_Numbers.Min(tmp);
          Multprec_Integer_Numbers.Copy(tmp,a(k,j));
          Multprec_Integer_Numbers.Clear(tmp);
          tmp := a_kj*aa;
          Multprec_Integer_Numbers.Add(a(k,j),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          Multprec_Integer_Numbers.Clear(a_kcol);
          Multprec_Integer_Numbers.Clear(a_kj);
          if Multprec_Integer_Numbers.Empty(a(k,j))
           then a(k,j) := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end loop;
        for k in v'range(1) loop
          Multprec_Integer_Numbers.Copy(v(k,col),a_kcol);
          Multprec_Integer_Numbers.Copy(v(k,j),a_kj);
         -- v(k,col) := a_kcol*ka    + a_kj*lb;
          tmp := a_kcol*ka;
          Multprec_Integer_Numbers.Copy(tmp,v(k,col));
          Multprec_Integer_Numbers.Clear(tmp);
          tmp := a_kj*lb;
          Multprec_Integer_Numbers.Add(v(k,col),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          if Multprec_Integer_Numbers.Empty(v(k,col))
           then v(k,col) := Multprec_Integer_Numbers.Create(integer(0));
          end if;
         -- v(k,j)   := a_kcol*(-bb) + a_kj*aa;
          tmp := a_kcol*bb;
          Multprec_Integer_Numbers.Min(tmp);
          Multprec_Integer_Numbers.Copy(tmp,v(k,j));
          tmp := a_kj*aa;
          Multprec_Integer_Numbers.Add(v(k,j),tmp);
          Multprec_Integer_Numbers.Clear(tmp);
          Multprec_Integer_Numbers.Clear(a_kcol);
          Multprec_Integer_Numbers.Clear(a_kj);
          if Multprec_Integer_NUmbers.Empty(v(k,j))
           then v(k,j) := Multprec_Integer_Numbers.Create(integer(0));
          end if;
        end loop;
        Multprec_Integer_Numbers.Clear(aa);
        Multprec_Integer_Numbers.Clear(bb);
        Multprec_Integer_Numbers.Clear(ka);
        Multprec_Integer_Numbers.Clear(lb);
        Multprec_Integer_Numbers.Clear(d);
      end if;
    end loop;
  end Lower_Triangulate;

-- PIVOTING ROUTINES :

  function Absolute_Value ( i : in Integer_Number ) return Integer_Number is

  -- DESCRIPTION :
  --   Returns the absolute value of the integer.

    res : Integer_Number;

  begin
    Copy(i,res);
    if i < 0
     then Min(res);
    end if;
    return res;
  end Absolute_Value;

  procedure Find_Pivots ( a : in Matrix; pivrow,pivcol : in out integer32 ) is

  -- DESCRIPTION :
  --   Finds the smallest nonzero entry in the matrix a, starting at
  --   the current pivrow and pivcol.

    row : integer32 := pivrow;
    col : integer32 := pivcol;
    first_time : boolean := true;
    min,tmp : Integer_Number;

  begin
    for i in pivrow..a'last(1) loop
      for j in pivcol..a'last(2) loop
        if not Multprec_Integer_Numbers.Equal(a(i,j),0) then
          if first_time then
            min := Absolute_Value(a(i,j));
            row := i; col := j;
            first_time := false;
          else
            tmp := Absolute_Value(a(i,j));
            if tmp < min 
             then Copy(tmp,min); row := i; col := j;
            end if;
            Clear(tmp);
          end if;
        end if;
      end loop;
    end loop;
    pivrow := row; pivcol := col;
  end Find_Pivots;

  procedure Swap ( x,y : in out Integer_Number ) is

  -- DESCRIPTION :
  --   Swaps the numbers x and y.

    temp : Integer_Number;

  begin
    Copy(x,temp);
    Copy(y,x);
    Copy(temp,y);
    Clear(temp);
  end Swap;

  procedure Switch_Rows_and_Columns
              ( a,u,v : in out Matrix;
                pivrow,row,pivcol,col : in integer32 ) is

  -- DESCRIPTION :
  --   Switches the pivot rows and columns in a to the current row and column.
  --   The unimodular matrices u and v are pivoted along.

  -- ON ENTRY :
  --   a        matrix diagonal up to (row,col), (row,col) not included;
  --   u        unimodular matrix for left multiplication;
  --   v        unimodular matrix for right multiplication;
  --   pivrow   pivot row;
  --   row      current row;
  --   pivcol   pivot column;
  --   col      current column.

  -- ON RETURN :
  --   a        matrix with (pivrow,pivcol) at place of (row,col);
  --   u        updated unimodular matrix for left multiplication;
  --   v        updated unimodular matrix for right multiplication.

  begin
    if pivrow /= row then                       -- interchange rows
      for k in a'range(2) loop
        Swap(a(row,k),a(pivrow,k));
      end loop;
      for k in u'range(2) loop
        Swap(u(row,k),u(pivrow,k));
      end loop;
    end if;
    if pivcol /= col then                       -- interchange columns
      for k in a'range(1) loop
        Swap(a(k,col),a(k,pivcol));
      end loop;
      for k in v'range(1) loop
        Swap(v(k,col),v(k,pivcol));
      end loop;
    end if;
  end Switch_Rows_and_Columns;

-- TARGET ROUTINES :

  function Identity ( n : natural32 ) return Matrix is

    res : Matrix(1..integer32(n),1..integer32(n));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := Multprec_Integer_Numbers.Create(integer(1));
         else res(i,j) := Multprec_Integer_Numbers.Create(integer(0));
        end if;
      end loop;
    end loop;
    return res;
  end Identity;

  function Diagonal ( mat : Matrix ) return boolean is
  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        if i /= j then
          if not Multprec_Integer_Numbers.Equal(mat(i,j),0)
           then return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Diagonal;

  function Rank_of_Diagonal_Matrix ( d : matrix ) return natural32 is

    res : natural32 := 0;

  begin
    for i in d'range(1) loop
      exit when i > d'last(2);
      if Multprec_Integer_Numbers.Equal(d(i,i),0)
       then return res;
       else res := res + 1;
      end if;
    end loop;
    return res;
  end Rank_of_Diagonal_Matrix;

  procedure Diagonalize ( u,a,v : in out Matrix ) is

    row : integer32 := a'first(1);
    col : integer32 := a'first(2);
    pivrow,pivcol : integer32;

  begin
    while not Diagonal(a) loop
      row := a'first(1);
      col := a'first(2);
      loop
        pivrow := row; pivcol := col;
        Find_Pivots(a,pivrow,pivcol);
        exit when Multprec_Integer_Numbers.Equal(a(pivrow,pivcol),0);
        Switch_Rows_and_Columns(a,u,v,pivrow,row,pivcol,col);
        Upper_Triangulate(u,a,row,col);
        Lower_Triangulate(a,v,row,col);
        row := row+1;
        col := col+1;
        exit when ((row > a'last(1)) or (col > a'last(2)));
      end loop;
    end loop;
  end Diagonalize;

end Multprec_Smith_Normal_Form;
