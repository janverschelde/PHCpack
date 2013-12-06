with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;

package body Standard_Complex_Row_Reduction is

  function Start_Pivots 
              ( n : integer32 ) return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n);

  begin
    for i in 1..n loop
      res(i) := i;
    end loop;
    return res;
  end Start_Pivots;

  function Pivot_in_Row ( A : Matrix; i,k : integer32;
                          tol : double_float ) return integer32 is
  begin
    for j in k..A'last(2) loop
      if AbsVal(A(i,j)) > tol
       then return j;
      end if;
    end loop;
    return 0;
  end Pivot_in_Row;

  procedure Swap_Columns
               ( A : in out Matrix; i,j : in integer32;
                 piv : in out Standard_Integer_Vectors.Vector ) is

    nattmp : constant integer32 := piv(i);
    cmptmp : Complex_Number;

  begin
    piv(i) := piv(j); piv(j) := nattmp;
    for k in A'first(1)..i loop
      cmptmp := A(k,i);
      A(k,i) := A(k,j);
      A(k,j) := cmptmp;
    end loop;
  end Swap_Columns;

  procedure Divide_Row_by_Pivot ( A : in out Matrix; i : in integer32 ) is
  begin
    for j in i+1..A'last(2) loop
      A(i,j) := A(i,j)/A(i,i);
    end loop;
    A(i,i) := Create(1.0);
  end Divide_Row_by_Pivot;

  procedure Eliminate ( A : in out Matrix; i : in integer32;
                        tol : in double_float ) is
  begin
    for k in 1..i-1 loop
      if AbsVal(A(i,k)) > tol then      -- eliminate A(i,k)
        for j in k+1..A'last(2) loop    -- subtract A(i,k) times row k
          A(i,j) := A(i,j) - A(i,k)*A(k,j);
        end loop;
        A(i,k) := Create(0.0);
      end if;
    end loop;
  end Eliminate;

  procedure Reduce_Row 
               ( A : in out Matrix; i : in integer32;
                 pivots : in out Standard_Integer_Vectors.Vector;
                 tol : in double_float; singular : out boolean ) is

    ind : integer32;

  begin
    singular := false;
    if i = 1 then
      ind := Pivot_in_Row(A,1,1,tol);
      if ind = 0 then
        singular := true;
      elsif ind /= i then
        Swap_Columns(A,1,ind,pivots);
      end if;
    else
      Eliminate(A,i,tol);
      ind := Pivot_in_Row(A,i,i,tol);
      if ind = 0 then
        singular := true;
      elsif ind /= i then
        Swap_Columns(A,i,ind,pivots);
      end if;
    end if;
    if not singular
     then Divide_Row_by_Pivot(A,i);
    end if;
  end Reduce_Row;

  procedure Reduce_Row 
               ( file : in file_type; A : in out Matrix; i : in integer32;
                 pivots : in out Standard_Integer_Vectors.Vector;
                 tol : in double_float; singular : out boolean ) is

    ind : integer32;

  begin
    singular := false;
    if i = 1 then
      ind := Pivot_in_Row(A,i,i,tol);
      put(file,"The pivot is row "); put(file,i,1); put(file," is ");
      put(file,ind,1); new_line(file);
      if ind = 0 then
        put(file,"Since row "); put(file,i,1);
        put(file," after column "); put(file,i,1);
        put_line(file," is zero, singular matrix.");
        singular := true;
      elsif ind /= i then
        Swap_Columns(A,i,ind,pivots);
      end if;
    else
      Eliminate(A,i,tol);
      ind := Pivot_in_Row(A,i,i,tol);
      if ind = 0 then
        put_line(file,"Zero row after elimination, singular matrix.");
        singular := true;
      else
        if ind /= i
         then Swap_Columns(A,i,ind,pivots);
        end if;
        put(file,"Row "); put(file,i,1);
        put_line(file," after elimination : ");
        for j in A'range(2) loop
          put(file,A(i,j),3);
        end loop;
        new_line(file);
      end if;
    end if;
    if not singular
     then Divide_Row_by_Pivot(A,i);
    end if;
  end Reduce_Row;

end Standard_Complex_Row_Reduction;
