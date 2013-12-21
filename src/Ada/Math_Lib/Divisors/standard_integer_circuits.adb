with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer_Linear_Solvers;

package body Standard_Integer_Circuits is

  function Submatrix ( A : Matrix; columns : Vector ) return Matrix is

    res : Matrix(A'range(1),columns'range);

  begin
    for j in columns'range loop
      for i in A'range(1) loop
        res(i,j) := A(i,columns(j));
      end loop;
    end loop;
    return res;
  end Submatrix;

  function Rank ( A : Matrix; columns : Vector ) return integer32 is

    B : constant Matrix(A'range(1),columns'range) := Submatrix(A,columns);

  begin
    return Standard_Integer_Linear_Solvers.rank(B);
  end Rank;

  function Circuit ( A : Matrix; basis : Vector; k : integer32 )
                   return Vector is

    res : Vector(basis'first..basis'last+1);
    cols : Vector(res'range);
    B : Matrix(A'range(1),res'range);

  begin
    cols(basis'range) := basis;
    cols(cols'last) := k;
    B := Submatrix(A,cols);
    Standard_Integer_Linear_Solvers.Upper_Triangulate(B);
    Standard_Integer_Linear_Solvers.Scale(B);
    res := (res'range => 0);
    Standard_Integer_Linear_Solvers.Solve0(B,res);
    return res;
  end Circuit;

  function Circuit ( file : file_type;
                     A : Matrix; basis : Vector; k : integer32 )
                   return Vector is

    res : Vector(basis'first..basis'last+1);
    cols : Vector(res'range);
    B : Matrix(A'range(1),res'range);

  begin
    cols(basis'range) := basis;
    cols(cols'last) := k;
    B := Submatrix(A,cols);
    put_line(file,"the submatrix B :"); put(file,B);
    Standard_Integer_Linear_Solvers.Upper_Triangulate(B);
    put_line(file,"the triangulated submatrix B :"); put(file,B);
    Standard_Integer_Linear_Solvers.Scale(B);
    put_line(file,"the scaled triangulated submatrix B :"); put(file,B);
    res := (res'range => 0);
    Standard_Integer_Linear_Solvers.Solve0(B,res);
    put(file,"the solution :"); put(file,res); new_line(file);
    return res;
  end Circuit;

  function Kernel_Vector ( n,k : integer32; b,c : Vector ) return Vector is

    res : Vector(1..n) := (1..n => 0);

  begin
    for i in b'range loop
      res(b(i)) := c(i);
    end loop;
    res(k) := c(c'last);
    return res;
  end Kernel_Vector;

  function Is_In ( b : Vector; k : integer32 ) return boolean is
  begin
    for i in b'range loop
      if b(i) = k
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_Zero ( v : Vector ) return boolean is
  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  procedure Enumerate_All_Circuits ( A : in Matrix; basis : in Vector ) is

    continue : boolean := true;

  begin
    for k in A'range(2) loop
      if not Is_In(basis,k) then
        declare
          c : constant Vector := Circuit(A,basis,k);
        begin
          if not Is_Zero(c) then
            Report(c,k,continue);
            exit when not continue;
          end if;
        end;
      end if;
    end loop;
  end Enumerate_All_Circuits;

  procedure Enumerate_All_Bases ( A : in Matrix; d : in integer32 ) is

    buffer : Vector(1..d);
    contin : boolean := true;

    procedure Enumerate ( k : in integer32 ) is

    -- DESCRIPTION :
    --   Sets the element at the k-th position in buffer.
  
      start : integer32;

    begin
      if k > buffer'last then
        if Rank(A,buffer) = d
         then Report(buffer,contin);
        end if;
      else
        if k = 1
         then start := 1;
         else start := buffer(k-1)+1;
        end if;
        for i in start..A'last(2) loop
          if not Is_In(buffer(1..k-1),i) then
            buffer(k) := i;
            Enumerate(k+1);
            exit when not contin;
          end if;
        end loop;
      end if;
    end Enumerate;

  begin
    Enumerate(1);
  end Enumerate_All_Bases;

  procedure Enumerate_Circuits ( A : in Matrix; d : in integer32 ) is

    procedure Circuit_at_Basis ( b : in Vector; cont : out boolean ) is

    -- DESCRIPTION :
    --   Computes all circuits with as basis the columns of A indexed by b.

      contin : boolean := true;
      cols : Vector(b'first..b'last+1);

    begin
      cols(b'range) := b;
      for k in b(b'last)+1..A'last(2) loop
        cols(cols'last) := k;
        Report(cols,contin);
        exit when not contin;
      end loop;
      cont := contin;
    end Circuit_at_Basis;
    procedure Enumerate is new Enumerate_All_Bases(Circuit_at_Basis);

  begin
    Enumerate(A,d);
  end Enumerate_Circuits;

  function Cardinality ( n,d : integer32) return integer32 is

    p : integer32 := n;
    q : integer32 := 1;
    f : constant integer32 := (n-d);

  begin
    for i in (d+1)..(n-1) loop
      p := p*i;
    end loop;
    for i in 2..(n-d) loop
      q := q*i;
    end loop;
    return ((p/q)*f);
  end Cardinality;

  function Circuits ( A : Matrix; d : integer32 ) return List is

    res,res_last : List;

    procedure Append_Circuit ( v : in Vector; continue : out boolean ) is

    -- DESCRIPTION :
    --   Appends the kernel vector of the columns defined by v to 
    --   the list res.

      c : constant Vector(v'range)
        := Circuit(A,v(v'first..v'last-1),v(v'last));
      w : constant Vector 
        := Kernel_Vector(A'last(2),v(v'last),v(v'first..v'last-1),c);

    begin
      if not Is_Zero(c)
       then Append(res,res_last,w);
      end if;
      continue := true;
    end Append_Circuit;
    procedure Enumerate is
      new Standard_Integer_Circuits.Enumerate_Circuits(Append_Circuit);

  begin
    Enumerate(A,d);
    return res;
  end Circuits;

  function List_to_Matrix ( L : List ) return Matrix is

    m : constant integer32 := integer32(Length_Of(L));
    lv : Link_to_Vector := Head_Of(L);
    n : constant integer32 := lv'last;
    res : Matrix(1..n,1..m);
    tmp : List := L;

  begin
    for j in 1..m loop
      lv := Head_Of(tmp);
      for i in lv'range loop
        res(i,j) := lv(i);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end List_to_Matrix;

  function Circuits ( A : Matrix; d : integer32 ) return Matrix is

    c : List := Circuits(A,d);
    B : constant Matrix := List_to_Matrix(c);

  begin
    Clear(c);
    return B;
  end Circuits;

end Standard_Integer_Circuits;
