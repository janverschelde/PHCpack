with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Integer_Circuits;          use Standard_Integer_Circuits;

package body Test_Integer_Circuits is

  function Random_Point_Configuration
             ( n,d : in integer32 ) return Matrix is

    res : Matrix(1..d,1..n);
    low,upp : integer32 := 0;
  
  begin
    put("Give lower bound on the coordinates : "); get(low);
    put("Give upper bound on the coordinates : "); get(upp);
    res := Random_Matrix(natural32(d),natural32(n),low,upp);
    return res;
  end Random_Point_Configuration;

  function Prompt_for_Columns ( d : in integer32 ) return Vector is

    k : integer32 := 0;

  begin
    put("Give the dimension in [2,"); put(d,1); put("] : "); get(k);
    declare
      res : Vector(1..k);
    begin
      put("Enter "); put(k,1); put(" indices : "); get(res);
      return res;
    end;
  end Prompt_for_Columns;

  procedure Check_Circuit ( A : in Matrix ) is

    cols : constant Vector := Prompt_for_Columns(A'last(1));
    rnk : constant integer32 := Rank(A,cols);
    ind : integer32 := 0;
    c : Vector(cols'first..cols'last+1);
    v : Vector(A'range(2));

  begin
    put("The rank of the columns "); put(cols);
    put(" is "); put(rnk,1); put_line(".");
    if rnk = cols'last then
      put("Give a column index : "); get(ind);
      c := Circuit(standard_output,A,cols,ind);
      put("The circuit defined by");
      put(cols); put(" and "); put(ind,1);
      put(" is"); put(c); new_line;
      v := Kernel_Vector(A'last(2),ind,cols,c);
      put("The kernel vector is"); put(v);
      put(" with A*v ="); put(A*v); put_line(".");
    end if;
  end Check_Circuit;

  procedure Enumerate_Circuits ( A : in Matrix ) is

    basis : constant Vector := Prompt_for_Columns(A'last(1));
    rnk : constant integer32 := Rank(A,basis);
    cnt : integer32 := 0;

    procedure Check_Circuit ( c : in Vector; k : in integer32;
                              continue : out boolean ) is
    begin
      cnt := cnt + 1;
      put("Circuit "); put(cnt,1); put(" defined by");
      put(basis); put(" and "); put(k,1);
      put(" is"); put(c); new_line;
      declare
        v : constant Vector := Kernel_Vector(A'last(2),k,basis,c);
      begin
        put("The kernel vector is "); put(v);
        put(" with A*v ="); put(A*v); put_line(".");
      end;
      continue := true;
    end Check_Circuit;
    procedure Enumerate is new Enumerate_All_Circuits(Check_Circuit);

  begin
    put("The rank of the columns"); put(basis);
    put(" is "); put(rnk,1); put_line(".");
    if rnk = basis'last
     then Enumerate(A,basis);
    end if;
  end Enumerate_Circuits;

  procedure Enumerate_Bases ( A : in Matrix; d : in integer32 ) is

    procedure Write_Basis ( v : in Vector; continue : out boolean ) is
    begin
      put(v); new_line;
      continue := true;
    end Write_Basis;
    procedure Enumerate is new Enumerate_All_Bases(Write_Basis);

  begin
    Enumerate(A,d);
  end Enumerate_Bases;

  procedure Enumerate_the_Circuits ( A : in Matrix; d : in integer32 ) is

    cnt : integer32 := 0;

    procedure Write_Circuit ( v : in Vector; continue : out boolean ) is

    -- DESCRIPTION :
    --   Writes the circuit defined by the columns in v.

      c : constant Vector(v'range)
        := Circuit(A,v(v'first..v'last-1),v(v'last));
      w : constant Vector 
        := Kernel_Vector(A'last(2),v(v'last),v(v'first..v'last-1),c);

    begin
      if not Is_Zero(c) then
        cnt := cnt + 1;
        put("circuit "); put(cnt,1); put(" :"); put(v); 
        put(" with coefficients"); put(c); new_line;
        put("the kernel vector is"); put(w);
        put(" with A*w ="); put(A*w); put_line(".");
      end if;
      continue := true;
    end Write_Circuit;
    procedure Enumerate is
      new Standard_Integer_Circuits.Enumerate_Circuits(Write_Circuit);

  begin
    Enumerate(A,d);
  end Enumerate_the_Circuits;

  procedure Matrix_of_Circuits ( A : in Matrix; d : in integer32 ) is

    n : constant integer32 := A'last(2);
    m : constant integer32 := Cardinality(n,d);
    B : constant Matrix := Circuits(A,d);

  begin
    put("Expected number of circuits : "); put(m,1); put_line(".");
    put("Found "); put(B'last(2),1); put_line(" circuits.");
    put_line("the matrix of kernel vectors defined by circuits : ");
    put(B,3);
    put_line("A*B = "); put(A*B);
  end Matrix_of_Circuits;

  procedure Circuits ( n,d : in integer32 ) is

    points : constant Matrix(1..d,1..n) := Random_Point_Configuration(n,d);
    ans : character;
    dim : integer32 := 0;

  begin
    put(n,1); put(" random points of dimension "); put(d,1);
    put_line(" :"); put(points);
    new_line;
    put_line("MENU for testing circuit computation :");
    put_line("  1. compute one circuit defined by column selection;");
    put_line("  2. enumerate all circuits sharing the same basis;");
    put_line("  3. show all column selections of full rank;");
    put_line("  4. enumerate all circuits of given dimension;");
    put_line("  5. compute matrix of all circuits of given dimension.");
    put("Type 1, 2, 3, 4, or 5 to select : ");
    Ask_Alternative(ans,"12345");
    new_line;
    if ans = '3' or ans = '4' or ans = '5'
     then put("Give the dimension : "); get(dim);
    end if;
    case ans is
      when '1' => Check_Circuit(points);
      when '2' => Enumerate_Circuits(points);
      when '3' => Enumerate_Bases(points,dim);
      when '4' => Enumerate_the_Circuits(points,dim);
      when '5' => Matrix_of_Circuits(points,dim);
      when others => null;
    end case;
  end Circuits;

  procedure Main is

    n,d : integer32 := 0;

  begin
    new_line;
    put("Give number of points : "); get(n);
    put("Give dimension of the points : "); get(d);
    Circuits(n,d);
  end Main;

end Test_Integer_Circuits;
