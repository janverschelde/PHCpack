with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer64_Vectors;
with Standard_Integer64_Vectors_io;      use Standard_Integer64_Vectors_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Lattice_Supports;          use Standard_Lattice_Supports;
with Standard_Integer_Orthogonals;       use Standard_Integer_Orthogonals;
--with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;
--with Multprec_Integer64_Numbers_io;      use Multprec_Integer64_Numbers_io;
--with Multprec_Integer64_Vectors;
--with Multprec_Integer64_Vectors_io;      use Multprec_Integer64_Vectors_io;
--with Multprec_Integer64_Matrices;
--with Multprec_Integer64_Matrices_io;     use Multprec_Integer64_Matrices_io;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Multprec_Integer_Vectors;
with Multprec_Integer_Vectors_io;        use Multprec_Integer_Vectors_io;
with Multprec_Integer_Matrices;
with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;
with Multprec_Lattice_Supports;          use Multprec_Lattice_Supports;
with Multprec_Integer_Orthogonals;       use Multprec_Integer_Orthogonals;

procedure ts_intorto is

-- DESCRIPTION :
--   Interactive development and test of Gram-Schmidt orthogonalization
--   applied to a set of integer vectors.

  function User_Input
             ( n,m : integer32 )
             return Standard_Integer64_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns an n-by-m integer matrix, entered by the user.

    res : Standard_Integer64_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix :"); get(res);
    return res;
  end User_Input;

  function User_Input
             ( n,m : integer32 )
             return Multprec_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns an n-by-m integer matrix, entered by the user.

    res : Multprec_Integer_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix :"); get(res);
    return res;
  end User_Input;

  function Random_Data
             ( n,m : integer32 )
             return Standard_Integer64_Matrices.Matrix is

  -- DESCRIPTION :
  --   Prompts the user for lower and upper bounds on the entries
  --   of a random n-by-m integer matrix.

    res : Standard_Integer64_Matrices.Matrix(1..n,1..m);
    lower,upper : integer64 := 0;

  begin
    put("Give a lower bound for the coordinates : "); get(lower);
    put("Give an upper bound for the coordinates : "); get(upper);
    res := Random_Matrix(natural32(n),natural32(m),lower,upper);
    return res;
  end Random_Data;

  function Convert ( A : Standard_Integer64_Matrices.Matrix )
                   return Multprec_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Converts a standard integer matrix into multiprecision format.

    res : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
    int : integer;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        int := integer(A(i,j));
        res(i,j) := Multprec_Integer_Numbers.Create(int);
      end loop;
    end loop;
    return res;
  end Convert;

  function Random_Data
             ( n,m : integer32 )
             return Multprec_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Prompts the user for lower and upper bounds on the entries
  --   of a random n-by-m integer matrix.

    res : Multprec_Integer_Matrices.Matrix(1..n,1..m);
    A : Standard_Integer64_Matrices.Matrix(1..n,1..m);

  begin
    A := Random_Data(n,m);
    res := Convert(A);
    return res;
  end Random_Data;

  procedure Standard_Test_Orthogonalization
               ( A : in Standard_Integer64_Matrices.Matrix ) is

    B : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2));
    p : integer64;

  begin
    put_line("The given matrix of vectors :"); put(A);
    B := Orthogonalize(A);
    put_line("The orthogonal matrix :"); put(B);
    put_line("Checking orthogonality : ");
    for i in B'range(2) loop
      for j in B'first(2)..i-1 loop
        p := Standard_Lattice_Supports.Inner_Product(B,i,j);
        put("  "); put(i,1); put(" x "); put(j,1);
        put(" = "); put(p,1); new_line;
      end loop;
    end loop;
  end Standard_Test_Orthogonalization;

  procedure Multprec_Test_Orthogonalization
               ( A : in Multprec_Integer_Matrices.Matrix ) is

    B : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
    p : Integer_Number;

  begin
    put_line("The given matrix of vectors :"); put(A);
    B := Orthogonalize(A);
    put_line("The orthogonal matrix :"); put(B);
    put_line("Checking orthogonality : ");
    for i in B'range(2) loop
      for j in B'first(2)..i-1 loop
        p := Multprec_Lattice_Supports.Inner_Product(B,i,j);
        put("  "); put(i,1); put(" x "); put(j,1);
        put(" = "); put(p,1); new_line;
        Clear(p);
      end loop;
    end loop;
  end Multprec_Test_Orthogonalization;

  procedure Standard_Test_Orthogonal_Complement
               ( A : in Standard_Integer64_Matrices.Matrix ) is

    w : Standard_Integer64_Vectors.Vector(A'range(1));

  begin
    put_line("Testing orthogonal complement ...");
    w := Complement(A);
    put("An orthogonal vector to A :"); put(w); new_line;
    put("Inner products with columns of A :");
    for k in A'range(2) loop
      put(" "); put(Standard_Lattice_Supports.Inner_Product(A,k,w),1);
    end loop;
    new_line;
  end Standard_Test_Orthogonal_Complement;

  procedure Multprec_Test_Orthogonal_Complement
               ( A : in Multprec_Integer_Matrices.Matrix ) is

    w : Multprec_Integer_Vectors.Vector(A'range(1));

  begin
    put_line("Testing orthogonal complement ...");
    w := Complement(A);
    put("An orthogonal vector to A :"); put(w); new_line;
    put("Inner products with columns of A :");
    for k in A'range(2) loop
      put(" "); put(Multprec_Lattice_Supports.Inner_Product(A,k,w),1);
    end loop;
    new_line;
  end Multprec_Test_Orthogonal_Complement;

  procedure Standard_Test ( n,m : in integer32 ) is

    A : Standard_Integer64_Matrices.Matrix(1..n,1..m);
    ans : character;

  begin
    put("-> generate random coordinates ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then A := Random_Data(n,m);
     else A := User_Input(n,m);
    end if;
    Standard_Test_Orthogonalization(A);
    Standard_Test_Orthogonal_Complement(A);
  end Standard_Test;

  procedure Multprec_Test ( n,m : in integer32 ) is

    A : Multprec_Integer_Matrices.Matrix(1..n,1..m);
    ans : character;

  begin
    put("-> generate random coordinates ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then A := Random_Data(n,m);
     else A := User_Input(n,m);
    end if;
    put("Skip orthogonalization test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n'
     then Multprec_Test_Orthogonalization(A);
    end if;
    Multprec_Test_Orthogonal_Complement(A);
  end Multprec_Test;

  procedure Main is

    n,m : integer32 := 0;
    ans : character;

  begin
    put("Give the dimension : "); get(n);
    put("Give the number of points : "); get(m);
    put("Use multiprecision arithmetic ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multprec_Test(n,m);
     else Standard_Test(n,m);
    end if;
  end Main;

begin
  Main;
end ts_intorto;
