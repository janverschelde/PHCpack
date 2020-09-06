with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Random_Matrices;
with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;
with Multprec_Random_Matrices;
with Standard_Smith_Normal_Form;
with Standard_Integer_Matrix_Inverse;
with Multprec_Smith_Normal_Form;
with Multprec_Integer_Matrix_Inverse;

package body Test_Integer_Inverse is

  procedure Test_Standard_Inverse
              ( A : in Standard_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

    use Standard_Integer_Matrices;
  
    d : integer32;
    B : Matrix(A'range(1),A'range(1));
  
  begin
    Standard_Integer_Matrix_Inverse.Inverse(A,d,B);
    put("The determinant : "); put(d,1); new_line;
    if d*d /= 1 then
      put_line("The matrix is not unimodular!"); bug := true;
    else
      if output then
        put_line("The inverse of the original matrix : "); put(B);
      end if;
      bug := not Standard_Integer_Matrix_Inverse.Is_Inverse_Pair(A,B);
      if bug then
        put_line("Product with inverse does not give identity!  Bug!");
      end if;
    end if;
  end Test_Standard_Inverse;

  procedure Test_Standard_Inverse
              ( A : in Standard_Integer64_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

    use Standard_Integer64_Matrices;
  
    d : integer64;
    B : Matrix(A'range(1),A'range(1));
  
  begin
    Standard_Integer_Matrix_Inverse.Inverse(A,d,B);
    put("The determinant : "); 
    Standard_Integer_Numbers_io.put(d,1); new_line;
    if d*d /= 1 then
      put_line("The matrix is not unimodular!"); bug := true;
    else
      if output then
        put_line("The inverse of the original matrix : "); put(B);
      end if;
      bug := not Standard_Integer_Matrix_Inverse.Is_Inverse_Pair(A,B);
      if bug then
        put_line("Product with inverse does not give identity!  Bug!");
      end if;
    end if;
  end Test_Standard_Inverse;

  procedure Test_Multprec_Inverse
              ( A : in Multprec_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

    use Multprec_Integer_Matrices;
  
    d : Integer_Number;
    B : Matrix(A'range(1),A'range(2));
  
  begin
    Multprec_Integer_Matrix_Inverse.Inverse(A,d,B);
    put("The determinant : "); put(d,1); new_line;
    if not Equal(d,1) and not Equal(d,-1) then
      put_line("The matrix is not unimodular!"); bug := true;
    else
      if output then
        put_line("The inverse of the original matrix : "); put(B);
      end if;
      bug := not Multprec_Integer_Matrix_Inverse.Is_Inverse_Pair(A,B);
      if bug then
        put_line("Product with inverse does not give identity!  Bug!");
      end if;
    end if;
  end Test_Multprec_Inverse;

  function Standard_Unimodular_Smith_Transformation
              ( A : Standard_Integer_Matrices.Matrix ) 
              return Standard_Integer_Matrices.Matrix is

    use Standard_Integer_Matrices;

    W : Matrix(A'range(1),A'range(2)) := A;
    U : Matrix(A'range(1),A'range(1))
      := Standard_Smith_Normal_Form.Identity(natural32(A'last(1)));
    V : Matrix(A'range(2),A'range(2))
      := Standard_Smith_Normal_Form.Identity(natural32(A'last(2)));

  begin
    Standard_Smith_Normal_Form.Diagonalize(U,W,V);
    return U*V;
  end Standard_Unimodular_Smith_Transformation;

  function Standard_Unimodular_Smith_Transformation
              ( A : Standard_Integer64_Matrices.Matrix ) 
              return Standard_Integer64_Matrices.Matrix is

    use Standard_Integer64_Matrices;

    W : Matrix(A'range(1),A'range(2)) := A;
    U : Matrix(A'range(1),A'range(1))
      := Standard_Smith_Normal_Form.Identity(natural32(A'last(1)));
    V : Matrix(A'range(2),A'range(2))
      := Standard_Smith_Normal_Form.Identity(natural32(A'last(2)));

  begin
    Standard_Smith_Normal_Form.Diagonalize(U,W,V);
    return U*V;
  end Standard_Unimodular_Smith_Transformation;

  function Multprec_Unimodular_Smith_Transformation
              ( A : Multprec_Integer_Matrices.Matrix ) 
              return Multprec_Integer_Matrices.Matrix is

    use Multprec_Integer_Matrices;

    W : Matrix(A'range(1),A'range(2)) := A;
    U : Matrix(A'range(1),A'range(1))
      := Multprec_Smith_Normal_Form.Identity(natural32(A'last(1)));
    V : Matrix(A'range(2),A'range(2))
      := Multprec_Smith_Normal_Form.Identity(natural32(A'last(2)));

  begin
    Multprec_Smith_Normal_Form.Diagonalize(U,W,V);
    return U*V;
  end Multprec_Unimodular_Smith_Transformation;

  procedure Inverse_of_Standard32_Random_Matrix ( n : in integer32 ) is

    m : integer32 := 0;
    A : Standard_Integer_Matrices.Matrix(1..n,1..n);
    U : Standard_Integer_Matrices.Matrix(1..n,1..n);
    ans : character;
    output,bug : boolean;
    low,upp : integer32 := 0;

  begin
    put("Give lower bound on numbers : "); get(low);
    put("Give upper bound on numbers : "); get(upp);
    new_line;
    put("Give the number of tests : "); get(m);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    for i in 1..m loop
      A := Standard_Random_Matrices.Random_Matrix
             (natural32(n),natural32(n),low,upp);
      U := Standard_Unimodular_Smith_Transformation(A);
      if output
       then put_line("A random unimodular matrix : "); put(U);
      end if;
      Test_Standard_Inverse(U,output,bug);
      if bug
       then put_line("bug found for the matrix "); put(U);
      end if;
      exit when bug;
    end loop;
    if not bug
     then put("Tested "); put(m,1); put(" cases successfully.");
    end if;
  end Inverse_of_Standard32_Random_Matrix;

  procedure Inverse_of_Standard64_Random_Matrix ( n : in integer32 ) is

    m : integer32 := 0;
    A : Standard_Integer64_Matrices.Matrix(1..n,1..n);
    U : Standard_Integer64_Matrices.Matrix(1..n,1..n);
    ans : character;
    output,bug : boolean;
    low,upp : integer64 := 0;

  begin
    put("Give lower bound on numbers : ");
    Standard_Integer_Numbers_io.get(low);
    put("Give upper bound on numbers : "); 
    Standard_Integer_Numbers_io.get(upp);
    new_line;
    put("Give the number of tests : "); get(m);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    for i in 1..m loop
      A := Standard_Random_Matrices.Random_Matrix
             (natural32(n),natural32(n),low,upp);
      U := Standard_Unimodular_Smith_Transformation(A);
      if output
       then put_line("A random unimodular matrix : "); put(U);
      end if;
      Test_Standard_Inverse(U,output,bug);
      if bug
       then put_line("bug found for the matrix "); put(U);
      end if;
      exit when bug;
    end loop;
    if not bug
     then put("Tested "); put(m,1); put(" cases successfully.");
    end if;
  end Inverse_of_Standard64_Random_Matrix;

  procedure Inverse_of_Standard_Random_Matrix ( n : in integer32 ) is

    ans : character;

  begin
    put("Use 64-bit arithmetic? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Inverse_of_Standard64_Random_Matrix(n);
     else Inverse_of_Standard32_Random_Matrix(n);
    end if;
  end Inverse_of_Standard_Random_Matrix;

  procedure Inverse_of_Multprec_Random_Matrix ( n : in integer32 ) is

    m : integer32 := 0;
    A : Multprec_Integer_Matrices.Matrix(1..n,1..n);
    U : Multprec_Integer_Matrices.Matrix(1..n,1..n);
    ans : character;
    output,bug : boolean;
    low,upp : integer32 := 0;

  begin
    put("Give lower bound on numbers : "); get(low);
    put("Give upper bound on numbers : "); get(upp);
    new_line;
    put("Give the number of tests : "); get(m);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    for i in 1..m loop
      A := Multprec_Random_Matrices.Random_Matrix
            (natural32(n),natural32(n),low,upp);
      U := Multprec_Unimodular_Smith_Transformation(A);
      if output
       then put_line("A random unimodular matrix : "); put(U);
      end if;
      Test_Multprec_Inverse(U,output,bug);
      if bug
       then put_line("bug found for the matrix "); put(U);
      end if;
      exit when bug;
    end loop;
    if not bug
     then put("Tested "); put(m,1); put(" cases successfully.");
    end if;
  end Inverse_of_Multprec_Random_Matrix;

  procedure Inverse_of_Given_Standard32_Matrix ( n : in integer32 ) is

    A : Standard_Integer_Matrices.Matrix(1..n,1..n);
    bug : boolean;

  begin
    put("Give a "); put(n,1); put("-by-"); put(n,1);
    put_line(" matrix : "); get(A);
    put_line("Your matrix :"); put(A);
    Test_Standard_Inverse(A,true,bug);
    if not bug
     then put_line("inverses checked out fine, okay");
     else put_line("bug detected");
    end if;
  end Inverse_of_Given_Standard32_Matrix;

  procedure Inverse_of_Given_Standard64_Matrix ( n : in integer32 ) is

    A : Standard_Integer64_Matrices.Matrix(1..n,1..n);
    bug : boolean;

  begin
    put("Give a "); put(n,1); put("-by-"); put(n,1);
    put_line(" matrix : "); get(A);
    put_line("Your matrix :"); put(A);
    Test_Standard_Inverse(A,true,bug);
    if not bug
     then put_line("inverses checked out fine, okay");
     else put_line("bug detected");
    end if;
  end Inverse_of_Given_Standard64_Matrix;

  procedure Inverse_of_Given_Standard_Matrix ( n : in integer32 ) is

    ans : character;

  begin
    put("Use 64-bit arithmetic? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Inverse_of_Given_Standard64_Matrix(n);
     else Inverse_of_Given_Standard32_Matrix(n);
    end if;
  end Inverse_of_Given_Standard_Matrix;

  procedure Inverse_of_Given_Multprec_Matrix ( n : in integer32 ) is

    A : Multprec_Integer_Matrices.Matrix(1..n,1..n);
    bug : boolean;

  begin
    put("Give a "); put(n,1); put("-by-"); put(n,1);
    put_line(" matrix : "); get(A);
    put_line("Your matrix :"); put(A);
    Test_Multprec_Inverse(A,true,bug);
    if not bug
     then put_line("inverses checked out fine, okay");
     else put_line("bug detected");
    end if;
  end Inverse_of_Given_Multprec_Matrix;

  procedure Main is

    n : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("MENU for testing the inverse of an integer matrix : ");
    put_line("  1. inverse of standard integer random unimodular matrix;");
    put_line("  2. inverse of multiprecision random unimodular matrix;");
    put_line("  3. give your own matrix, with standard integer arithmetic;");
    put_line("  4. give your own matrix, with multiprecision arithmetic.");
    put("Type 1, 2, 3, or 4 to select: "); Ask_Alternative(ans,"1234");
    new_line;
    put("Give the dimension : "); get(n);
    case ans is
      when '1' => Inverse_of_Standard_Random_Matrix(n);
      when '2' => Inverse_of_Multprec_Random_Matrix(n);
      when '3' => Inverse_of_Given_Standard_Matrix(n);
      when '4' => Inverse_of_Given_Multprec_Matrix(n);
      when others => null;
    end case;
  end Main;

end Test_Integer_Inverse;
