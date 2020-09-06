with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Random_Matrices;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;  
with Standard_Smith_Normal_Form;
with Multprec_Natural_Numbers;
with Multprec_Integer_Numbers;
with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;
with Multprec_Random_Matrices;
with Multprec_Smith_Normal_Form;

package body Test_Normal_Forms is

  procedure Test_Hermite_Normal_Form 
              ( mat : in Standard_Integer_Matrices.Matrix ) is

    use Standard_Integer_Matrices;

    wrk : Matrix(mat'range(1),mat'range(2));
    u : Matrix(mat'range(1),mat'range(1));
    v : Matrix(mat'range(2),mat'range(2));

  begin
    wrk := mat;
    Upper_Triangulate(u,wrk);
    put_line("The left U-matrix used in Hermite normal form : "); put(u);
    put_line("The matrix in upper triangular form : "); put(wrk);
    put_line("The left U-matrix times original matrix : "); put(u*mat);
    wrk := mat;
    Lower_Triangulate(wrk,v);
    put_line("The right U-matrix used in Hermite normal form : "); put(v);
    put_line("The matrix in lower triangular form : "); put(wrk);
    put_line("The right U-matrix times original matrix : "); put(mat*v);
  end Test_Hermite_Normal_Form;

  procedure Interactive_Test_Hermite_Normal_Form ( n,m : in integer32 ) is

    use Standard_Integer_Matrices;

    mat : Matrix(1..n,1..m);

  begin
    put("Give "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix : "); get(mat);
    put_line("Your matrix :"); put(mat);
    Test_Hermite_Normal_Form(mat);
  end Interactive_Test_Hermite_Normal_Form;

  procedure Test_Standard32_Smith_Normal_Form
              ( mat : in Standard_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

    use Standard_Integer_Matrices;
    use Standard_Smith_Normal_Form;

    wrk : Matrix(mat'range(1),mat'range(2)) := mat;
    n : constant natural32 := natural32(mat'last(1));
    m : constant natural32 := natural32(mat'last(2));
    u : Matrix(mat'range(1),mat'range(1)) := Identity(n);
    v : Matrix(mat'range(2),mat'range(2)) := Identity(m);
    multwrk : Matrix(mat'range(1),mat'range(2));
    rnk : natural32;

  begin
    if output
     then put_line("The matrix to diagonalize : "); put(mat);
    end if;
    Diagonalize(u,wrk,v);
    multwrk := u*mat;
    multwrk := multwrk*v;  
    bug := not ((multwrk = wrk) and Diagonal(wrk));
    if output or bug then
      if bug and not output
       then put_line("The bug input matrix : "); put(mat);
      end if;
      put_line("The left U-matrix used in Smith normal form : "); put(u);
      put_line("The right U-matrix used in Smith normal form : "); put(v);
      put_line("The matrix in diagonal form : "); put(wrk);
      put_line("The left U times original times right V : "); put(multwrk);
    end if;
    if not bug then
      rnk := Rank_of_Diagonal_Matrix(wrk);
      put("rank = "); put(rnk,1); put_line("  okay");
    else
      put_line("Bug detected!");
    end if;
  end Test_Standard32_Smith_Normal_Form;

  procedure Test_Standard64_Smith_Normal_Form
              ( mat : in Standard_Integer64_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

    use Standard_Integer64_Matrices;
    use Standard_Smith_Normal_Form;

    wrk : Matrix(mat'range(1),mat'range(2)) := mat;
    n : constant natural32 := natural32(mat'last(1));
    m : constant natural32 := natural32(mat'last(2));
    u : Matrix(mat'range(1),mat'range(1)) := Identity(n);
    v : Matrix(mat'range(2),mat'range(2)) := Identity(m);
    multwrk : Matrix(mat'range(1),mat'range(2));
    rnk : natural32;

  begin
    if output
     then put_line("The matrix to diagonalize : "); put(mat);
    end if;
    Diagonalize(u,wrk,v);
    multwrk := u*mat;
    multwrk := multwrk*v;  
    bug := not ((multwrk = wrk) and Diagonal(wrk));
    if output or bug then
      if bug and not output
       then put_line("The bug input matrix : "); put(mat);
      end if;
      put_line("The left U-matrix used in Smith normal form : "); put(u);
      put_line("The right U-matrix used in Smith normal form : "); put(v);
      put_line("The matrix in diagonal form : "); put(wrk);
      put_line("The left U times original times right V : "); put(multwrk);
    end if;
    if not bug then
      rnk := Rank_of_Diagonal_Matrix(wrk);
      put("rank = "); put(rnk,1); put_line("  okay");
    else
      put_line("Bug detected!");
    end if;
  end Test_Standard64_Smith_Normal_Form;

  --procedure Check_Empty_Entries 
  --            ( mat : in Multprec_Integer_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   For every entry in the matrix that is empty,
  --   the tuple (i,j) defining that empty entry is written.

  --begin
  --  for i in mat'range(1) loop
  --    for j in mat'range(2) loop
  --      if Multprec_Integer_Numbers.Empty(mat(i,j)) then
  --        put("empty entry at (");
  --        put(i,1); put(","); put(j,1);
  --        put_line(")");
  --      end if;
  --    end loop;
  --  end loop;
  --end Check_Empty_Entries;

  procedure Fix_Empty_Entries
              ( mat : in out Multprec_Integer_Matrices.Matrix;
                output : in boolean ) is
  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        if Multprec_Integer_Numbers.Empty(mat(i,j)) then
          mat(i,j) := Multprec_Integer_Numbers.Create(integer(0));
          if output
           then put("fixed ("); put(i,1); put(","); put(j,1); put_line(")");
          end if;
        else
          declare
            u : constant Multprec_Natural_Numbers.Natural_Number
              := Multprec_Integer_Numbers.Unsigned(mat(i,j));
          begin
            if Multprec_Natural_Numbers.Equal(u,0) then
              -- if Multprec_Integer_Numbers.Sign(mat(i,j)) < 0 then
              --  Multprec_Integer_Numbers.Min(mat(i,j));
              mat(i,j) := Multprec_Integer_Numbers.Create(integer(0));
              if output then
                put("fixed sign at ("); put(i,1); put(","); put(j,1);
                put_line(")");
              end if;
              -- end if;
            end if;
          end;
        end if;
      end loop;
    end loop;
  end Fix_Empty_Entries;

  procedure Test_Multprec_Smith_Normal_Form
              ( mat : in Multprec_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

    use Multprec_Integer_Matrices;
    use Multprec_Smith_Normal_Form;

    wrk : Matrix(mat'range(1),mat'range(2));
    n : constant natural32 := natural32(mat'last(1));
    m : constant natural32 := natural32(mat'last(2));
    u : Matrix(mat'range(1),mat'range(1)) := Identity(n);
    v : Matrix(mat'range(2),mat'range(2)) := Identity(m);
    multwrk : Matrix(mat'range(1),mat'range(2));
    rnk : natural32;

  begin
    Multprec_Integer_Matrices.Copy(mat,wrk);
    if output
     then put_line("The matrix to diagonalize : "); put(mat);
    end if;
    Diagonalize(u,wrk,v);
    Fix_Empty_Entries(u,output);
    Fix_Empty_Entries(v,output);
    Fix_Empty_Entries(wrk,output);
    if output then
      put_line("The left U-matrix used in Smith normal form : "); put(u);
     -- Check_Empty_Entries(u); Fix_Empty_Entries(u);
     -- put_line("after fixing u :"); put(u);
      put_line("The right U-matrix used in Smith normal form : "); put(v);
     -- Check_Empty_Entries(v); Fix_Empty_Entries(v);
     -- put_line("after fixing v :"); put(v);
      put_line("The matrix in diagonal form : "); put(wrk);
     -- Check_Empty_Entries(wrk); Fix_Empty_Entries(wrk);
     -- put_line("after fixing entries :"); put(wrk);
    end if;
    put_line("multiplying u with original matrix ...");
    multwrk := u*mat;     Fix_Empty_Entries(multwrk,output);
    put_line("after fixing, multwrk :"); put(multwrk);
    put_line("multiplying with v ...");
    multwrk := multwrk*v; Fix_Empty_Entries(multwrk,output);
    put_line("after fixing, multwrk :"); put(multwrk);
    put_line("checking equality of matrices and diagonality ...");
    bug := not ((Equal(multwrk,wrk)) and Diagonal(wrk));
    if bug then
      put_line("The bug input matrix : "); put(mat);
      put_line("The matrix after diagonalization : "); put(wrk);
      put_line("The left U times original times right V : "); put(multwrk);
    end if;
    if not bug then
      rnk := Rank_of_Diagonal_Matrix(wrk);
      put("rank = "); put(rnk,1); put_line("  okay");
    else
      put_line("Bug detected!");
    end if;
  end Test_Multprec_Smith_Normal_Form;

  procedure Interactive_Test_Standard_Smith_Normal_Form
              ( n,m : in integer32 ) is

    use Standard_Integer_Matrices;

    mat : Matrix(1..n,1..m);
    bug : boolean;

  begin
    put("Give "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix : "); get(mat);
    Test_Standard32_Smith_Normal_Form(mat,true,bug);
    if bug
     then put_line("Reported a bug...");
    end if;
  end Interactive_Test_Standard_Smith_Normal_Form;

  procedure Interactive_Test_Multprec_Smith_Normal_Form
              ( n,m : in integer32 ) is

    use Multprec_Integer_Matrices;

    mat : Matrix(1..n,1..m);
    bug : boolean;

  begin
    put("Give "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix : "); get(mat);
    Test_Multprec_Smith_Normal_Form(mat,true,bug);
    if bug
     then put_line("Reported a bug...");
    end if;
  end Interactive_Test_Multprec_Smith_Normal_Form;

  procedure Random_Test_Standard32_Smith_Normal_Form ( n,m : in integer32 ) is

    use Standard_Integer_Matrices;
    use Standard_Random_Matrices;

    low,upp : integer32 := 0;
    nbs : natural32 := 0;
    mat : Matrix(1..n,1..m);
    output,bug : boolean;
    ans : character;

  begin
    put("Give lower bound on numbers : "); get(low);
    put("Give upper bound on numbers : "); get(upp);
    put("Give the number of samples : "); get(nbs);
    put("Do you want output during tests ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    for i in 1..nbs loop
      mat := Random_Matrix(natural32(n),natural32(m),low,upp);
      Test_Standard32_Smith_Normal_Form(mat,output,bug);
      exit when bug;
    end loop;
    if not bug
     then put("Successfully tested "); put(nbs,1); put_line(" samples.");
     else put_line("Generated a case that shows there is a bug.");
    end if;
  end Random_Test_Standard32_Smith_Normal_Form;

  procedure Random_Test_Standard64_Smith_Normal_Form ( n,m : in integer32 ) is

    use Standard_Integer64_Matrices;
    use Standard_Random_Matrices;

    low,upp : integer64 := 0;
    nbs : natural32 := 0;
    mat : Matrix(1..n,1..m);
    output,bug : boolean;
    ans : character;

  begin
    put("Give lower bound on numbers : ");
    Standard_Integer_Numbers_io.get(low);
    put("Give upper bound on numbers : ");
    Standard_Integer_Numbers_io.get(upp);
    put("Give the number of samples : "); get(nbs);
    put("Do you want output during tests ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    for i in 1..nbs loop
      mat := Random_Matrix(natural32(n),natural32(m),low,upp);
      Test_Standard64_Smith_Normal_Form(mat,output,bug);
      exit when bug;
    end loop;
    if not bug
     then put("Successfully tested "); put(nbs,1); put_line(" samples.");
     else put_line("Generated a case that shows there is a bug.");
    end if;
  end Random_Test_Standard64_Smith_Normal_Form;

  procedure Random_Test_Standard_Smith_Normal_Form ( n,m : in integer32 ) is

    ans : character;

  begin
    put("Use 64-bit arithmetic? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Random_Test_Standard32_Smith_Normal_Form(n,m);
     else Random_Test_Standard64_Smith_Normal_Form(n,m);
    end if;
  end Random_Test_Standard_Smith_Normal_Form;

  procedure Random_Test_Multprec_Smith_Normal_Form ( n,m : in integer32 ) is

    use Multprec_Integer_Matrices;
    use Multprec_Random_Matrices;

    low,upp : integer32 := 0;
    nbs : natural32 := 0;
    mat : Matrix(1..n,1..m);
    output,bug : boolean;
    ans : character;

  begin
    put("Give lower bound on numbers : "); get(low);
    put("Give upper bound on numbers : "); get(upp);
    put("Give the number of samples : "); get(nbs);
    put("Do you want output during tests ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    for i in 1..nbs loop
      mat := Random_Matrix(natural32(n),natural32(m),low,upp);
      Test_Multprec_Smith_Normal_Form(mat,output,bug);
      exit when bug;
    end loop;
    if not bug
     then put("Successfully tested "); put(nbs,1); put_line(" samples.");
     else put_line("Generated a case that shows there is a bug.");
    end if;
  end Random_Test_Multprec_Smith_Normal_Form;

  procedure Main is

    n,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Test on Hermite and Smith normal form of integer matrices.");
    loop
      new_line;
      put_line("MENU for testing Hermite and Smith normal forms : ");
      put_line("  0. Exit this program.");
      put_line("  1. Hermite normal form on a given matrix.");
      put_line("  2. Smith normal form on a given standard integer matrix.");
      put_line("  3. Smith normal form on a given multiprecision matrix.");
      put_line("  4. Smith normal form on a random standard integer matrix.");
      put_line("  5. Smith normal form on a random multiprecision matrix.");
      put("Type 0, 1, 2, 3, 4, or 5 : "); Ask_Alternative(ans,"012345");
      exit when ans = '0';
      new_line;
      put("Give number of rows : "); get(n);
      put("Give number of columns : "); get(m); 
      case ans is
        when '1' => Interactive_Test_Hermite_Normal_Form(n,m);
        when '2' => Interactive_Test_Standard_Smith_Normal_Form(n,m);
        when '3' => Interactive_Test_Multprec_Smith_Normal_Form(n,m);
        when '4' => Random_Test_Standard_Smith_Normal_Form(n,m);
        when '5' => Random_Test_Multprec_Smith_Normal_Form(n,m);
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Normal_Forms;
