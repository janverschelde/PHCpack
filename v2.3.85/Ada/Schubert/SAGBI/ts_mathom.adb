with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Complex_Matrices;           use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;        use Standard_Complex_Matrices_io;
with Matrix_Homotopies;

procedure ts_mathom is

-- DESCRIPTION :
--   Test of the package Matrix_Homotopies.

  procedure Read is

    n,m : integer32 := 0;

  begin
    new_line;
    put_line("Adding a new matrix homotopy");
    new_line;
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    declare
      start,target : Matrix(1..n,1..m);
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" complex matrix as start : "); get(start);
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" complex matrix as target : "); get(target);
      Matrix_Homotopies.Add(start,target);
      put_line("Matrix homotopy added...");
	end;
  end Read;

  procedure Modify is

    k : natural32 := 0;
    n,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Modifying or adding a matrix homotopy.");
    new_line;
    put("Give the number of the map : "); get(k);
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    declare
      start,target : Matrix(1..n,1..m);
    begin
      put("Do you want to modify the start ? (y/n) "); get(ans);
      if ans = 'y' then
        put("Give "); put(n,1); put("x"); put(m,1);
        put_line(" complex matrix as start : "); get(start);
        Matrix_Homotopies.Add_Start(k,start);
      end if;
      put("Do you want to modify the target ? (y/n) "); get(ans);
      if ans = 'y' then
        put("Give "); put(n,1); put("x"); put(m,1);
        put_line(" complex matrix as target : "); get(target);
        Matrix_Homotopies.Add_Target(k,target);
      end if;
      put_line("Matrix homotopy modified...");
    end;
  end Modify;

  procedure Eval is

    k : natural32 := 0;
    t : Complex_Number;
    ans : character;

  begin
    new_line;
    put_line("Evaluating a matrix homotopy.");
    new_line;
    loop
      put("Give the number of the map : "); get(k);
      loop
        put("Give the continuation parameter t : "); get(t);
        declare
          eva : constant Matrix := Matrix_Homotopies.Eval(k,t);
        begin
          put_line("The evaluated map : "); put(eva);
        end;
        put("Do you want to evaluate this maps for other values ? (y/n) ");
        get(ans);
        exit when (ans /= 'y');
      end loop;
      put("Do you want to evaluate other maps ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Eval;

  procedure Main is

    nb : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing the matrix homotopies.");
    new_line;
    put("Give the number of matrix homotopies : "); get(nb);
    Matrix_Homotopies.Init(nb);
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. Exit this program;");
      put_line("  1. Read and Add a new matrix homotopy;");
      put_line("  2. Modify or Add a matrix homotopy;");
      put_line("  3. Evaluate a matrix homotopy.");
      put("Type 0,1,2 or 3 to make your choice : "); get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Read;
        when '2' => Modify;
        when '3' => Eval;
        when others => null;
      end case;
    end loop;
    Matrix_Homotopies.Clear;
  end Main;

begin
  Main;
end ts_mathom;
