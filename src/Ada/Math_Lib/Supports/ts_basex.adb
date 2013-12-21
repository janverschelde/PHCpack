with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;     use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;      use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;     use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;    use Standard_Floating_Matrices_io;
with Standard_Random_Numbers;          use Standard_Random_Numbers;
with Standard_Random_Vectors;          use Standard_Random_Vectors;
with Standard_Random_Matrices;         use Standard_Random_Matrices;
with Basis_Exchanges;                  use Basis_Exchanges;

procedure ts_basex is 

-- DESCRIPTION :
--   Test on the operations in the package Basis_Exchanges.

  tol : constant double_float := 10.0**(-10);

  procedure Initialize
               ( nrows,ncols : in integer32;
                 active : out Standard_Integer_Vectors.Vector;
                 cffmat,invbas : out Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Initializes the coefficient matrix and inverse of the basis so that
  --   they both contain the identity matrix in their first ncols columns.
  --   The rest of the coefficient matrix is left open.
  --   The active indices are the first nrows.

  -- REQUIRED :
  --   cffmat'range(1) = 1..nrows; cffmat'range(2) = 1..ncols;
  --   invbas'range(1) = invbas'range(2) = active'range = 1..nrows.

  begin
    for i in 1..nrows loop
      active(i) := i;
      for j in 1..nrows loop
        if i = j
         then cffmat(i,j) := 1.0; invbas(i,j) := 1.0;
         else cffmat(i,j) := 0.0; invbas(i,j) := 0.0;
        end if;
      end loop;
    end loop;
  end Initialize;

  procedure Random_Initialize
               ( nrows,ncols : in integer32;
                 active : out Standard_Integer_Vectors.Vector;
                 cffmat,invbas : out Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Initializes the coefficient matrix and inverse of the basis so that
  --   they both contain the identity matrix in their first ncols columns.
  --   The rest of the coefficient matrix is filled with random numbers.
  --   The active indices are the first nrows.

  -- REQUIRED :
  --   cffmat'range(1) = 1..nrows; cffmat'range(2) = 1..ncols;
  --   invbas'range(1) = invbas'range(2) = active'range = 1..nrows.

  begin
    Initialize(nrows,ncols,active,cffmat,invbas);
    for i in 1..nrows loop
      for j in 1..ncols-nrows loop
        cffmat(i,nrows+j) := Random;
      end loop;
    end loop;
  end Random_Initialize;

  procedure Interactive_Initialize
               ( nrows,ncols : in integer32;
                 active : out Standard_Integer_Vectors.Vector;
                 cffmat,invbas : out Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Initializes the coefficient matrix and inverse of the basis so that
  --   they both contain the identity matrix in their first ncols columns.
  --   The rest of the coefficient matrix is filled with user input.
  --   The active indices are the first nrows.

  -- REQUIRED :
  --   cffmat'range(1) = 1..nrows; cffmat'range(2) = 1..ncols;
  --   invbas'range(1) = invbas'range(2) = active'range = 1..nrows.

  begin
    Initialize(nrows,ncols,active,cffmat,invbas);
    for j in nrows+1..ncols loop
      put("Give "); put(nrows,1); put(" floats for column ");
      put(j-nrows,1); put_line(" :");
      for i in 1..nrows loop
        get(cffmat(i,j));
      end loop;
    end loop;
  end Interactive_Initialize;

  procedure Write_Product 
              ( binv,cff : in Standard_Floating_Matrices.Matrix;
                cff_cols : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the product of the matrix binv with the coefficients
  --   from the indicated columns on screen.

    accu : double_float;

  begin
    for i in binv'range(1) loop
      for j in binv'range(2) loop
        accu := 0.0;
        for k in binv'range(1) loop
          accu := accu + binv(i,k)*cff(k,cff_cols(j));
        end loop;
        put(accu,3);
      end loop;
      new_line;
    end loop;
  end Write_Product;

  procedure Check_Product
              ( binv,cff : in Standard_Floating_Matrices.Matrix;
                cff_cols : in Standard_Integer_Vectors.Vector;
                bug : out boolean ) is

  -- DESCRIPTION :
  --   This procedure checks whether the product of the matrix binv
  --   with the columns of the coefficient matrix forms the identity.
  --   If this is not the case, then bug is true on return.

    val : double_float;

  begin
    for i in binv'range(1) loop
      for j in binv'range(2) loop 
        val := 0.0;
        for k in binv'range(2) loop
          val := val + binv(i,k)*cff(k,cff_cols(j));
        end loop;
        if i = j
         then bug := (abs(val - 1.0) > tol);
         else bug := (abs(val) > tol);
        end if;
        exit when bug;
      end loop;
      exit when bug;
    end loop;
  end Check_Product;

  procedure Check_and_Report_Product
              ( binv,cff : in Standard_Floating_Matrices.Matrix;
                cff_cols : in Standard_Integer_Vectors.Vector;
                output : in boolean; bug : out boolean ) is

  -- DESCRIPTION :
  --   Does the product check and the reporting of a bug.

  begin
    Check_Product(binv,cff,cff_cols,bug);
    if not bug then
      put("  OK.");
      if output then
        put_line("   The product is ");
        Write_Product(binv,cff,cff_cols);
        put_line("The inverse of the basis : ");
        put(binv,3);
      else
       new_line;
      end if;
    else
      put_line("  BUG found!!!  The product is ");
      Write_Product(binv,cff,cff_cols);
    end if;
  end Check_and_Report_Product;

  procedure Test_Basis ( binv,mat : in Standard_Floating_Matrices.Matrix;
                         mat_cols : in Standard_Integer_Vectors.Vector;
                         output : in boolean; bug : out boolean ) is

  -- DESCRIPTION :
  --   Tests whether the computed matrix binv represents truely the
  --   inverse of the basis for the columns in the coefficient matrix.
  --   The flag bug is set to true if this is not the case.

    colbinv : Standard_Floating_Matrices.Matrix(binv'range(1),binv'range(2));
    sing : boolean;

  begin
    put("The active columns :"); put(mat_cols);
    Check_and_Report_Product(binv,mat,mat_cols,output,bug);
    put("Recompute the basis from the columns :"); put(mat_cols);
    Column_Basis(binv'last(1),mat,mat_cols,colbinv,sing);
    if sing
     then put_line("Active columns form singular matrix!!!"); bug := true;
     else Check_and_Report_Product(colbinv,mat,mat_cols,output,sing);
    end if;
  end Test_Basis;

  procedure Test_Solve ( n,m : in integer32;
                         binv,mat : in Standard_Floating_Matrices.Matrix;
                         rhs : in Standard_Floating_Vectors.Vector;
                         mat_cols : in Standard_Integer_Vectors.Vector;
                         output : in boolean; bug : out boolean ) is

  -- DESCRIPTION :
  --   Tests the solve procedure with the computed basis on selected
  --   columns of the coefficient matrix and right-hand side vector.

    sol : Standard_Floating_Vectors.Vector(1..n) := Solve(n,binv,rhs,mat_cols);
    sing : boolean;

    procedure Test_Solution is

      prd : double_float;

    begin
      for i in sol'range loop
        prd := 0.0;
        for j in mat'range(1) loop
          prd := prd + mat(j,mat_cols(i))*sol(j);
        end loop;
        bug := (abs(prd - rhs(mat_cols(i))) > tol);
        if output or bug then
          put("eval("); put(i,1); put(") : "); put(prd);
          put("  rhs("); put(mat_cols(i),1);
          put(") : "); put(rhs(mat_cols(i)));
          if bug
           then put_line("  BUG!!!");
           else put_line("  Ok");
          end if;
        end if;
        exit when bug;
      end loop;
    end Test_Solution;

  begin
    bug := false;
    if output then
      put("Computing the solution for selected columns");
      put(mat_cols); put_line(" :");
    end if;
    Test_Solution;
    if output then
      put("Recompute the solution for the columns ");
      put(mat_cols); put_line(" :");
    end if;
    Column_Solve(n,mat,mat_cols,rhs,sol,sing);
    if sing
     then put_line("  singular matrix reported.");
     else Test_Solution;
    end if;
  end Test_Solve;

  procedure Enumerate_Bases ( nrows,ncols : in integer32;
                              random_input : in boolean ) is

    cffmat : Standard_Floating_Matrices.Matrix(1..nrows,1..ncols);
    basinv : Standard_Floating_Matrices.Matrix(1..nrows,1..nrows);
    active : Standard_Integer_Vectors.Vector(1..nrows);
    rhs : constant Standard_Floating_Vectors.Vector(1..ncols)
        := Random_Vector(1,ncols);

    function Value ( binv : Standard_Floating_Matrices.Matrix;
                     cols : Standard_Integer_Vectors.Vector )
                   return double_float is

    -- DESCRIPTION :
    --   Returns the sum of binv*rhs, for the selected columns.

      eva : constant Standard_Floating_Vectors.Vector
          := Solve(nrows,binv,rhs,cols);

    begin
      return Standard_Floating_Vectors.Sum(eva);
    end Value;

    procedure Write ( binv : in Standard_Floating_Matrices.Matrix;
                      cols : in Standard_Integer_Vectors.Vector;
                      level : in integer32; continue : out boolean ) is
    begin
      for i in 1..level loop
        put("   ");
      end loop;
      put(level,1); put(" : ");
      put(cols);
      for i in 1..(nrows+1-level) loop
        put("   ");
      end loop;
      put(Value(binv,cols),3); new_line;
      continue := true;
    end Write;
    procedure Write_Bases is
      new Basis_Exchanges.Enumerate_Basis_Inverses(Write);

  begin
    if random_input
     then Random_Initialize(nrows,ncols,active,cffmat,basinv);
     else Interactive_Initialize(nrows,ncols,active,cffmat,basinv);
    end if;
    put_line("The coefficient matrix :"); put(cffmat,3);
    put_line("The initial basis :"); put(basinv,3);
    put("root at "); put(active); put(" value :");
    put(Value(basinv,active),3); new_line;
    Write_Bases(nrows,ncols,cffmat,basinv,active,tol);
  end Enumerate_Bases;

  procedure Test_Update ( nrows,ncols : in integer32;
                          random_input : in boolean ) is

    cffmat : Standard_Floating_Matrices.Matrix(1..nrows,1..ncols);
    basinv : Standard_Floating_Matrices.Matrix(1..nrows,1..nrows);
    active : Standard_Integer_Vectors.Vector(1..nrows);
    bug : boolean := false;

    procedure Test ( binv : in Standard_Floating_Matrices.Matrix;
                     cols : in Standard_Integer_Vectors.Vector;
                     level : in integer32; continue : out boolean ) is

    begin
      Test_Basis(binv,cffmat,cols,false,bug);
      continue := not bug;
    end Test;
    procedure Test_Enumeration is
      new Basis_Exchanges.Enumerate_Basis_Inverses(Test);

  begin
    if random_input
     then Random_Initialize(nrows,ncols,active,cffmat,basinv);
     else Interactive_Initialize(nrows,ncols,active,cffmat,basinv);
    end if;
    put_line("The coefficient matrix :"); put(cffmat,3);
    put_line("The initial basis :"); put(basinv,3);
    put("The active columns :"); put(active); new_line;
    Test_Enumeration(nrows,ncols,cffmat,basinv,active,tol);
    if not bug
     then put_line("Enumeration was successful, no bug found.");
    end if;
  end Test_Update;

  procedure Interactive_Test_Update ( nrows,ncols : in integer32 ) is

    cffmat : Standard_Floating_Matrices.Matrix(1..nrows,1..ncols);
    basinv : Standard_Floating_Matrices.Matrix(1..nrows,1..nrows);
    active,newact : Standard_Integer_Vectors.Vector(1..nrows);
    bug,fail : boolean := false;
    leave_ind,enter_ind : integer32 := 0;

  begin
    Interactive_Initialize(nrows,ncols,active,cffmat,basinv);
    put_line("The coefficient matrix :"); put(cffmat,3);
    put_line("The initial basis :"); put(basinv,3);
    put("The active columns :"); put(active); new_line;
    loop
      put("Give an active variable that will leave (0 to exit loop) : ");
      get(leave_ind);
      exit when (leave_ind = 0);
      put("Give a passive variable that will enter (0 to exit loop) : ");
      get(enter_ind);
      exit when (leave_ind = 0);
      newact := active;
      active(leave_ind) := enter_ind;
      Update(nrows,ncols,basinv,cffmat,newact,leave_ind,enter_ind,tol,fail);
      if fail then
        put_line("Index could not enter the basis.");
        put("The active rows : "); put(active); new_line;
      else
        Test_Basis(basinv,cffmat,newact,true,bug);
      end if;
    end loop;
  end Interactive_Test_Update;

  procedure Random_Test_Initial_Basis ( nrows,ncols : in integer32 ) is

    cffmat : constant Standard_Floating_Matrices.Matrix(1..nrows,1..ncols)
           := Random_Matrix(natural32(nrows),natural32(ncols));
    rhs : constant Standard_Floating_Vectors.Vector(1..ncols)
        := Random_Vector(1,ncols);
    binv : Standard_Floating_Matrices.Matrix(1..nrows,1..nrows);
    cols : Standard_Integer_Vectors.Vector(1..nrows);
    fail : boolean;

  begin
    Initial_Basis(nrows,ncols,cffmat,tol,binv,cols,fail);
    if fail then
      put_line("Failure to produce an initial basis.");
    else
      put_line("An initial basis has been found.");
      Test_Basis(binv,cffmat,cols,true,fail);
      Test_Solve(nrows,ncols,binv,cffmat,rhs,cols,true,fail);
    end if;
  end Random_Test_Initial_Basis;

  procedure Interactive_Test_Initial_Basis ( nrows,ncols : in integer32 ) is

    cffmat : Standard_Floating_Matrices.Matrix(1..nrows,1..ncols);
    rhs : Standard_Floating_Vectors.Vector(1..ncols);
    binv : Standard_Floating_Matrices.Matrix(1..nrows,1..nrows);
    cols : Standard_Integer_Vectors.Vector(1..nrows);
    fail : boolean;

  begin
    put("Give "); put(nrows,1); put("-by-"); put(ncols,1);
    put_line(" coefficient matrix :"); get(cffmat);
    put("Give "); put(ncols,1);
    put_line("-vector as right-handsides :"); get(rhs);
    Initial_Basis(nrows,ncols,cffmat,tol,binv,cols,fail);
    if fail then
      put_line("Failure to produce an initial basis.");
    else
      put_line("An initial basis has been found.");
      put_line("The inverse of the basis : "); put(binv,3);
      Test_Basis(binv,cffmat,cols,true,fail);
      Test_Solve(nrows,ncols,binv,cffmat,rhs,cols,true,fail);
    end if;
  end Interactive_Test_Initial_Basis;

  procedure Main is

    nvars,ncons : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Basis generation from a coefficient matrix.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. list all columns of enumerated bases for random input");
      put_line("  2. list all columns of enumerated bases for user input");
      put_line("  3. enumerate + test all basis inverses of a random matrix");
      put_line("  4. enumerate + test all basis inverses of a given matrix");
      put_line("  5. test interactively the update of basis inverses");
      put_line("  6. random test on the finding of one basis + solve");
      put_line("  7. test interactively the finding of one basis + solve");
      put("Type your answer (0,1,2,3,4,5,6 or 7) : "); get(ans);
      exit when (ans = '0');
      loop
        new_line;
        put("Give the number of rows (= #variables, dim) : "); get(nvars);
        put("Give the number of columns (= #constraints) : "); get(ncons);
        exit when (ncons >= nvars);
        new_line;
        put("#constraints = "); put(ncons,1); put(" < "); put(nvars,1);
        put_line(" = dimension.  Please try again."); 
      end loop;
      case ans is
        when '1' => Enumerate_Bases(nvars,ncons,true);
        when '2' => Enumerate_Bases(nvars,ncons,false);
        when '3' => Test_Update(nvars,ncons,true);
        when '4' => Test_Update(nvars,ncons,false);
        when '5' => Interactive_Test_Update(nvars,ncons);
        when '6' => Random_Test_Initial_Basis(nvars,ncons);
        when '7' => Interactive_Test_Initial_Basis(nvars,ncons);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_basex;
