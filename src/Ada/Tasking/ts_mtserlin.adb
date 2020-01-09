with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_VecMats;
with Standard_Complex_VecMats_io;        use Standard_Complex_VecMats_io;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Vectors_io; use Standard_Complex_Series_Vectors_io;
with Standard_Complex_Series_Matrices;
with Standard_Complex_Vector_Series;
with Standard_Complex_Vector_Series_io;  use Standard_Complex_Vector_Series_io;
with Standard_Complex_Matrix_Series;
with Standard_Complex_Matrix_Series_io;  use Standard_Complex_Matrix_Series_io;
with Standard_Random_Series_Vectors;
with Standard_Random_Series_Matrices;
with Standard_Series_Matrix_Solvers;
with Series_Coefficient_Vectors;

procedure ts_mtserlin is

-- DESCRIPTION :
--   Tests the linearization of solving linear systems of truncated series
--   with multitasking.

  procedure Multitasked_Solve_Next_by_lufac
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Applies multitasking for the backsubstitution to solve 
  --   the matrix series equation
  --   defined by the matrix series in A and right hand side in b.

  -- REQUIRED :
  --   A'last >= 0 and b'last >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series,
  --            A(0) contains the LU factorization of A(0);
  --   b        the right hand side as a vector series;
  --   wrk      allocated work space for the nbt tasks.

  -- ON RETURN :
  --   b        all coefficients of the solution series up to b.deg,
  --            provided info = 0.

    use Standard_Series_Matrix_Solvers;

  begin
    null;
  end Multitasked_Solve_Next_by_lufac;

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b.

  -- REQUIRED :
  --   A'last >= 0 and b'last >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   info     if info /= 0, then the lead coefficient matrix of A
  --            was deemed singular and x is undefined,
  --            if info = 0, then the system is regular;
  --   x        all coefficients of the solution series up to b.deg,
  --            provided info = 0.

    use Standard_Series_Matrix_Solvers;

    wrk : Standard_Complex_VecVecs.VecVec(1..nbt);

  begin
    Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0 then
      for k in wrk'range loop
        declare
          cff : Standard_Complex_Vectors.Vector(0..A'last);
        begin
          wrk(k) := new Standard_Complex_Vectors.Vector'(cff);
        end;
      end loop;
      Multitasked_Solve_Next_by_lufac(nbt,A,b,ipvt,wrk);
    end if;
  end Multitasked_Solve_by_lufac;

  procedure Standard_Test ( nbt,n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients in standard double precision.
  --   Converts an n-by-n matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.

    use Standard_Complex_Series_Matrices;
    use Standard_Series_Matrix_Solvers;

    sA : constant Standard_Complex_Series_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant Standard_Complex_Matrix_Series.Matrix 
       := Standard_Complex_Matrix_Series.Create(sA); 
    vm : constant Standard_Complex_VecMats.VecMat(0..As.deg)
       := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
    sx : constant Standard_Complex_Series_Vectors.Vector(1..n)
       := Standard_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    xs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sx);
    sb : constant Standard_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sb);
    sbcff : constant Standard_Complex_VecVecs.VecVec(1..n)
          := Series_Coefficient_Vectors.Standard_Series_Coefficients(sb);
    bscff : constant Standard_Complex_VecVecs.VecVec(0..bs.deg)
          := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..n);
    info : integer32;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The coefficient matrices : "); put(vm);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of b : "); put_line(sbcff);
    put_line("The coefficients of the vector series b :"); put(bs);
    put_line("The coefficients of the vector series b :"); put_line(bscff);
    if nbt > 1 
     then Multitasked_Solve_by_lufac(nbt,vm,bscff,ipvt,info);
     else Solve_by_lufac(vm,bscff,ipvt,info,wrk);
    end if;
    put("info : "); put(info,1); new_line;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The computed leading vector series of the solution :");
    put_line(bscff(0));
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(bscff(k));
    end loop;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system,
  --   the degrees of the series in the system, and the number of tasks.

    nbt,dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the number of equations and variables : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    put("  Give the number of tasks : "); get(nbt);
    new_line;
    Standard_Test(nbt,dim,deg);
  end Main;

begin
  Main;
end ts_mtserlin;
