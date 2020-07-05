with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Vector_Splitters;
with Standard_Matrix_Splitters;
with Standard_Inlined_Linear_Solvers;
with Standard_Complex_Vector_Series;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Matrix_Series;
with Standard_Complex_Series_Matrices;
with Standard_Random_Series_Vectors;
with Standard_Random_Series_Matrices;
with Standard_Series_Matrix_Solvers;
with Series_Coefficient_Vectors;

procedure ts_perfserlin is

-- DESCRIPTION :
--   Test the performance on solving linear systems of power series with
--   linearization, flat data structures, and inlined linear solvers.

  procedure Row_Matrix_Multiply
              ( rArows,iArows : in Standard_Floating_VecVecs.Link_to_VecVec;
                rx,ix,ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes the matrix-vector product y = A*x,
  --   with A given as the real and imaginary parts of its rows,
  --   and x and y as vectors of real and imaginary parts.

  -- REQUIRED :
  --   ry'range = iy'range = rArows'range = iArows'range, and 
  --   rx'range = ix'range = rArows(k)'range = iArows(k)'range,
  --   for all k in rArows'range.

  -- ON ENTRY :
  --   rArows   real parts of the complex numbers on the rows of A;
  --   iArows   imaginary parts of the complex numbers on the rows of A;
  --   rx       real parts of the numbers of the complex vector x;
  --   ix       imaginary parts of the numbers of the complex vector x.
  --   ry       space allocated for the real parts of y;
  --   iy       space allocated for the imaginary parts of y.

  -- ON RETURN :
  --   ry       real parts of the complex vector y = A*x;
  --   iy       imaginary parts of the complex vector y = A*x.

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;
    pr,pi,qr,qi : double_float;

  begin
    for k in rArows'range loop
      rlnk := rArows(k); ilnk := iArows(k);
      ry(k) := 0.0; iy(k) := 0.0;
      for j in rx'range loop
        pr := rlnk(j); pi := ilnk(j);
        qr := rx(j);   qi := ix(j);
        ry(k) := ry(k) + pr*qr - pi*qi;
        iy(k) := iy(k) + pr*qi + pi*qr;
      end loop;
    end loop;
  end Row_Matrix_Multiply;

  procedure Split_Rows
              ( A : in Standard_Complex_Matrices.Link_to_Matrix;
                rArows : in Standard_Floating_VecVecs.Link_to_VecVec;
                iArows : in Standard_Floating_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Splits the rows of A in vectors of real and imaginary parts of
  --   the complex numbers on the rows of A.

  -- ON ENTRY :
  --   A        a matrix of complex numbers;
  --   rArows   space allocated for the real parts of all numbers of A;
  --   iArows   space allocated for the imaginary parts of all numbers of A.

  -- ON RETURN :
  --   rArows   rArows(i)(j) contains the real part of A(i,j);
  --   iArows   iArows(i)(j) contains the imaginary part of A(i,j).

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in A'range(1) loop
      rlnk := rArows(i);
      ilnk := iArows(i);
      for j in A'range(2) loop
        rlnk(j) := Standard_Complex_Numbers.REAL_PART(A(i,j));
        ilnk(j) := Standard_Complex_Numbers.IMAG_PART(A(i,j));
      end loop;
    end loop;
  end Split_Rows;

  procedure Split_Rows
              ( vm : in Standard_Complex_VecMats.VecMat;
                rv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec ) is

  -- DESCRIPTION :
  --   Splits the rows of the matrix vm(k) into vectors of real and
  --   imaginary parts of the complex numbers on the rows of vm(k),
  --   for k in rv'range = iv'range.

  -- REQUIRED :
  --   rv'range fits within vm'range and 
  --   rv(k)'range = vm(k)'range(1), for all k in rv'range, and
  --   rv(k)(j)'range = vm(k)'range(2), for all j in rv(k)'range;
  --   iv has the same dimensions as rv.

  -- ON ENTRY :
  --   vm       a vector of complex matrices;
  --   rv       rv(k) has space allocated for the real parts 
  --            of all numbers of vm(k);
  --   iv       iv(k) has space allocated for the imaginary parts 
  --            of all numbers of vm(k).

  -- ON RETURN :
  --   rv       rv(k) stores the real parts of vm(k), for k in rv'range;
  --   iv       iv(k) stores the imaginary parts of vm(k), for k in iv'range.

  begin
    for k in rv'range loop
      Split_Rows(vm(k),rv(k),iv(k));
    end loop;
  end Split_Rows;

  procedure Inlined_Solve_by_lufac
              ( dim : in integer32; 
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, without condition number estimate,
  --   where the matrices A(k) are given as vector 
  --   using an inlined linear system solver,
  --   with allocated real work space vectors.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   dim      dimension of the vectors;
  --   b        the right hand side coefficients of a vector series;
  --   rv       all real pars of all A(k), for k in 1..degree;
  --   iv       all imaginary pars of all A(k), for k in 1..degree;
  --   rc       real parts of the columns of A(0);
  --   ic       imaginary parts of the columns of A(0);
  --   rb       allocated work space vector for all real parts
  --            of the solution vectors;
  --   ib       allocated work space vector for all imaginary parts
  --            of the solution vectors;
  --   ry       allocated work space vector of range 1..dim;
  --   iy       allocated work space vector of range 1..dim.

  -- ON RETURN :
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   b        b contains the coefficients of the solution series x,
  --            provided info = 0, otherwise b is unchanged;
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k);
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    Standard_Inlined_Linear_Solvers.lufac(rc,ic,dim,ipvt,info);
    if info = 0 then
      Standard_Vector_Splitters.Complex_Parts(b(0),rb(0),ib(0));
      Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb(0),ib(0));
      Standard_Vector_Splitters.Complex_Merge(rb(0),ib(0),b(0));
      for k in 1..b'last loop                        -- loop to compute x(k)
        Row_Matrix_Multiply(rv(k),iv(k),rb(0),ib(0),ry,iy); -- y = A(k)*x(0)
        Standard_Vector_Splitters.Complex_Parts(b(k),rb(k),ib(k));
        rlnk := rb(k); ilnk := ib(k);
        for j in rlnk'range loop          -- compute b(k) = b(k) - A(k)*x(0)
          rlnk(j) := rlnk(j) - ry(j);
          ilnk(j) := ilnk(j) - iy(j);
        end loop;
        for i in 1..(k-1) loop
          Row_Matrix_Multiply(rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
          for j in rlnk'range loop        -- substract A(k-1)*x(k) from b(k)
            rlnk(j) := rlnk(j) - ry(j);
            ilnk(j) := ilnk(j) - iy(j);
          end loop;
        end loop;
        Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb(k),ib(k));
        Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),b(k));
      end loop;
    end if;
  end Inlined_Solve_by_lufac;

  procedure Inlined_Solve_by_lufac
              ( A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, without condition number estimate,
  --   using an inlined linear system solver.
  --   Allocates and deallocates all real work space vectors.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series.

  -- ON RETURN :
  --   A        A(0) contains the output of lufac on A(0);
  --   b        b contains the coefficients of the solution series x,
  --            provided info = 0, otherwise b is unchanged;
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

    a0lu : constant Standard_Complex_Matrices.Link_to_Matrix := A(0);
    dim : constant integer32 := a0lu'last(1);
    deg : constant integer32 := A'last;
    rcols,icols,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;

  begin
    Standard_Floating_VecVecVecs.Allocate(rv,1,deg,1,dim,1,dim);
    Standard_Floating_VecVecVecs.Allocate(iv,1,deg,1,dim,1,dim);
    Split_Rows(A,rv,iv);
    rcols := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    icols := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rb := Standard_Vector_Splitters.Allocate(dim,dim,0,1);
    ib := Standard_Vector_Splitters.Allocate(dim,dim,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    Standard_Matrix_Splitters.Complex_Parts(a0lu.all,rcols,icols);
    Inlined_Solve_by_lufac(dim,b,ipvt,info,rv,iv,rcols,icols,rb,ib,ry,iy);
    Standard_Floating_VecVecs.Deep_Clear(rcols);
    Standard_Floating_VecVecs.Deep_Clear(icols);
    Standard_Floating_VecVecs.Deep_Clear(rb);
    Standard_Floating_VecVecs.Deep_Clear(ib);
    Standard_Floating_Vectors.Clear(ry);
    Standard_Floating_Vectors.Clear(iy);
    Standard_Floating_VecVecVecs.Clear(rv);
    Standard_Floating_VecVecVecs.Clear(iv);
  end Inlined_Solve_by_lufac;

  procedure Standard_Test ( n,d : in integer32 ) is

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
    vm1 : constant Standard_Complex_VecMats.VecMat(0..As.deg)
        := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..As.deg)
        := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
    sx : constant Standard_Complex_Series_Vectors.Vector(1..n)
       := Standard_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    xs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sx);
    sb : constant Standard_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sb);
    bscff1 : constant Standard_Complex_VecVecs.VecVec(0..bs.deg)
           := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    bscff2 : constant Standard_Complex_VecVecs.VecVec(0..bs.deg)
           := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..n);

  begin
    Solve_by_lufac(vm1,bscff1,ipvt,info,wrk);
    new_line;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The computed leading vector series of the solution :");
    put_line(bscff1(0));
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(bscff1(k));
    end loop;
    Inlined_Solve_by_lufac(vm2,bscff2,ipvt,info);
    new_line;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The recomputed leading vector series of the solution :");
    put_line(bscff2(0));
    for k in 1..bs.deg loop
      put("The term "); put(k,1); put_line(" in the solution :");
      put_line(xs.cff(k));
      put("Recomputed term "); put(k,1); put_line(" in the solution :");
      put_line(bscff2(k));
    end loop;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degrees of the series in the system.

    dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Performance of solving linear systems of power series ...");
    put("  Give the dimension of the system : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    Standard_Test(dim,deg);
  end Main;

begin
  Main;
end ts_perfserlin;
