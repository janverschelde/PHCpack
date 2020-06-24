with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Vector_Splitters;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;
with Standard_Matrix_Splitters;
with Standard_Complex_Linear_Solvers;
with Standard_Inlined_Linear_Solvers;

procedure ts_perfdlu is

-- DESCRIPTION :
--   Development of a better performing LU factorization of complex matrices.

  procedure Test_lufac ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the correctness of a LU factorization on a splitted matrix,
  --   randomly generated for dimension dim.

    ndm : constant natural32 := natural32(dim);
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Standard_Random_Matrices.Random_Matrix(ndm,ndm);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    rvv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ivv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    vvfac : Standard_Complex_Matrices.Matrix(1..dim,1..dim);

  begin
    put("A random matrix of dimension "); put(dim,1); put_line(" :");
    put(mat,3);
    Standard_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
    for k in mat'range(2) loop
      put("Elements in column "); put(k,1); put_line(" :");
      for i in mat'range(1) loop
        put(mat(i,k)); new_line;
        put(rvv(k)(i)); put("  "); put(ivv(k)(i)); new_line;
      end loop;
    end loop;
    Standard_Inlined_Linear_Solvers.lufac(rvv,ivv,dim,ipvt,info);
    put("info after the vectorized lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Standard_Matrix_Splitters.Complex_Merge(rvv,ivv,vvfac);
    put_line("The complex LU factorization :"); put(wrk,3);
    put_line("The recomputed LU factorization :"); put(vvfac,3);
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Test_lufac;

  procedure Test_lusolve ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the correctness of solving a linear system on a splitted matrix,
  --   randomly generated for dimension dim.
  --   The right hand side of the linear system is computed so the exact
  --   solution to the linear system is a vector of ones.

    use Standard_Complex_Matrices; -- for computation of rhs

    ndm : constant natural32 := natural32(dim);
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Standard_Random_Matrices.Random_Matrix(ndm,ndm);
    sol : constant Standard_Complex_Vectors.Vector(1..dim)
        := (1..dim => Standard_Complex_Numbers.Create(1.0));
    rhs : constant Standard_Complex_Vectors.Vector(1..dim) := mat*sol;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    wrkrhs : Standard_Complex_Vectors.Vector(1..dim) := rhs;
    rvv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ivv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rb : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    ib : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);

  begin
    put("A random matrix of dimension "); put(dim,1); put_line(" :");
    put(mat,3);
    Standard_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Standard_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,wrkrhs);
    put_line("The solution :"); put_line(wrkrhs);
    Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
    Standard_Inlined_Linear_Solvers.lufac(rvv,ivv,dim,ipvt,info);
    put("info after the vectorized lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Standard_Vector_Splitters.Complex_Parts(rhs,rb,ib);
    Standard_Inlined_Linear_Solvers.lusolve(rvv,ivv,dim,ipvt,rb,ib);
    Standard_Vector_Splitters.Complex_Merge(rb,ib,wrkrhs);
    put_line("The recomputed solution :"); put_line(wrkrhs);
    Standard_Floating_Vectors.Clear(rb);
    Standard_Floating_Vectors.Clear(ib);
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Test_lusolve;

  function Diagonal ( dim : integer32; first : double_float )
                    return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a diagonal matrix with ones on its diagonal,
  --   except for the first element which equals first,
  --   to generate random matrices with a controlled condition number.

    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        if i = j
         then res(i,j) := Standard_Complex_Numbers.Create(1.0);
         else res(i,j) := Standard_Complex_Numbers.Create(0.0);
        end if;
      end loop;
    end loop;
    res(1,1) := Standard_Complex_Numbers.Create(first);
    return res;
  end Diagonal;

  procedure Test_estco ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the correctness of estimating the condition number
  --   of a splitted matrix, randomly generated for dimension dim.

    use Standard_Complex_Matrices; -- for computation of the matrix

    ndm : constant natural32 := natural32(dim);
    mt1 : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Standard_Random_Matrices.Random_Matrix(ndm,ndm);
    mt2 : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Diagonal(dim,1.0E-6);
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := mt1*mt2;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    rvv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ivv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rz : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    iz : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    nrm1,rco : double_float;

  begin
    put("A random matrix of dimension "); put(dim,1); put_line(" :");
    put(mat,3);
    Standard_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    nrm1 := Standard_Complex_Linear_Solvers.Norm1(mat);
    put("1-norm of the matrix :"); put(nrm1); new_line;
    Standard_Complex_Linear_Solvers.estco(wrk,dim,ipvt,nrm1,rco);
    put("estimated inverse condition :"); put(rco); new_line;
    Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
    nrm1 := Standard_Inlined_Linear_Solvers.Norm1(rvv,ivv);
    put("   recomputed 1-norm :"); put(nrm1); new_line;
    Standard_Inlined_Linear_Solvers.lufac(rvv,ivv,dim,ipvt,info);
    put("info after the vectorized lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Standard_Inlined_Linear_Solvers.estco(rvv,ivv,dim,ipvt,nrm1,rz,iz,rco);
    put("       recomputed condition :"); put(rco); new_line;
    Standard_Floating_Vectors.Clear(rz);
    Standard_Floating_Vectors.Clear(iz);
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Test_estco;

  procedure Time_lufac ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   For randomly generated matrices of dimension dim,
  --   does as many LU factorizations as the value of frq,
  --   on complex matrices and splitted matrices.

    timer : Timing_Widget;
    ndm : constant natural32 := natural32(dim);
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Standard_Random_Matrices.Random_Matrix(ndm,ndm);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    rvv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ivv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    vvfac : Standard_Complex_Matrices.Matrix(1..dim,1..dim);

  begin
    tstart(timer);
    for k in 1..frq loop
      wrk := mat;
      Standard_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex LU factorization");
    tstart(timer);
    for k in 1..frq loop
      Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
      Standard_Inlined_Linear_Solvers.lufac(rvv,ivv,dim,ipvt,info);
      Standard_Matrix_Splitters.Complex_Merge(rvv,ivv,vvfac);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"vectorized LU factorization");
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Time_lufac;

  procedure Time_lusolve ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   For randomly generated matrices of dimension dim,
  --   does as many calls to lufac and lusolve as the value of frq,
  --   on complex matrices and splitted matrices.

    use Standard_Complex_Matrices;

    timer : Timing_Widget;
    ndm : constant natural32 := natural32(dim);
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Standard_Random_Matrices.Random_Matrix(ndm,ndm);
    sol : constant Standard_Complex_Vectors.Vector(1..dim)
        := (1..dim => Standard_Complex_Numbers.Create(1.0));
    rhs : constant Standard_Complex_Vectors.Vector(1..dim) := mat*sol;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    wrkrhs : Standard_Complex_Vectors.Vector(1..dim);
    rvv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ivv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rb : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    ib : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);

  begin
    tstart(timer);
    for k in 1..frq loop
      wrk := mat;
      Standard_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
      wrkrhs := rhs;
      Standard_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,wrkrhs);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex lufac + lusolve");
    tstart(timer);
    for k in 1..frq loop
      Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
      Standard_Inlined_Linear_Solvers.lufac(rvv,ivv,dim,ipvt,info);
      Standard_Vector_Splitters.Complex_Parts(rhs,rb,ib);
      Standard_Inlined_Linear_Solvers.lusolve(rvv,ivv,dim,ipvt,rb,ib);
      Standard_Vector_Splitters.Complex_Merge(rb,ib,wrkrhs);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"vectorized lufac + lusolve");
    Standard_Floating_Vectors.Clear(rb);
    Standard_Floating_Vectors.Clear(ib);
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Time_lusolve;

  procedure Time_estco ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   For randomly generated matrices of dimension dim,
  --   does as many calls to lufac and lusolve as the value of frq,
  --   on complex matrices and splitted matrices.

    timer : Timing_Widget;
    ndm : constant natural32 := natural32(dim);
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Standard_Random_Matrices.Random_Matrix(ndm,ndm);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    rvv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ivv : Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rz : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    iz : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    nrm,rco : double_float;

  begin
    tstart(timer);
    for k in 1..frq loop
      wrk := mat;
      nrm := Standard_Complex_Linear_Solvers.Norm1(mat);
      Standard_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
      Standard_Complex_Linear_Solvers.estco(wrk,dim,ipvt,nrm,rco);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex lufac + estco");
    tstart(timer);
    for k in 1..frq loop
      Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
      nrm := Standard_Inlined_Linear_Solvers.Norm1(rvv,ivv);
      Standard_Inlined_Linear_Solvers.lufac(rvv,ivv,dim,ipvt,info);
      Standard_Inlined_Linear_Solvers.estco(rvv,ivv,dim,ipvt,nrm,rz,iz,rco);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"vectorized lufac + estco");
    Standard_Floating_Vectors.Clear(rz);
    Standard_Floating_Vectors.Clear(iz);
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Time_estco;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension,
  --   generates a random matrix and runs some tests.

    dim,frq : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    new_line;
    put_line("MENU to test LU on splitted matrices :");
    put_line("  1. test lufac on random matrices");
    put_line("  2. test lusolve on random matrices");
    put_line("  3. test estco on random matrices");
    put("Type 1, 2, or 3 to select a test : ");
    Ask_Alternative(ans,"13");
    new_line;
    put("Give the frequency (0 for interactive) : "); get(frq);
    case ans is
      when '1' => if frq = 0
                   then Test_lufac(dim);
                   else Time_lufac(dim,frq);
                  end if;
      when '2' => if frq = 0
                   then Test_lusolve(dim);
                   else Time_lusolve(dim,frq);
                  end if;
      when '3' => if frq = 0
                   then Test_estco(dim);
                   else Time_estco(dim,frq);
                  end if;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_perfdlu;
