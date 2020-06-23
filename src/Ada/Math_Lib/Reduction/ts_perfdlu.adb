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
with Standard_Vector_Splitters;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

procedure ts_perfdlu is

-- DESCRIPTION :
--   Development of a better performing LU factorization of complex matrices.

  procedure Complex_Parts
              ( mat : in Standard_Complex_Matrices.Matrix;
                rvv,ivv : in Standard_Floating_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Returns in rvv the real parts of the elements in mat
  --   and in ivv the imaginary parts of the elements in mat.
  --   The vector representation of the matrix mat is column wise.

  -- REQUIRED :
  --   rvv'range = ivv'range = mat'range(2), and for all k in rvv'range:
  --   rvv(k)'range = ivv(k)'range = mat'range(1).

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in mat'range(2) loop -- split k-th column of the matrix
      rlnk := rvv(k);
      ilnk := ivv(k);
      for i in mat'range(1) loop
        rlnk(i) := Standard_Complex_Numbers.REAL_PART(mat(i,k));
        ilnk(i) := Standard_Complex_Numbers.IMAG_PART(mat(i,k));
      end loop;
    end loop;
  end Complex_Parts;

  procedure Complex_Merge
              ( rvv,ivv : in Standard_Floating_VecVecs.Link_to_VecVec;
                mat : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in rvv and ivv are the real and imaginary parts of the
  --   columns of a complex matrix, with values defined on return.
  --   Mirrors the Complex_Parts procedure.

  -- REQUIRED :
  --   rvv'range = ivv'range = mat'range(2), and for all k in rvv'range:
  --   rvv(k)'range = ivv(k)'range = mat'range(1).

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in rvv'range loop
      rlnk := rvv(k);
      ilnk := ivv(k);
      for i in rlnk'range loop
        mat(i,k) := Standard_Complex_Numbers.Create(rlnk(i),ilnk(i));
      end loop;
    end loop;
  end Complex_Merge;

  procedure lufac ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    dim : in integer32;
                    ipvt : in out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

  -- DESCRIPTION :
  --   Computes the LU factorization of the complex matrix stored as
  --   a pair of real and imaginary parts of the column coefficients.
  --   Applies partial row pivoting.

    kp1,ell,nm1 : integer32;
    smax,tmp,pr,pi,zr,zi,qr,qi : double_float;
    rak,iak,raj,iaj : Standard_Floating_Vectors.Link_to_Vector;

  begin
    info := 0;
    nm1 := integer32(dim - 1);
    if nm1 >= 1 then
      for k in 1..nm1 loop
        rak := rcols(k); iak := icols(k);
        kp1 := k + 1;
        ell := k;                               -- find the pivot index ell
        smax := abs(rak(k)) + abs(iak(k));
        for i in kp1..dim loop
          tmp := abs(rak(i)) + abs(iak(i));
          if tmp > smax 
           then ell := i; smax := tmp;
          end if;
        end loop;
        ipvt(k) := ell;
        if smax = 0.0 then
          info := k;               -- this column is already triangulated
        else
          if ell /= k then                    -- interchange if necessary
            tmp := rak(k); rak(k) := rak(ell); rak(ell) := tmp;
            tmp := iak(k); iak(k) := iak(ell); iak(ell) := tmp; 
          end if;
          pr := rak(k); pi := iak(k);              -- compute multipliers
          tmp := pr*pr + pi*pi;                -- square of norm of ak(k)
          pr := -pr/tmp; pi := pi/tmp;           -- denote acc = -1/ak(k)
          for i in kp1..dim loop
            zr := rak(i); zi := iak(i);             -- ak(i) := ak(i)*acc
            rak(i) := zr*pr - zi*pi;
            iak(i) := zr*pi + zi*pr;
          end loop;
          for j in kp1..dim loop                       -- row elimination
            raj := rcols(j); iaj := icols(j);
            if ell /= k then -- Swap(aj(ell),aj(k));
              tmp := raj(k); raj(k) := raj(ell); raj(ell) := tmp;
              tmp := iaj(k); iaj(k) := iaj(ell); iaj(ell) := tmp;
            end if;
            for i in kp1..dim loop       --  aj(i) := aj(i) + aj(k)*ak(i)
              pr := raj(k); pi := iaj(k);
              qr := rak(i); qi := iak(i);
              zr := pr*qr - pi*qi;
              zi := pr*qi + pi*qr;
              raj(i) := raj(i) + zr;
              iaj(i) := iaj(i) + zi;
            end loop;
          end loop;
        end if;
      end loop;
    end if;
    ipvt(dim) := dim;
    if rak(dim) = 0.0 and iak(dim) = 0.0
     then info := dim;
    end if;
  end lufac;

  procedure Test ( dim : in integer32 ) is

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
    lufac(wrk,dim,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Complex_Parts(mat,rvv,ivv);
    for k in mat'range(2) loop
      put("Elements in column "); put(k,1); put_line(" :");
      for i in mat'range(1) loop
        put(mat(i,k)); new_line;
        put(rvv(k)(i)); put("  "); put(ivv(k)(i)); new_line;
      end loop;
    end loop;
    lufac(rvv,ivv,dim,ipvt,info);
    put("info after the vectorized lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Complex_Merge(rvv,ivv,vvfac);
    put_line("The complex LU factorization :"); put(wrk,3);
    put_line("The recomputed LU factorization :"); put(vvfac,3);
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Test;

  procedure Timed_Test ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the correctness of a LU factorization on a splitted matrix,
  --   randomly generated for dimension dim,
  --   and times the runs for as many times as the value of frq.

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
      lufac(wrk,dim,ipvt,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex LU factorization");
    tstart(timer);
    for k in 1..frq loop
      Complex_Parts(mat,rvv,ivv);
      lufac(rvv,ivv,dim,ipvt,info);
      Complex_Merge(rvv,ivv,vvfac);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"vectorized LU factorization");
    Standard_Floating_VecVecs.Deep_Clear(rvv);
    Standard_Floating_VecVecs.Deep_Clear(ivv);
  end Timed_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension,
  --   generates a random matrix and runs some tests.

    dim,frq : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Interactive test ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test(dim);
    else
      put("Give the frequency : "); get(frq);
      Timed_Test(dim,frq);
    end if;
  end Main;

begin
  Main;
end ts_perfdlu;
