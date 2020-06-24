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
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

procedure ts_perfdlu is

-- DESCRIPTION :
--   Development of a better performing LU factorization of complex matrices.

  procedure lufac ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    dim : in integer32;
                    ipvt : in out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

  -- DESCRIPTION :
  --   Computes the LU factorization of the complex matrix stored as
  --   a pair of real and imaginary parts of the column coefficients.
  --   Applies partial row pivoting.

  -- REQUIRED : rcols'range = icols'range = ipvt'range = 1..dim,
  --   and for k in 1..dim: rcols(k)'range = icols(k)'range.

  -- ON ENTRY :
  --   rcols    real parts of the columns of a complex matrix;
  --   icols    imaginary parts of the columns of a complex matrix;
  --   dim      dimension of the matrix;
  --   ipvt     work space for the pivoting information.

  -- ON RETURN :
  --   rcols    real parts of the LU factorization of the matrix;
  --   icols    imaginary parts of the LU factorization of the matrix;
  --   ipvt     pivoting information;
  --   info     0 if all went well, otherwise info indicates the first
  --            column where a zero pivot occurred.

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

  procedure lusolve ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec;
                      icols : in Standard_Floating_VecVecs.Link_to_VecVec;
                      dim : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      rb : in Standard_Floating_Vectors.Link_to_Vector;
                      ib : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Uses the output of lufac to solve a linear system.

  -- REQUIRED : rcols'range = icols'range = ipvt'range = 1..dim,
  --   and for k in 1..dim: rcols(k)'range = icols(k)'range;
  --   and rb'range = ib'range = 1..dim.

  -- ON ENTRY :
  --   rcols    real parts of the LU factorization of the matrix;
  --   icols    imaginary parts of the LU factorization of the matrix;
  --   dim      dimension of the matrix;
  --   ipvt     pivoting information, as computed by lufac;
  --   rb       real parts of the right hand side vector;
  --   ib       imaginary parts of the right hand side vector.

  -- ON RETURN :
  --   rb       real parts of the solution vector;
  --   ib       imaginary parts of the solution vector.

    ell,nm1,kb : integer32;
    rtmp,itmp,pr,pi,zr,zi : double_float;
    rak,iak : Standard_Floating_Vectors.Link_to_Vector;
 
  begin
    nm1 := dim-1;
    if nm1 >= 1 then                                       -- solve L*y = b
      for k in 1..nm1 loop
        ell := ipvt(k);
        rtmp := rb(ell);
        itmp := ib(ell);
        if ell /= k then
          rb(ell) := rb(k); rb(k) := rtmp; 
          ib(ell) := ib(k); ib(k) := itmp;
        end if;
        rak := rcols(k);
        iak := icols(k);
        for i in (k+1)..dim loop                -- b(i) := b(i) + tmp*ak(i)
          pr := rak(i); pi := iak(i);
          zr := rtmp*pr - itmp*pi;
          zi := rtmp*pi + itmp*pr;
          rb(i) := rb(i) + zr;
          ib(i) := ib(i) + zi;
        end loop;
      end loop;
    end if;
    for k in 1..dim loop                                   -- solve U*x = y
      kb := dim+1-k;
      rak := rcols(kb);
      iak := icols(kb);
      pr := rak(kb); pi := iak(kb);                -- b(kb) := b(kb)/ak(kb)
      rtmp := pr*pr + pi*pi; 
      pr := pr/rtmp; pi := -pi/rtmp;            -- p = pr + i*pi = 1/ak(kb)
      zr := rb(kb); zi := ib(kb);
      rb(kb) := zr*pr - zi*pi;                     -- multiply b(kb) with p
      ib(kb) := zr*pi + zi*pr;
      rtmp := -rb(kb);
      itmp := -ib(kb);
      for j in 1..(kb-1) loop                   -- b(j) := b(j) + tmp*ak(j)
        pr := rak(j); pi := iak(j);
        zr := rtmp*pr - itmp*pi;
        zi := rtmp*pi + itmp*pr;
        rb(j) := rb(j) + zr;
        ib(j) := ib(j) + zi;
      end loop;
    end loop;
  end lusolve;

  function Norm1 ( rcols : Standard_Floating_VecVecs.Link_to_VecVec;
                   icols : Standard_Floating_VecVecs.Link_to_VecVec )
                 return double_float is

  -- DESCRIPTION :
  --   Returns the 1-norm of the complex matrix with real parts of 
  --   its columns in rcols and the imaginary parts in icols.

  -- REQUIRED : rcols'range = icols'range, and for all k in rcols'range:
  --   rcols(k)'range = icols(k)'range.

    res : double_float := 0.0;
    sum : double_float;
    raj,iaj : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for j in rcols'range loop
      sum := 0.0;
      raj := rcols(j);
      iaj := icols(j);
      for i in raj'range loop
        sum := sum + abs(raj(i)) + abs(iaj(i));
      end loop;
      if sum > res
       then res := sum;
      end if; 
    end loop;
    return res;
  end Norm1;

  procedure estco ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    dim : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in double_float;
                    rz : in Standard_Floating_Vectors.Link_to_Vector;
                    iz : in Standard_Floating_Vectors.Link_to_Vector;
                    rcond : out double_float ) is

  -- DESCRIPTION :
  --   Estimates the condition number, given a LU factorization.

  -- REQUIRED : rcols'range = icols'range = ipvt'range = 1..dim,
  --   and for k in 1..dim: rcols(k)'range = icols(k)'range;
  --   and rz'range = iz'range = 1..dim.

  -- ON ENTRY :
  --   rcols    columns with the real parts of the outcome of lufac;
  --   icols    columns with the imaginary parts of the outcome of lufac;
  --   dim      dimension of the matrix;
  --   ipvt     pivoting information, as computed by lufac;
  --   anorm    the 1-norm of the matrix, computed by Norm1,
  --            on the original matrix before the lufac;
  --   rz       work space vector of 1..dim for real parts;
  --   iz       work space vector of 1..dim for imaginary parts.

  -- ON RETURN :
  --   rcond    estimated for the inverse of the condition number.

    rak,iak,raj,iaj : Standard_Floating_Vectors.Link_to_Vector;
    kb,kp1,ell : integer32;
    s,sm,sum,ynorm,pr,pi,zr,zi : double_float;
    rek,iek,rt,it,rwk,iwk,rwkm,iwkm : double_float;

  begin
    rek := 1.0; iek := 0.0;                          -- solve ctrans(u)*w = e
    for j in 1..dim loop
      rz(j) := 0.0; iz(j) := 0.0;
    end loop;
    for k in 1..dim loop
      rak := rcols(k); iak := icols(k);
      sum := abs(rz(k)) + abs(iz(k));
      if sum /= 0.0 then                             -- ek := csign(ek,-z(k))
        s := abs(rek) + abs(iek);
        rek := -s*rz(k)/sum; iek := -s*iz(k)/sum;
      end if;
      sum := abs(rek - rz(k)) + abs(iek - iz(k));
      sm := abs(rak(k)) + abs(iak(k));
     -- if cabs(ek-z(k)) > cabs(ak(k)) then
      if sum > sm then
        s := sm/sum;                        -- s := cabs(ak(k))/cabs(ek-z(k))
        for i in 1..dim loop  
          rz(i) := s*rz(i);                            -- z := Create(s) * z
          iz(i) := s*iz(i);
        end loop;
        rek := s*rek; iek := s*iek;                  -- ek := Create(s) * ek
      end if;
      rwk := rek - rz(k);                               -- wk := ek - z(k)
      iwk := iek - iz(k); 
      rwkm := -rek - rz(k);                             -- wkm := -ek - z(k)
      iwkm := -iek - iz(k);
      s := abs(rwk) + abs(iwk);                         -- s := cabs(wk)
      sm := abs(rwkm) + abs(iwkm);                      -- sm := cabs(wkm)
      sum := abs(rak(k)) + abs(iak(k));
      if sum = 0.0  then
        rwk := 1.0; iwk := 0.0;                         -- wk := Create(1.0)
        rwkm := 1.0; iwkm := 0.0;                       -- wkm := Create(1.0)
      else
        pr := rak(k); pi := -iak(k);       -- p = pr + i*pi = dconjg(ak(k))
        zr := pr*pr + pi*pi;
        pr := pr/zr; pi := -pi/zr;         -- p = pr + i*pi = 1/dconjg(ak(k))
        zr := rwk*pr - iwk*pi;
        zi := rwk*pi + iwk*pr;
        rwk := zr; iwk := zi;              -- wk := wk / dconjg(ak(k));
        zr := rwkm*pr - iwkm*pi;
        zi := rwkm*pi + iwkm*pr;
        rwkm := zr; iwkm := zi;            -- wkm := wkm / dconjg(ak(k));
      end if;
      kp1 := k + 1;
      if kp1 <= dim then
        for j in kp1..dim loop
          raj := rcols(j); iaj := icols(j);
          pr := raj(k); pi := -iaj(k);      -- p = pr + i*pi = dconjg(aj(k))
          zr := rwkm*pr - iwkm*pi;
          zi := rwkm*pi + iwkm*pr;      -- z = zr * i*zi = wkm*dconjg(aj(k))
          sm := sm + abs(rz(j) + zr) + abs(iz(j) + zi);
                                -- sm := sm + cabs(z(j) + wkm*dconjg(aj(k)))
          pr := raj(k); pi := -iaj(k);      -- p = pr + i*pi = dconjg(aj(k))
          zr := rwk*pr - iwk*pi;
          zi := rwk*pi + iwk*pr;         -- z = zr + i*zi = wk*dconjg(aj(k))
          rz(j) := rz(j) + zr;
          iz(j) := iz(j) + zi;           -- z(j) := z(j) + wk*dconjg(aj(k));
          s := s + abs(rz(j)) + abs(iz(j));           -- s := s + cabs(z(j))
        end loop;
        if s < sm then
          rt := rwkm - rwk;
          it := iwkm - iwk;
          rwk := rwkm;
          iwk := iwkm;
          for j in kp1..dim loop
            raj := rcols(j); iaj := icols(j);
            pr := raj(k); pi := -iaj(k);    -- p = pr + i*pi = dconjg(aj(k))
            zr := rt*pr - it*pi;
            zi := rt*pi + it*pr;          -- z = zr + i*zi = t*dconjg(aj(k))
            rz(j) := rz(j) + zr;
            iz(j) := iz(j) + zi;           -- z(j) := z(j) + t*dconjg(aj(k))
          end loop;
        end if;
      end if;
      rz(k) := rwk; iz(k) := iwk;
    end loop;
    sum := 0.0;
    for i in 1..dim loop
      sum := sum + abs(rz(i)) + abs(iz(i));     --  sum := sum + cabs(z(i))
    end loop;
    s := 1.0/sum;
    for i in 1..dim loop
      rz(i) := s*rz(i);
      iz(i) := s*iz(i);                               -- z := Create(s) * z
    end loop;
    for k in 1..dim loop                           -- solve ctrans(L)*y = w
      kb := dim+1-k;
      if kb < dim then
        rak := rcols(kb);
        iak := icols(kb);
        rt := 0.0;
        it := 0.0;
        for i in (kb+1)..dim loop
          pr := rak(i); pi := -iak(i);     -- p = pr + i*pi = dconjg(ak(i))
          zr := rz(i); zi := iz(i);
          rt := rt + pr*zr - pi*zi;
          it := it + pr*zi + pi*zr;          -- t := t + dconjg(ak(i))*z(i)
        end loop;
        rz(kb) := rz(kb) + rt;
        iz(kb) := iz(kb) + it;
      end if;
      sum := abs(rz(kb)) + abs(iz(kb));
      if sum > 1.0 then
        s := 1.0/sum;
        for i in 1..dim loop
          rz(i) := s*rz(i);
          iz(i) := s*iz(i);                           -- z := Create(s) * z
        end loop;
      end if;
      ell := ipvt(kb);
      rt := rz(ell); rz(ell) := rz(kb); rz(kb) := rt;
      it := iz(ell); iz(ell) := iz(kb); iz(kb) := it;
    end loop;
    sum := 0.0;
    for i in 1..dim loop
      sum := sum + abs(rz(i)) + abs(iz(i));
    end loop;
    s := 1.0/sum;
    for i in 1..dim loop
      rz(i) := s*rz(i);
      iz(i) := s*iz(i);                              --  z := Create(s) * z
    end loop;
    ynorm := 1.0;
    for k in 1..dim loop                                   -- solve L*v = y
      ell := ipvt(k);
      rt := rz(ell); rz(ell) := rz(k); rz(k) := rt;
      it := iz(ell); iz(ell) := iz(k); iz(k) := it;
      if k < dim then
        rak := rcols(k);
        iak := icols(k);
        for i in (k+1)..dim loop
          pr := rak(i); pi := iak(i);
          zr := rt*pr - it*pi;
          zi := rt*pi + it*pr;
          rz(i) := rz(i) + zr;
          iz(i) := iz(i) + zi;                  -- z(i) := z(i) + t * ak(i)
        end loop;
      end if;
      sum := abs(rz(k)) + abs(iz(k));
      if sum > 1.0 then
        s := 1.0/sum;
        for i in 1..dim loop
          rz(i) := s*rz(i);
          iz(i) := s*iz(i);                          --  z := Create(s) * z
        end loop;
        ynorm := s*ynorm;
      end if;
    end loop;
    sum := 0.0;
    for i in 1..dim loop
      sum := sum + abs(rz(i)) + abs(iz(i));
    end loop;
    s := 1.0/sum;
    for i in 1..dim loop
      rz(i) := s*rz(i);
      iz(i) := s*iz(i);                             --  z := Create(s) * z
    end loop;
    ynorm := s*ynorm;
    for k in 1..dim loop                                  -- solve u*z = v
      kb := dim+1-k;
      rak := rcols(kb);
      iak := icols(kb);
      sum := abs(rz(kb)) + abs(iz(kb));
      sm := abs(rak(kb)) + abs(iak(kb));
      if sum > sm then                    -- if cabs(z(kb)) > cabs(ak(kb))
        s := sm/sum;                    -- s := cabs(ak(kb)) / cabs(z(kb))
        for i in 1..dim loop
          rz(i) := s*rz(i);
          iz(i) := s*iz(i);                         --  z := Create(s) * z
        end loop;
        ynorm := s*ynorm;
      end if;
      sum := abs(rak(kb)) + abs(iak(kb));
      if sum = 0.0 then
        rz(kb) := 1.0; iz(kb) := 0.0;
      else
        pr := rak(kb); pi := iak(kb);
        zr := pr*pr + pi*pi;
        pr := pr/zr; pi := -pi/zr;             -- p = pr + i*pi = 1/ak(kb)
        zr := rz(kb); zi := iz(kb);
        rz(kb) := zr*pr - zi*pi;
        iz(kb) := zr*pi + zi*pr;               -- z(kb) := z(kb) / ak(kb)
      end if;
      rt := -rz(kb);
      it := -iz(kb);
      for i in 1..(kb-1) loop
        pr := rak(i); pi := iak(i);
        zr := rt*pr - it*pi;
        zi := rt*pi + it*pr;
        rz(i) := rz(i) + zr;
        iz(i) := iz(i) + zi;                   -- z(i) := z(i) + t * ak(i)
      end loop;
    end loop;
    sum := 0.0;                                        -- make znorm = 1.0
    for i in 1..dim loop
      sum := sum + abs(rz(i)) + abs(iz(i));
    end loop;
    s := 1.0/sum;
    ynorm := s*ynorm;
    if anorm = 0.0
     then rcond := 0.0;
     else rcond := ynorm/anorm;
    end if;
  end estco;

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
    lufac(wrk,dim,ipvt,info);
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
    lufac(rvv,ivv,dim,ipvt,info);
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
    lufac(wrk,dim,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    lusolve(wrk,dim,ipvt,wrkrhs);
    put_line("The solution :"); put_line(wrkrhs);
    Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
    lufac(rvv,ivv,dim,ipvt,info);
    put("info after the vectorized lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    Standard_Vector_Splitters.Complex_Parts(rhs,rb,ib);
    lusolve(rvv,ivv,dim,ipvt,rb,ib);
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
    lufac(wrk,dim,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    nrm1 := Norm1(mat);
    put("1-norm of the matrix :"); put(nrm1); new_line;
    estco(wrk,dim,ipvt,nrm1,rco);
    put("estimated inverse condition :"); put(rco); new_line;
    Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
    nrm1 := Norm1(rvv,ivv);
    put("   recomputed 1-norm :"); put(nrm1); new_line;
    lufac(rvv,ivv,dim,ipvt,info);
    put("info after the vectorized lufac : "); put(info,1); new_line;
    put("ipvt : "); put(ipvt); new_line;
    estco(rvv,ivv,dim,ipvt,nrm1,rz,iz,rco);
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
      lufac(wrk,dim,ipvt,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex LU factorization");
    tstart(timer);
    for k in 1..frq loop
      Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
      lufac(rvv,ivv,dim,ipvt,info);
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
      lufac(wrk,dim,ipvt,info);
      wrkrhs := rhs;
      lusolve(wrk,dim,ipvt,wrkrhs);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex lufac + lusolve");
    tstart(timer);
    for k in 1..frq loop
      Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
      lufac(rvv,ivv,dim,ipvt,info);
      Standard_Vector_Splitters.Complex_Parts(rhs,rb,ib);
      lusolve(rvv,ivv,dim,ipvt,rb,ib);
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
      nrm := Norm1(mat);
      lufac(wrk,dim,ipvt,info);
      estco(wrk,dim,ipvt,nrm,rco);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex lufac + estco");
    tstart(timer);
    for k in 1..frq loop
      Standard_Matrix_Splitters.Complex_Parts(mat,rvv,ivv);
      nrm := Norm1(rvv,ivv);
      lufac(rvv,ivv,dim,ipvt,info);
      estco(rvv,ivv,dim,ipvt,nrm,rz,iz,rco);
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
