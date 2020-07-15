package body Standard_Inlined_Linear_Solvers is

  procedure lufac ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    dim : in integer32;
                    ipvt : in out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

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
    rak := rcols(dim); iak := icols(dim);
    smax := abs(rak(dim)) + abs(iak(dim));
    if smax = 0.0
     then info := dim;
    end if;
  end lufac;

  procedure lusolve ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec;
                      icols : in Standard_Floating_VecVecs.Link_to_VecVec;
                      dim : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      rb : in Standard_Floating_Vectors.Link_to_Vector;
                      ib : in Standard_Floating_Vectors.Link_to_Vector ) is

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

  procedure lufco ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    dim : in integer32;
                    ipvt : in out Standard_Integer_Vectors.Vector;
                    rz : in Standard_Floating_Vectors.Link_to_Vector;
                    iz : in Standard_Floating_Vectors.Link_to_Vector;
                    rcond : out double_float ) is

    info : integer32;
    anorm : constant double_float := Norm1(rcols,icols);

  begin
    lufac(rcols,icols,dim,ipvt,info);
    if info /= 0
     then rcond := 0.0;
     else estco(rcols,icols,dim,ipvt,anorm,rz,iz,rcond);
    end if;
  end lufco;

end Standard_Inlined_Linear_Solvers;
