with unchecked_deallocation;

package body Standard_Coefficient_Circuits is

  function Allocate ( nbr,dim : integer32 ) return Circuit is

    res : Circuit(nbr);
    rforward : constant Standard_Floating_Vectors.Vector(1..dim-1)
             := (1..dim-1 => 0.0);
    iforward : constant Standard_Floating_Vectors.Vector(1..dim-1)
             := (1..dim-1 => 0.0);
    rbackward : constant Standard_Floating_Vectors.Vector(1..dim-2)
              := (1..dim-2 => 0.0);
    ibackward : constant Standard_Floating_Vectors.Vector(1..dim-2)
              := (1..dim-2 => 0.0);
    rcross : constant Standard_Floating_Vectors.Vector(1..dim-2)
           := (1..dim-2 => 0.0);
    icross : constant Standard_Floating_Vectors.Vector(1..dim-2)
           := (1..dim-2 => 0.0);

  begin
    res.dim := dim;
    res.rfwd := new Standard_Floating_Vectors.Vector'(rforward);
    res.ifwd := new Standard_Floating_Vectors.Vector'(iforward);
    res.rbck := new Standard_Floating_Vectors.Vector'(rbackward);
    res.ibck := new Standard_Floating_Vectors.Vector'(ibackward);
    res.rcrs := new Standard_Floating_Vectors.Vector'(rcross);
    res.icrs := new Standard_Floating_Vectors.Vector'(icross);
    return res;
  end Allocate;

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUIT :

  procedure Speel ( c : in Circuit;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Speel(c.xps,c.rcf,c.icf,c.rcst,c.icst,xr,xi,ryd,iyd,
          c.rfwd,c.ifwd,c.rbck,c.ibck,c.rcrs,c.icrs);
  end Speel;

  procedure Speel ( c : in Circuit;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec ) is
  begin
    Speel(c.xps,c.idx,c.fac,c.rcf,c.icf,c.rcst,c.icst,xr,xi,ryd,iyd,
          c.rfwd,c.ifwd,c.rbck,c.ibck,c.rcrs,c.icrs,rpwt,ipwt);
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rcf : in Standard_Floating_Vectors.Vector;
                    icf : in Standard_Floating_Vectors.Vector;
                    rcst,icst : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector ) is

    idk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    rkcff,ikcff : double_float;
    zr,zi,pr,pi : double_float;

    use Standard_Integer_Vectors;

  begin
    ryd(0) := rcst; iyd(0) := icst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          idx1 := idk(1); rkcff := rcf(k); ikcff := icf(k);
         -- yd(0) := yd(0) + cff*x(idx1);
          pr := xr(idx1); pi := xr(idx1);
          zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
          ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
         -- yd(idx1) := yd(idx1) + cff;
          ryd(idx1) := ryd(idx1) + rkcff; iyd(idx1) := iyd(idx1) + rkcff;
        else
          Fused_Forward_Backward_Cross
            (idk.all,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
         -- cff := c.cff(k); yd(0) := yd(0) + cff*c.fwd(idk'last-1); 
          rkcff := rcf(k); ikcff := icf(k);
          idx1 := idk'last-1; pr := rfwd(idx1); pi := ifwd(idx1);
          zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
          ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
          if idk'last = 2 then
            idx1 := idk(1); idx2 := idk(2);
           -- yd(idx2) := yd(idx2) + cff*x(idx1);
            pr := xr(idx1); pi := xi(idx1);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx2) := ryd(idx2) + zr; iyd(idx2) := iyd(idx2) + zi;
           -- yd(idx1) := yd(idx1) + cff*x(idx2);
            pr := xr(idx2); pi := xi(idx2);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
          else -- idk'last > 2
            idx1 := idk(1); 
           -- yd(idx1) := yd(idx1) + cff*c.bck(idk'last-2);
            idx2 := idk'last-2; pr := rbck(idx2); pi := ibck(idx2);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
            for j in idk'first+1..idk'last-1 loop
              idx1 := idk(j);
             -- yd(idx1) := yd(idx1) + cff*c.crs(j-1);
              idx2 := j-1; pr := rcrs(idx2); pi := icrs(idx2);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
            end loop;
            idx1 := idk(idk'last);
           -- yd(idx1) := yd(idx1) + cff*c.fwd(idk'last-2);
            idx2 := idk'last-2; pr := rfwd(idx2); pi := ifwd(idx2);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcf : in Standard_Floating_Vectors.Vector;
                    icf : in Standard_Floating_Vectors.Vector;
                    rcst,icst : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec ) is

    idk,fck,xpk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    rkcff,ikcff,racc,iacc,factor : double_float;
    zr,zi,pr,pi : double_float;

    use Standard_Integer_Vectors;

  begin
    ryd(0) := rcst; iyd(0) := icst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        idx1 := idk(1); rkcff := rcf(k); ikcff := icf(k);
        fck := fac(k);
        if fck = null then
          if idk'last = 1 then
           -- yd(0) := yd(0) + cff*x(idx1);
            pr := xr(idx1); pi := xr(idx1);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
           -- yd(idx1) := yd(idx1) + cff;
            ryd(idx1) := ryd(idx1) + rkcff; iyd(idx1) := iyd(idx1) + rkcff;
          else
            Fused_Forward_Backward_Cross
              (idk.all,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
           -- cff := c.cff(k); yd(0) := yd(0) + cff*c.fwd(idk'last-1); 
            rkcff := rcf(k); ikcff := icf(k);
            idx1 := idk'last-1; pr := rfwd(idx1); pi := ifwd(idx1);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
            if idk'last = 2 then
              idx1 := idk(1); idx2 := idk(2);
             -- yd(idx2) := yd(idx2) + cff*x(idx1);
              pr := xr(idx1); pi := xi(idx1);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(idx2) := ryd(idx2) + zr; iyd(idx2) := iyd(idx2) + zi;
             -- yd(idx1) := yd(idx1) + cff*x(idx2);
              pr := xr(idx2); pi := xi(idx2);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
            else -- idk'last > 2
              idx1 := idk(1); 
             -- yd(idx1) := yd(idx1) + cff*c.bck(idk'last-2);
              idx2 := idk'last-2; pr := rbck(idx2); pi := ibck(idx2);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
              for j in idk'first+1..idk'last-1 loop
                idx1 := idk(j);
               -- yd(idx1) := yd(idx1) + cff*c.crs(j-1);
                idx2 := j-1; pr := rcrs(idx2); pi := icrs(idx2);
                zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
                ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
              end loop;
              idx1 := idk(idk'last);
             -- yd(idx1) := yd(idx1) + cff*c.fwd(idk'last-2);
              idx2 := idk'last-2; pr := rfwd(idx2); pi := ifwd(idx2);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
            end if;
          end if;
        else
          xpk := xps(k);
          if idk'last = 1 then
            Multiply_Factor(xpk,fck,xr,xi,rkcff,ikcff,rpwt,ipwt,racc,iacc);
           -- yd(0) := yd(0) + acc*x(idx1);
            pr := xr(idx1); pi := xr(idx1);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
           -- yd(idx1) := yd(idx1) + factor*acc;
            factor := Create(xpk(idx1));
            ryd(idx1) := ryd(idx1) + factor*racc;
            iyd(idx1) := iyd(idx1) + factor*iacc;
          else
            Fused_Forward_Backward_Cross
              (idk.all,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
            Multiply_Factor(xpk,fck,xr,xi,rkcff,ikcff,rpwt,ipwt,racc,iacc);
            pr := rfwd(idk'last-1); pi := ifwd(idk'last-1);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
            if idk'last = 2 then
              idx2 := idk(2);
             -- wrk := acc*x(idx1); factor := Create(xpk(idx2));
             -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
              pr := xr(idx1); pi := xi(idx1);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              factor := Create(xpk(idx2));
              ryd(idx2) := ryd(idx2) + factor*zr;
              iyd(idx2) := iyd(idx2) + factor*zi;
             -- wrk := acc*x(idx2); factor := Create(xpk(idx1));
             -- wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
              pr := xr(idx2); pi := xi(idx2);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              factor := Create(xpk(idx1));
              ryd(idx1) := ryd(idx1) + factor*zr;
              iyd(idx1) := iyd(idx1) + factor*zi;
            else -- idk'last > 2
             -- wrk := acc*bck(idk'last-2);
              pr := rbck(idk'last-2); pi := ibck(idk'last-2);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              factor := Create(xpk(idx1));
             -- wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
              ryd(idx1) := ryd(idx1) + factor*zr;
              iyd(idx1) := iyd(idx1) + factor*zi;
              for j in idk'first+1..idk'last-1 loop
               -- wrk := acc*crs(j-1);
                pr := rcrs(j-1); pi := icrs(j-1);
                zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
                idx2 := idk(j); factor := Create(xpk(idx2));
               -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
                ryd(idx2) := ryd(idx2) + factor*zr;
                iyd(idx2) := iyd(idx2) + factor*zi;
              end loop;
             -- wrk := acc*fwd(idk'last-2);
              idx2 := idk(idk'last);
              pr := rfwd(idk'last-2); pi := ifwd(idk'last-2);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              factor := Create(xpk(idx2));
             -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
              ryd(idx2) := ryd(idx2) + factor*zr;
              iyd(idx2) := iyd(idx2) + factor*zi;
            end if;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

-- AUXILIARY PROCEDURES :

  procedure Forward ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                      xi : in Standard_Floating_Vectors.Link_to_Vector;
                      fr : in Standard_Floating_Vectors.Link_to_Vector;
                      fi : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
  end Forward;

  procedure Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
   -- b(b'first) := x(x'last)*x(x'last-1);
    pr := xr(xr'last); pi := xi(xr'last);
    idx := xi'last-1;
    qr := xr(idx); qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    br(1) := zr; bi(1) := zi;
    for k in 2..xr'last-2 loop
     -- b(k) := b(k-1)*x(x'last-k);
      pr := zr; pi := zi;
      idx := xr'last-k;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      br(k) := zr; bi(k) := zi;
    end loop;
  end Forward_Backward;

  procedure Fused_Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr1,zi1,pr1,pi1,qr1,qi1 : double_float;
    zr2,zi2,pr2,pi2,qr2,qi2 : double_float;
    idx1,idx2 : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr1 := xr(1); pi1 := xi(1);
    idx1 := xr'first+1;
    qr1 := xr(idx1);  qi1 := xi(idx1);
    zr1 := pr1*qr1 - pi1*qi1;
    zi1 := pr1*qi1 + pi1*qr1;
    fr(1) := zr1; fi(1) := zi1;
   -- b(b'first) := x(x'last)*x(x'last-1);
    pr2 := xr(xr'last); pi2 := xi(xr'last);
    idx2 := xi'last-1;
    qr2 := xr(idx2); qi2 := xi(idx2);
    zr2 := pr2*qr2 - pi2*qi2;
    zi2 := pr2*qi2 + pi2*qr2;
    br(1) := zr2; bi(1) := zi2;
    for k in 2..dim-1 loop 
     -- f(k) := f(k-1)*x(k+1);
      pr1 := zr1; pi1 := zi1;
      idx1 := k+1;
      qr1 := xr(idx1); qi1 := xi(idx1);
      zr1 := pr1*qr1 - pi1*qi1;
      zi1 := pr1*qi1 + pi1*qr1;
      fr(k) := zr1; fi(k) := zi1;
     -- b(k) := b(k-1)*x(x'last-k);
      pr2 := zr2; pi2 := zi2;
      idx2 := xr'last-k;
      qr2 := xr(idx2); qi2 := xi(idx2);
      zr2 := pr2*qr2 - pi2*qi2;
      zi2 := pr2*qi2 + pi2*qr2;
      br(k) := zr2; bi(k) := zi2;
    end loop;
    if dim > 1 then
      pr1 := zr1; pi1 := zi1;
      idx1 := dim+1;
      qr1 := xr(idx1); qi1 := xi(idx1);
      zr1 := pr1*qr1 - pi1*qi1;
      zi1 := pr1*qi1 + pi1*qr1;
      fr(dim) := zr1; fi(dim) := zi1;
    end if;
  end Fused_Forward_Backward;

  procedure Forward_Backward_Cross
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
   -- b(b'first) := x(x'last)*x(x'last-1);
    pr := xr(xr'last); pi := xi(xr'last);
    idx := xi'last-1;
    qr := xr(idx); qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    br(1) := zr; bi(1) := zi;
    for k in 2..xr'last-2 loop
     -- b(k) := b(k-1)*x(x'last-k);
      pr := zr; pi := zi;
      idx := xr'last-k;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      br(k) := zr; bi(k) := zi;
    end loop;
    if xr'last > 2 then
      if xr'last = 3 then
       -- c(1) := x(1)*x(3)
        pr := xr(1); pi := xi(1);
        qr := xr(3); qi := xi(3);
        zr := pr*qr - pi*qi;
        zi := pr*qi + pi*qr;
        cr(1) := zr; ci(1) := zi;
      else
       -- c(1) := x(1)*b(x'last-3);
        pr := xr(1); pi := xi(1);
        idx := xr'last-3;
        qr := br(idx); qi := bi(idx);
        zr := pr*qr - pi*qi;
        zi := pr*qi + pi*qr;
        cr(1) := zr; ci(1) := zi;
        for k in 2..xr'last-3 loop
         -- c(k) := f(k-1)*b(x'last-2-k);
          idx := k-1;
          pr := fr(idx); pi := fi(idx);
          idx := xr'last-2-k;
          qr := br(idx); qi := bi(idx);
          zr := pr*qr - pi*qi;
          zi := pr*qi + pi*qr;
          cr(k) := zr; ci(k) := zi;
        end loop;
        pr := xr(xr'last); pi := xi(xi'last);
        idx := xr'last-3;
        qr := fr(idx); qi := fi(idx);
        zr := pr*qr - pi*qi;
        zi := pr*qi + pi*qr;
        idx := xr'last-2;
        cr(idx) := zr; ci(idx) := zi;
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Fused_Forward_Backward_Cross
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr1,zi1,pr1,pi1,qr1,qi1 : double_float;
    zr2,zi2,pr2,pi2,qr2,qi2 : double_float;
    zr3,zi3,pr3,pi3,qr3,qi3 : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;
    firstend,lastend,plusidx,minidx : integer32;

  begin
    if xr'last >= 8 then
      if xr'last mod 2 = 0 then
        lastend := xr'last-4;
        firstend := lastend/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        pr1 := xr(1); pi1 := xi(1);
        idx := xr'first+1;
        qr1 := xr(idx);  qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        pr2 := xr(xr'last); pi2 := xi(xr'last);
        idx := xi'last-1;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        minidx := firstend+1; plusidx := minidx+1;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          idx := plusidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-plusidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
          plusidx := plusidx + 1;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          idx := minidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-minidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
          minidx := minidx - 1;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        idx := dim+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        pr3 := xr(1); pi3 := xi(1);
        idx := xr'last-3;
        qr3 := br(idx); qi3 := bi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        pr3 := xr(xr'last); pi3 := xi(xi'last);
        idx := xr'last-3;
        qr3 := fr(idx); qi3 := fi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        idx := xr'last-2;
        cr(idx) := zr3; ci(idx) := zi3;
      else
        lastend := xr'last-4;
        firstend := (xr'last-3)/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        pr1 := xr(1); pi1 := xi(1);
        idx := xr'first+1;
        qr1 := xr(idx);  qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        pr2 := xr(xr'last); pi2 := xi(xr'last);
        idx := xi'last-1;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        plusidx := firstend+1;
       -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
       -- idx := plusidx-1;
       -- pr3 := fr(idx); pi3 := fi(idx);
       -- f(plusidx-1) is still in zr1 and zi1
        idx := xr'last-2-plusidx;
        qr3 := br(idx); qi3 := bi(idx);
       -- zr3 := pr3*qr3 - pi3*qi3;
       -- zi3 := pr3*qi3 + pi3*qr3;
        zr3 := zr1*qr3 - zi1*qi3;
        zi3 := zr1*qi3 + zi1*qr3;
        cr(plusidx) := zr3; ci(plusidx) := zi3;
        minidx := plusidx;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          plusidx := plusidx + 1;
          idx := plusidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-plusidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          minidx := minidx - 1;
          idx := minidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-minidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        idx := dim+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        pr3 := xr(1); pi3 := xi(1);
        idx := xr'last-3;
        qr3 := br(idx); qi3 := bi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        pr3 := xr(xr'last); pi3 := xi(xi'last);
        idx := xr'last-3;
        qr3 := fr(idx); qi3 := fi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        idx := xr'last-2;
        cr(idx) := zr3; ci(idx) := zi3;
      end if;
    else
     -- f(f'first) := x(x'first)*x(x'first+1);
      pr1 := xr(1); pi1 := xi(1);
      idx := xr'first+1;
      qr1 := xr(idx);  qi1 := xi(idx);
      zr1 := pr1*qr1 - pi1*qi1;
      zi1 := pr1*qi1 + pi1*qr1;
      fr(1) := zr1; fi(1) := zi1;
     -- b(b'first) := x(x'last)*x(x'last-1);
      pr2 := xr(xr'last); pi2 := xi(xr'last);
      idx := xi'last-1;
      qr2 := xr(idx); qi2 := xi(idx);
      zr2 := pr2*qr2 - pi2*qi2;
      zi2 := pr2*qi2 + pi2*qr2;
      br(1) := zr2; bi(1) := zi2;
      for k in 2..xr'last-2 loop 
       -- f(k) := f(k-1)*x(k+1);
        pr1 := zr1; pi1 := zi1;
        idx := k+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(k) := zr1; fi(k) := zi1;
       -- b(k) := b(k-1)*x(x'last-k);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-k;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(k) := zr2; bi(k) := zi2;
      end loop;
      if dim > 1 then
        pr1 := zr1; pi1 := zi1;
        idx := dim+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
      end if;
      if xr'last > 2 then
        if xr'last = 3 then
         -- c(1) := x(1)*x(3)
          pr3 := xr(1); pi3 := xi(1);
          qr3 := xr(3); qi3 := xi(3);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(1) := zr3; ci(1) := zi3;
        else
         -- c(1) := x(1)*b(x'last-3);
          pr3 := xr(1); pi3 := xi(1);
          idx := xr'last-3;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(1) := zr3; ci(1) := zi3;
          for k in 2..xr'last-3 loop
           -- c(k) := f(k-1)*b(x'last-2-k);
            idx := k-1;
            pr3 := fr(idx); pi3 := fi(idx);
            idx := xr'last-2-k;
            qr3 := br(idx); qi3 := bi(idx);
            zr3 := pr3*qr3 - pi3*qi3;
            zi3 := pr3*qi3 + pi3*qr3;
            cr(k) := zr3; ci(k) := zi3;
          end loop;
          pr3 := xr(xr'last); pi3 := xi(xi'last);
          idx := xr'last-3;
          qr3 := fr(idx); qi3 := fi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          idx := xr'last-2;
          cr(idx) := zr3; ci(idx) := zi3;
        end if;
      end if;
    end if;
  end Fused_Forward_Backward_Cross;

  procedure Fused_Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr1,zi1,pr1,pi1,qr1,qi1 : double_float;
    zr2,zi2,pr2,pi2,qr2,qi2 : double_float;
    zr3,zi3,pr3,pi3,qr3,qi3 : double_float;
    ptr : integer32;
    dim : constant integer32 := idx'last-1;
    firstend,lastend,plusidx,minidx : integer32;

  begin
    if idx'last >= 8 then
      if idx'last mod 2 = 0 then
        lastend := idx'last-4;
        firstend := lastend/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        ptr := idx(1); pr1 := xr(ptr); pi1 := xi(ptr);
        ptr := idx(idx'first+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        ptr := idx(idx'last); pr2 := xr(ptr); pi2 := xi(ptr);
        ptr := idx(idx'last-1); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        minidx := firstend+1; plusidx := minidx+1;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          ptr := plusidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-plusidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
          plusidx := plusidx+1;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          ptr := minidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-minidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
          minidx := minidx - 1;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(dim+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := br(ptr); qi3 := bi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        ptr := idx(idx'last); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := fr(ptr); qi3 := fi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        ptr := idx'last-2; cr(ptr) := zr3; ci(ptr) := zi3;
      else
        lastend := idx'last-4;
        firstend := (idx'last-3)/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        ptr := idx(1); pr1 := xr(ptr); pi1 := xi(ptr);
        ptr := idx(idx'first+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        ptr := idx(idx'last); pr2 := xr(ptr); pi2 := xi(ptr);
        ptr := idx(idx'last-1); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        plusidx := firstend+1;
       -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
        ptr := plusidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
        ptr := idx'last-2-plusidx; qr3 := br(ptr); qi3 := bi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        cr(plusidx) := zr3; ci(plusidx) := zi3;
        minidx := plusidx;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          plusidx := plusidx + 1;
          ptr := plusidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-plusidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          minidx := minidx - 1;
          ptr := minidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-minidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(dim+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := br(ptr); qi3 := bi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        ptr := idx(idx'last); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := fr(ptr); qi3 := fi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        ptr := idx'last-2;
        cr(ptr) := zr3; ci(ptr) := zi3;
      end if;
    else
     -- f(f'first) := x(x'first)*x(x'first+1);
      ptr := idx(1); pr1 := xr(ptr); pi1 := xi(ptr);
      ptr := idx(idx'first+1); qr1 := xr(ptr); qi1 := xi(ptr);
      zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
      fr(1) := zr1; fi(1) := zi1;
      if idx'last > 2 then
       -- b(b'first) := x(x'last)*x(x'last-1);
        ptr := idx(idx'last); pr2 := xr(ptr); pi2 := xi(ptr);
        ptr := idx(idx'last-1); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..idx'last-2 loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        if dim > 1 then
          pr1 := zr1; pi1 := zi1;
          ptr := idx(dim+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(dim) := zr1; fi(dim) := zi1;
        end if;
        if idx'last = 3 then
         -- c(1) := x(1)*x(3)
          ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
          ptr := idx(3); qr3 := xr(ptr); qi3 := xi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(1) := zr3; ci(1) := zi3;
        else
         -- c(1) := x(1)*b(x'last-3);
          ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
          ptr := idx'last-3; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
         -- put("Assigning "); put(zr3); put_line(" to cr(1)");
         -- put("Assigning "); put(zi3); put_line(" to ci(1)");
          cr(1) := zr3; ci(1) := zi3;
          for k in 2..idx'last-3 loop
           -- c(k) := f(k-1)*b(x'last-2-k);
            ptr := k-1; pr3 := fr(ptr); pi3 := fi(ptr);
            ptr := idx'last-2-k; qr3 := br(ptr); qi3 := bi(ptr);
            zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
            cr(k) := zr3; ci(k) := zi3;
          end loop;
          ptr := idx(idx'last); pr3 := xr(ptr); pi3 := xi(ptr);
          ptr := idx'last-3; qr3 := fr(ptr); qi3 := fi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          ptr := idx'last-2; cr(ptr) := zr3; ci(ptr) := zi3;
        end if;
      end if;
    end if;
  end Fused_Forward_Backward_Cross;

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(mxe'range);

  begin
    for k in mxe'range loop
      if mxe(k) > 1 then
        res(k) := new Standard_Floating_Vectors.Vector'(1..mxe(k)-1 => 0.0);
      end if;
    end loop;
    return res;
  end Allocate;

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec ) is

    rlnk : Standard_Floating_Vectors.Link_to_Vector;
    ilnk : Standard_Floating_Vectors.Link_to_Vector;
    zr,zi,yr,yi,xrk,xik : double_float;

  begin
    for k in xr'range loop
      if mxe(k) > 1 then
        rlnk := rpwt(k); ilnk := ipwt(k);
       -- lnk(1) := x(k)*x(k);
        xrk := xr(k); xik := xi(k);
        zr := xrk*xrk - xik*xik;
        zi := 2.0*xrk*xik;
        rlnk(1) := zr; ilnk(1) := zi;
        for i in 2..mxe(k)-1 loop
         -- lnk(i) := lnk(i-1)*x(k);
          yr := zr; yi := zi;
          zr := xrk*yr - xik*yi;
          zi := xrk*yi + xik*yr;
          rlnk(i) := zr; ilnk(i) := zi;
        end loop;
      end if;
    end loop;
  end Power_Table;

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                rcf,icf : in double_float;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                rpf,ipf : out double_float ) is

    rpwx : Standard_Floating_Vectors.Link_to_Vector;
    ipwx : Standard_Floating_Vectors.Link_to_Vector;
    idx,powidx : integer32;
    zr,zi,xrk,xik : double_float;

  begin
    idx := fac(fac'first); powidx := xps(idx);
    if powidx = 2 then
     -- res := cff*x(idx);
      xrk := xr(idx); xik := xi(idx);
      zr := xrk*rcf - xik*icf;
      zi := xrk*icf + xik*rcf;
      rpf := zr; ipf := zi;
    else
      rpwx := rpwt(idx); ipwx := ipwt(idx);
     -- res := cff*pwx(powidx-2);
      xrk := rpwx(powidx-2); xik := ipwx(powidx-2);
      zr := xrk*rcf - xik*icf;
      zi := xrk*icf + xik*rcf;
      rpf := zr; ipf := zi;
    end if;
    for k in fac'first+1..fac'last loop
      idx := fac(k); powidx := xps(idx);
      if powidx = 2 then
       -- res := res*x(idx);
        xrk := xr(idx); xik := xi(idx);
        zr := xrk*rpf - xik*ipf;
        zi := xrk*ipf + xik*rpf;
        rpf := zr; ipf := zi;
      else
        rpwx := rpwt(idx); ipwx := ipwt(idx);
       -- res := res*pwx(powidx-2);
        xrk := rpwx(powidx-2); xik := ipwx(powidx-2);
        zr := xrk*rpf - xik*ipf;
        zi := xrk*ipf + xik*rpf;
        rpf := zr; ipf := zi;
      end if;
    end loop;
  end Multiply_Factor;

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit ) is
  begin
    Standard_Floating_Vectors.Clear(c.rfwd);
    Standard_Floating_Vectors.Clear(c.ifwd);
    Standard_Floating_Vectors.Clear(c.rbck);
    Standard_Floating_Vectors.Clear(c.ibck);
    Standard_Floating_Vectors.Clear(c.rcrs);
    Standard_Floating_Vectors.Clear(c.icrs);
  end Clear;

  procedure Clear ( c : in out Link_to_Circuit ) is

    procedure free is new unchecked_deallocation(Circuit,Link_to_Circuit);

  begin
    if c /= null then
      Clear(c.all);
      free(c);
    end if;
  end Clear;

  procedure Clear ( c : in out Circuits ) is
  begin
    for k in c'range loop
      Clear(c(k));
    end loop;
  end Clear;

end Standard_Coefficient_Circuits;
