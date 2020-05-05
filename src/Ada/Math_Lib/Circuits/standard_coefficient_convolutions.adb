with unchecked_deallocation;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Vector_Splitters;

package body Standard_Coefficient_Convolutions is

  function Allocate ( mxe : Standard_Integer_Vectors.Vector;
                      deg : integer32 )
                    return Link_to_VecVecVec is

    res : Link_to_VecVecVec;
    pwt : VecVecVec(mxe'range);

    use Standard_Vector_Splitters;

  begin
    for i in mxe'range loop
      if mxe(i) > 2
       then pwt(i) := Allocate_Floating_Coefficients(mxe(i)-2,deg);
      end if;
    end loop;
    res := new VecVecVec'(pwt);
    return res;
  end Allocate;

  procedure Create ( rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                     mxe : in Standard_Integer_Vectors.Vector;
                     deg : in integer32;
                     rpwt,ipwt : out Link_to_VecVecVec ) is
  begin
    rpwt := Allocate(mxe,deg);
    ipwt := Allocate(mxe,deg);
    Compute(rpwt,ipwt,mxe,rx,ix);
  end Create;

-- BASIC COMPUTATIONAL PROCEDURES :

  procedure Update ( rvl,ivl : in Standard_Floating_Vectors.Link_to_Vector;
                     rnc,inc : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for i in rvl'range loop
      exit when (i > rnc'last);
      rvl(i) := rvl(i) + rnc(i);
      ivl(i) := ivl(i) + inc(i);
    end loop;
  end Update;

  procedure Multiply
              ( xr,xi,yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr,zi : in Standard_Floating_Vectors.Link_to_Vector ) is

    deg : constant integer32 := xr'last;
    rpa,ipa : double_float; -- accumulates real and imaginary parts
    xr0,xi0 : double_float; -- to hold values in xr and xi
    yr0,yi0 : double_float; -- to hold values in yr and yi
    idx : integer32;

  begin
   -- product(0) := first(0)*second(0);
    xr0 := xr(0); xi0 := xi(0); yr0 := yr(0); yi0 := yi(0);
    zr(0) := xr0*yr0 - xi0*yi0; zi(0) := xi0*yr0 + xr0*yi0;
    for k in 1..deg loop
     -- product(k) := first(0)*second(k);
      xr0 := xr(0); xi0 := xi(0); yr0 := yr(k); yi0 := yi(k);
      rpa := xr0*yr0 - xi0*yi0; ipa := xi0*yr0 + xr0*yi0;
      for i in 1..k loop
       -- product(k) := product(k) + first(i)*second(k-i);
        xr0 := xr(i); xi0 := xi(i);
        idx := k-i;
        yr0 := yr(idx); yi0 := yi(idx);
        rpa := rpa + xr0*yr0 - xi0*yi0; ipa := ipa + xi0*yr0 + xr0*yi0;
      end loop;
      zr(k) := rpa; zi(k) := ipa;
    end loop;
  end Multiply;

-- COMPUTING THE POWER TABLE :

  procedure Compute ( rpwt,ipwt : in Link_to_VecVecVec;
                      mxe : in Standard_Integer_Vectors.Vector;
                      rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec ) is

    rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;

  begin
    for i in rx'range loop
      if mxe(i) > 2 then
        rxpw := rpwt(i); ixpw := ipwt(i);
        Multiply(rx(i),ix(i),rx(i),ix(i),rxpw(1),ixpw(1));
        for k in 2..(mxe(i)-2) loop
          Multiply(rxpw(k-1),ixpw(k-1),rx(i),ix(i),rxpw(k),ixpw(k));
        end loop;
      end if;
    end loop;
  end Compute;

-- REVERSE MODE OF ALGORITHMIC DIFFERENTIATION :

  procedure Speel
              ( rx,ix : in Standard_Floating_VecVecs.VecVec;
                rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                rcrs,icrs : in Standard_Floating_VecVecs.VecVec ) is
  begin
    Multiply(rx(1),ix(1),rx(2),ix(2),rfwd(1),ifwd(1));
    for k in 3..rx'last loop
      Multiply(rfwd(k-2),ifwd(k-2),rx(k),ix(k),rfwd(k-1),ifwd(k-1));
    end loop;
    if rx'last > 2 then
      Multiply(rx(rx'last),ix(ix'last),
               rx(rx'last-1),ix(ix'last-1),rbck(1),ibck(1));
      for k in 2..rx'last-2 loop
        Multiply(rbck(k-1),ibck(k-1),rx(rx'last-k),ix(ix'last-k),
                 rbck(k),ibck(k));
      end loop;
      if rx'last = 3 then
        Multiply(rx(1),ix(1),rx(3),ix(3),rcrs(1),icrs(1));
      else
        Multiply(rx(1),ix(1),rbck(rx'last-3),ibck(ix'last-3),
                 rcrs(1),icrs(1));
        for k in 2..rx'last-3 loop
          Multiply(rfwd(k-1),ifwd(k-1),rbck(rx'last-2-k),ibck(ix'last-2-k),
                   rcrs(k),icrs(k));
        end loop;
        Multiply(rfwd(rx'last-3),ifwd(ix'last-3),rx(rx'last),ix(ix'last),
                 rcrs(rx'last-2),icrs(ix'last-2));
      end if;
    end if;
  end Speel;

  procedure Speel
              ( rx,ix : in Standard_Floating_VecVecs.VecVec;
                idx : in Standard_Integer_Vectors.Vector;
                rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                rcrs,icrs : in Standard_Floating_VecVecs.VecVec ) is

    p,q,r : integer32;

  begin
    p := idx(1); q := idx(2);
    Multiply(rx(p),ix(p),rx(q),ix(q),rfwd(1),ifwd(1));
    for k in 3..idx'last loop
      p := k-2; q := idx(k); r := k-1;
      Multiply(rfwd(p),ifwd(p),rx(q),ix(q),rfwd(r),ifwd(r));
    end loop;
    if idx'last > 2 then
      p := idx(idx'last); q := idx(idx'last-1);
      Multiply(rx(p),ix(p),rx(q),ix(q),rbck(1),ibck(1));
      for k in 2..idx'last-2 loop
        p := k-1; q := idx(idx'last-k);
        Multiply(rbck(p),ibck(p),rx(q),ix(q),rbck(k),ibck(k));
      end loop;
      if idx'last = 3 then
        p := idx(1); q := idx(3);
        Multiply(rx(p),ix(p),rx(q),ix(q),rcrs(1),icrs(1));
      else
        p := idx(1); q := idx'last-3;
        Multiply(rx(p),ix(p),rbck(q),ibck(q),rcrs(1),icrs(1));
        for k in 2..idx'last-3 loop
          p := k-1; q := idx'last-2-k;
          Multiply(rfwd(p),ifwd(p),rbck(q),ibck(q),rcrs(k),icrs(k));
        end loop;
        p := idx'last-3; q := idx(idx'last); r := idx'last-2;
        Multiply(rfwd(p),ifwd(p),rx(q),ix(q),rcrs(r),icrs(r));
      end if;
    end if;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec ) is

    use Standard_Integer_Vectors;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    ryptr : constant Standard_Floating_Vectors.Link_to_Vector := ryd(ryd'last);
    iyptr : constant Standard_Floating_Vectors.Link_to_Vector := iyd(iyd'last);
    p : Standard_Floating_Vectors.Link_to_Vector;
    q,r : integer32;

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          q := idk(1);
          Update(ryptr,iyptr,rx(q),ix(q));
          p := ryd(q); p(0) := p(0) + 1.0;
        else
          Speel(rx,ix,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs);
          q := idk'last-1;
          Update(ryptr,iyptr,rfwd(q),ifwd(q));
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Update(ryd(q),iyd(q),rx(r),ix(r));
            Update(ryd(r),iyd(r),rx(q),ix(q));
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Update(ryd(q),iyd(q),rbck(r),ibck(r));
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Update(ryd(q),iyd(q),rcrs(r),icrs(r));
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Update(ryd(q),iyd(q),rfwd(r),ifwd(r));
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    iwrk : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Integer_Vectors;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    ryptr : constant Standard_Floating_Vectors.Link_to_Vector := ryd(ryd'last);
    iyptr : constant Standard_Floating_Vectors.Link_to_Vector := iyd(iyd'last);
    rp,ip : Standard_Floating_Vectors.Link_to_Vector;
    q,r : integer32;

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        rp := rcff(k); ip := icff(k);
        if idk'last = 1 then
          q := idk(1);
          Multiply(rp,ip,rx(q),ix(q),rwrk,iwrk);
          Update(ryptr,iyptr,rwrk,iwrk);
          Update(ryd(q),iyd(q),rp,ip);
        else
          Speel(rx,ix,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs);
          q := idk'last-1;
          Multiply(rp,ip,rfwd(q),ifwd(q),rwrk,iwrk);
          Update(ryptr,iyptr,rwrk,iwrk);
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Multiply(rp,ip,rx(r),ix(r),rwrk,iwrk);
            Update(ryd(q),iyd(q),rwrk,iwrk);
            Multiply(rp,ip,rx(q),ix(q),rwrk,iwrk);
            Update(ryd(r),iyd(r),rwrk,iwrk);
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Multiply(rp,ip,rbck(r),ibck(r),rwrk,iwrk);
            Update(ryd(q),iyd(q),rwrk,iwrk);
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Multiply(rp,ip,rcrs(r),icrs(r),rwrk,iwrk);
              Update(ryd(q),iyd(q),rwrk,iwrk);
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Multiply(rp,ip,rfwd(r),ifwd(r),rwrk,iwrk);
            Update(ryd(q),iyd(q),rwrk,iwrk);
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Multiply_Factor
                  ( xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_Vectors.Link_to_Vector;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec ) is

    rpwx,ipwx : Standard_Floating_VecVecs.Link_to_VecVec;
    rlpw,ilpw : Standard_Floating_Vectors.Link_to_Vector;
    powidx,fptr : integer32;

  begin
    fptr := facidx(facidx'first);
    rpwx := rpwt(fptr);
    ipwx := ipwt(fptr);
    powidx := xpk(fptr);   -- power in power table
    if powidx = 2 then
      Multiply(rcff,icff,rx(fptr),ix(fptr),racc,iacc);
    else
      rlpw := rpwx(powidx-2);  -- coefficients of higher powers
      ilpw := ipwx(powidx-2);
      Multiply(rcff,icff,rlpw,ilpw,racc,iacc);
    end if;
    for k in facidx'first+1..facidx'last loop
      for i in rwrk'range loop
        rwrk(i) := racc(i); iwrk(i) := iacc(i);
      end loop;
      fptr := facidx(k);
      rpwx := rpwt(fptr);
      ipwx := ipwt(fptr);
      powidx := xpk(fptr);   -- power in power table
      if powidx = 2 then
        Multiply(rwrk,iwrk,rx(fptr),ix(fptr),racc,iacc);
      else
        rlpw := rpwx(powidx-2);  -- coefficients of higher powers
        ilpw := ipwx(powidx-2);
        Multiply(rwrk,iwrk,rlpw,ilpw,racc,iacc);
      end if;
    end loop;
  end Multiply_Factor;

  procedure Multiply_Power
                  ( multiplier : in integer32;
                    rcff : in Standard_Floating_Vectors.Link_to_Vector; 
                    icff : in Standard_Floating_Vectors.Link_to_Vector ) is

    factor : constant double_float := Create(multiplier);

  begin
    for i in rcff'range loop
      rcff(i) := factor*rcff(i);
      icff(i) := factor*icff(i);
    end loop;
  end Multiply_Power;

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec ) is

    use Standard_Integer_Vectors;

    idk,xpk,fck : Standard_Integer_Vectors.Link_to_Vector;
    ryptr : constant Standard_Floating_Vectors.Link_to_Vector := ryd(ryd'last);
    iyptr : constant Standard_Floating_Vectors.Link_to_Vector := iyd(iyd'last);
    rpcf,ipcf : Standard_Floating_Vectors.Link_to_Vector;
    pidx,qidx : integer32;

  begin
    for k in idx'range loop
      idk := idx(k);           -- the k-th exponent index 
      if idk /= null then
        xpk := xps(k);         -- the k-th exponent vector
        fck := fac(k);         -- the k-th factor index
        rpcf := rcff(k); ipcf := icff(k);
        if idk'last = 1 then
          pidx := idk(1);
          if fck = null then
            Multiply(rpcf,ipcf,rx(pidx),ix(pidx),rwrk,iwrk);
            Update(ryptr,iyptr,rwrk,iwrk);
            Update(ryd(pidx),iyd(pidx),rpcf,ipcf);
          else
            Multiply_Factor(xpk,fck,rx,ix,rpcf,ipcf,rwrk,iwrk,
                            racc,iacc,rpwt,ipwt);
            Multiply(racc,iacc,rx(pidx),ix(pidx),rwrk,iwrk);
            Update(ryptr,iyptr,rwrk,iwrk);
            Multiply_Power(xpk(pidx),racc,iacc);
            Update(ryd(pidx),iyd(pidx),racc,iacc);
          end if;
        else
          Speel(rx,ix,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs);
          pidx := idk'last-1;
          if fck = null then
            Multiply(rpcf,ipcf,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
          else
            Multiply_Factor(xpk,fck,rx,ix,rpcf,ipcf,rwrk,iwrk,
                            racc,iacc,rpwt,ipwt);
            Multiply(racc,iacc,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
          end if;
          Update(ryptr,iyptr,rwrk,iwrk);
          if idk'last = 2 then
            pidx := idk(1); qidx := idk(2);
            if fck = null then
              Multiply(rpcf,ipcf,rx(pidx),ix(pidx),rwrk,iwrk);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
              Multiply(rpcf,ipcf,rx(qidx),ix(qidx),rwrk,iwrk);
              Update(ryd(pidx),iyd(pidx),rwrk,iwrk);
            else -- use the common factor in acc
              Multiply(racc,iacc,rx(pidx),ix(pidx),rwrk,iwrk);
              if xpk(qidx) > 1
               then Multiply_Power(xpk(qidx),rwrk,iwrk);
              end if;
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
              Multiply(racc,iacc,rx(qidx),ix(qidx),rwrk,iwrk);
              if xpk(pidx) > 1
               then Multiply_Power(xpk(pidx),rwrk,iwrk);
              end if;
              Update(ryd(pidx),iyd(pidx),rwrk,iwrk);
            end if;
          else -- idk'last > 2 
            pidx := idk'last-2; qidx := idk(1);
            if fck = null then
              Multiply(rpcf,ipcf,rbck(pidx),ibck(pidx),rwrk,iwrk);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
              for j in idk'first+1..idk'last-1 loop
                Multiply(rpcf,ipcf,rcrs(j-1),icrs(j-1),rwrk,iwrk);
                qidx := idk(j);
                Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
              end loop;
              Multiply(rpcf,ipcf,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
              qidx := idk(idk'last);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
            else
              Multiply(racc,iacc,rbck(pidx),ibck(pidx),rwrk,iwrk);
              Multiply_Power(xpk(qidx),rwrk,iwrk);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
              for j in idk'first+1..idk'last-1 loop
                Multiply(racc,iacc,rcrs(j-1),icrs(j-1),rwrk,iwrk);
                qidx := idk(j);
                Multiply_Power(xpk(qidx),rwrk,iwrk);
                Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
              end loop;
              Multiply(racc,iacc,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
              qidx := idk(idk'last);
              Multiply_Power(xpk(qidx),rwrk,iwrk);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk);
            end if;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

-- EVALUATION AND DIFFERENTIATION ON CIRCUITS :

  procedure EvalDiff ( c : in Circuit;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec ) is

    use Standard_Floating_Vectors;

  begin
    Speel(c.xps,c.idx,c.fac,c.rcf,c.icf,rx,ix,c.rfwd,c.ifwd,c.rbck,c.ibck,
          c.rcrs,c.icrs,ryd,iyd,c.rwrk,c.iwrk,c.racc,c.iacc,rpwt,ipwt);
    if c.rct /= null and c.ict /= null
     then Update(ryd(ryd'last),iyd(iyd'last),c.rct,c.ict);
    end if;
  end EvalDiff;

-- DEALLOCATORS :

  procedure Clear ( pwt : in out VecVecVec ) is
  begin
    for k in pwt'range loop
      Standard_Floating_VecVecs.Deep_Clear(pwt(k));
    end loop;
  end Clear;

  procedure Clear ( pwt : in out Link_to_VecVecVec ) is

    procedure free is new unchecked_deallocation(VecVecVec,Link_to_VecVecVec);

  begin
    if pwt /= null then
      Clear(pwt.all);
      free(pwt);
    end if;
  end Clear;

  procedure Clear ( pwt : in out VecVecVec_Array ) is
  begin
    for k in pwt'range loop
      Clear(pwt(k));
    end loop;
  end Clear;

  procedure Clear ( c : in out Circuit ) is
  begin
    Standard_Integer_VecVecs.Clear(c.xps);
    Standard_Integer_VecVecs.Clear(c.idx);
    Standard_Integer_VecVecs.Clear(c.fac);
    Standard_Floating_VecVecs.Clear(c.rcf);
    Standard_Floating_VecVecs.Clear(c.icf);
    Standard_Floating_Vectors.Clear(c.rct);
    Standard_Floating_Vectors.Clear(c.ict);
    Standard_Floating_VecVecs.Clear(c.rfwd);
    Standard_Floating_VecVecs.Clear(c.ifwd);
    Standard_Floating_VecVecs.Clear(c.rbck);
    Standard_Floating_VecVecs.Clear(c.ibck);
    Standard_Floating_VecVecs.Clear(c.rcrs);
    Standard_Floating_VecVecs.Clear(c.icrs);
    Standard_Floating_Vectors.Clear(c.rwrk);
    Standard_Floating_Vectors.Clear(c.iwrk);
    Standard_Floating_Vectors.Clear(c.racc);
    Standard_Floating_Vectors.Clear(c.iacc);
  end Clear;

  procedure Clear ( c : in out Link_to_Circuit ) is

    procedure free is
      new unchecked_deallocation(Circuit,Link_to_Circuit);

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

  procedure Clear ( c : in out Link_to_Circuits ) is

    procedure free is
      new unchecked_deallocation(Circuits,Link_to_Circuits);

  begin
    if c /= null then
      Clear(c.all);
      free(c);
    end if;
  end Clear;


end Standard_Coefficient_Convolutions;
