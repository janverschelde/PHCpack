with unchecked_deallocation;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
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

  procedure Create ( rx,ix : in Standard_Floating_VecVecs.VecVec;
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

-- REVERSE MODE OF ALGORITHMIC DIFFERENTIATION :

  procedure Speel ( rx,ix : in Standard_Floating_VecVecs.VecVec;
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

  procedure Speel ( rx,ix : in Standard_Floating_VecVecs.VecVec;
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

-- COMPUTING THE POWER TABLE :

  procedure Compute ( rpwt,ipwt : in Link_to_VecVecVec;
                      mxe : in Standard_Integer_Vectors.Vector;
                      rx,ix : in Standard_Floating_VecVecs.VecVec ) is

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

end Standard_Coefficient_Convolutions;
