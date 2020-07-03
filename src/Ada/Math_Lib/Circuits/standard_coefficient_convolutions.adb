with unchecked_deallocation;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Vector_Splitters;
with Exponent_Indices;

package body Standard_Coefficient_Convolutions is

-- ALLOCATORS AND CONSTRUCTORS :

  function Exponent_Maxima
             ( c : Circuits; dim : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Indices.Maxima(c(c'first).xps);

  begin
    for k in c'first+1..c'last loop
      declare
        mxe : constant Standard_Integer_Vectors.Vector(1..dim)
            := Exponent_Indices.Maxima(c(k).xps);
      begin
        for i in mxe'range loop
          if mxe(i) > res(i)
           then res(i) := mxe(i);
          end if;
        end loop;
      end;
    end loop;
    return res;
  end Exponent_Maxima;

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

  function Linearized_Allocation
             ( dim,deg : integer32 )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(0..deg);

  begin
    for k in 0..deg loop
      declare
        cff : constant Standard_Complex_Vectors.Vector(1..dim)
            := (1..dim => Standard_Complex_Numbers.Create(0.0));
      begin
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Linearized_Allocation;

  function Allocate_Coefficients
             ( nbq,nvr,deg : integer32 )
             return Standard_Complex_VecMats.VecMat is

    res : Standard_Complex_VecMats.VecMat(0..deg);

  begin
    for k in res'range loop
      declare
        mat : Standard_Complex_Matrices.Matrix(1..nbq,1..nvr);
      begin
        for i in 1..nbq loop
          for j in 1..nvr loop
            mat(i,j) := Standard_Complex_Numbers.Create(0.0);
          end loop;
        end loop;
        res(k) := new Standard_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Allocate_Coefficients;

  function Create ( c : Circuits; dim,deg : integer32 ) return System is

    neq : constant integer32 := c'last;
    res : System(neq,neq+1,dim,dim+1,deg);

    use Standard_Vector_Splitters;

  begin
    res.crc := c;
    res.mxe := Exponent_Maxima(c,dim);
    res.rpwt := Allocate(res.mxe,deg);
    res.ipwt := Allocate(res.mxe,deg);
    res.ryd := Allocate_Floating_Coefficients(dim+1,deg);
    res.iyd := Allocate_Floating_Coefficients(dim+1,deg);
    res.vy := Linearized_Allocation(neq,deg);
    res.yv := Allocate_Complex_Coefficients(neq,deg);
    res.vm := Allocate_Coefficients(neq,dim,deg);
    return res;
  end Create;

  function Create ( c : Circuits;
                    dim,deg : integer32 ) return Link_to_System is

    res_rep : constant System(c'last,c'last+1,dim,dim+1,deg)
            := Create(c,dim,deg);
    res : constant Link_to_System := new System'(res_rep);

  begin
    return res;
  end Create;

-- COPY PROCEDURES :

  procedure Copy ( c_from : in Circuit; c_to : out Circuit ) is

    use Standard_Floating_Vectors;

  begin
    Standard_Integer_VecVecs.Copy(c_from.xps,c_to.xps);
    Standard_Integer_VecVecs.Copy(c_from.idx,c_to.idx);
    Standard_Integer_VecVecs.Copy(c_from.fac,c_to.fac);
    Standard_Floating_VecVecs.Copy(c_from.rcf,c_to.rcf);
    Standard_Floating_VecVecs.Copy(c_from.icf,c_to.icf);
    if c_from.rct /= null
     then c_to.rct := new Standard_Floating_Vectors.Vector'(c_from.rct.all);
    end if;
    if c_from.ict /= null
     then c_to.ict := new Standard_Floating_Vectors.Vector'(c_from.ict.all);
    end if;
    Standard_Floating_VecVecs.Copy(c_from.rfwd,c_to.rfwd);
    Standard_Floating_VecVecs.Copy(c_from.ifwd,c_to.ifwd);
    Standard_Floating_VecVecs.Copy(c_from.rbck,c_to.rbck);
    Standard_Floating_VecVecs.Copy(c_from.ibck,c_to.ibck);
    Standard_Floating_VecVecs.Copy(c_from.rcrs,c_to.rcrs);
    Standard_Floating_VecVecs.Copy(c_from.icrs,c_to.icrs);
    if c_from.rwrk /= null
     then c_to.rwrk := new Standard_Floating_Vectors.Vector'(c_from.rwrk.all);
    end if;
    if c_from.iwrk /= null
     then c_to.iwrk := new Standard_Floating_Vectors.Vector'(c_from.iwrk.all);
    end if;
    if c_from.racc /= null
     then c_to.racc := new Standard_Floating_Vectors.Vector'(c_from.racc.all);
    end if;
    if c_from.iacc /= null
     then c_to.iacc := new Standard_Floating_Vectors.Vector'(c_from.iacc.all);
    end if;
  end Copy;

  procedure Copy ( c_from : in Link_to_Circuit; c_to : out Link_to_Circuit ) is
  begin
    Clear(c_to); -- if c_from = null, then c_to becomes null as well
    if c_from /= null then
      declare
        crc : Circuit(c_from.nbr,c_from.dim,c_from.dim1,c_from.dim2);
      begin
        Copy(c_from.all,crc);
        c_to := new Circuit'(crc);
      end;
    end if;
  end Copy;

  procedure Copy ( c_from : in Circuits; c_to : out Circuits ) is
  begin
    for k in c_from'range loop
      Copy(c_from(k),c_to(k));
    end loop;
  end Copy;

  procedure Copy ( s_from : in System; s_to : out System ) is
  begin
    Copy(s_from.crc,s_to.crc);
    s_to.mxe := s_from.mxe;
    Copy(s_from.rpwt,s_to.rpwt);
    Copy(s_from.ipwt,s_to.ipwt);
    Standard_Floating_VecVecs.Copy(s_from.ryd,s_to.ryd);
    Standard_Floating_VecVecs.Copy(s_from.iyd,s_to.iyd);
    Standard_Complex_VecVecs.Copy(s_from.vy,s_to.vy);
    Standard_Complex_VecVecs.Copy(s_from.yv,s_to.yv);
    Standard_Complex_VecMats.Copy(s_from.vm,s_to.vm);
  end Copy;

  procedure Copy ( s_from : in Link_to_System; s_to : out Link_to_System ) is
  begin
    Clear(s_to);
    if s_from /= null then
      declare
        s : System(s_from.neq,s_from.neq1,s_from.dim,s_from.dim1,s_from.deg);
      begin
        Copy(s_from.all,s);
        s_to := new System'(s);
      end;
    end if;
  end Copy;

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

  procedure Update ( deg : in integer32;
                     rvl,ivl : in Standard_Floating_Vectors.Link_to_Vector;
                     rnc,inc : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for i in rvl'first..deg loop
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

  procedure Multiply
              ( deg : in integer32;
                xr,xi,yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr,zi : in Standard_Floating_Vectors.Link_to_Vector ) is

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

  procedure Compute ( deg : in integer32; rpwt,ipwt : in Link_to_VecVecVec;
                      mxe : in Standard_Integer_Vectors.Vector;
                      rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec ) is

    rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;

  begin
    for i in rx'range loop
      if mxe(i) > 2 then
        rxpw := rpwt(i); ixpw := ipwt(i);
        Multiply(deg,rx(i),ix(i),rx(i),ix(i),rxpw(1),ixpw(1));
        for k in 2..(mxe(i)-2) loop
          Multiply(deg,rxpw(k-1),ixpw(k-1),rx(i),ix(i),rxpw(k),ixpw(k));
        end loop;
      end if;
    end loop;
  end Compute;

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

  procedure Speel ( deg : in integer32;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec ) is
  begin
    Multiply(deg,rx(1),ix(1),rx(2),ix(2),rfwd(1),ifwd(1));
    for k in 3..rx'last loop
      Multiply(deg,rfwd(k-2),ifwd(k-2),rx(k),ix(k),rfwd(k-1),ifwd(k-1));
    end loop;
    if rx'last > 2 then
      Multiply(deg,rx(rx'last),ix(ix'last),
               rx(rx'last-1),ix(ix'last-1),rbck(1),ibck(1));
      for k in 2..rx'last-2 loop
        Multiply(deg,rbck(k-1),ibck(k-1),rx(rx'last-k),ix(ix'last-k),
                 rbck(k),ibck(k));
      end loop;
      if rx'last = 3 then
        Multiply(deg,rx(1),ix(1),rx(3),ix(3),rcrs(1),icrs(1));
      else
        Multiply(deg,rx(1),ix(1),rbck(rx'last-3),ibck(ix'last-3),
                 rcrs(1),icrs(1));
        for k in 2..rx'last-3 loop
          Multiply(deg,rfwd(k-1),ifwd(k-1),rbck(rx'last-2-k),
                   ibck(ix'last-2-k),rcrs(k),icrs(k));
        end loop;
        Multiply(deg,rfwd(rx'last-3),ifwd(ix'last-3),rx(rx'last),ix(ix'last),
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

  procedure Speel ( deg : in integer32;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec ) is

    p,q,r : integer32;

  begin
    p := idx(1); q := idx(2);
    Multiply(deg,rx(p),ix(p),rx(q),ix(q),rfwd(1),ifwd(1));
    for k in 3..idx'last loop
      p := k-2; q := idx(k); r := k-1;
      Multiply(deg,rfwd(p),ifwd(p),rx(q),ix(q),rfwd(r),ifwd(r));
    end loop;
    if idx'last > 2 then
      p := idx(idx'last); q := idx(idx'last-1);
      Multiply(deg,rx(p),ix(p),rx(q),ix(q),rbck(1),ibck(1));
      for k in 2..idx'last-2 loop
        p := k-1; q := idx(idx'last-k);
        Multiply(deg,rbck(p),ibck(p),rx(q),ix(q),rbck(k),ibck(k));
      end loop;
      if idx'last = 3 then
        p := idx(1); q := idx(3);
        Multiply(deg,rx(p),ix(p),rx(q),ix(q),rcrs(1),icrs(1));
      else
        p := idx(1); q := idx'last-3;
        Multiply(deg,rx(p),ix(p),rbck(q),ibck(q),rcrs(1),icrs(1));
        for k in 2..idx'last-3 loop
          p := k-1; q := idx'last-2-k;
          Multiply(deg,rfwd(p),ifwd(p),rbck(q),ibck(q),rcrs(k),icrs(k));
        end loop;
        p := idx'last-3; q := idx(idx'last); r := idx'last-2;
        Multiply(deg,rfwd(p),ifwd(p),rx(q),ix(q),rcrs(r),icrs(r));
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

  procedure Speel ( deg : in integer32;
                    idx : in Standard_Integer_VecVecs.VecVec;
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
          Update(deg,ryptr,iyptr,rx(q),ix(q));
          p := ryd(q); p(0) := p(0) + 1.0;
        else
          Speel(deg,rx,ix,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs);
          q := idk'last-1;
          Update(deg,ryptr,iyptr,rfwd(q),ifwd(q));
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Update(deg,ryd(q),iyd(q),rx(r),ix(r));
            Update(deg,ryd(r),iyd(r),rx(q),ix(q));
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Update(deg,ryd(q),iyd(q),rbck(r),ibck(r));
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Update(deg,ryd(q),iyd(q),rcrs(r),icrs(r));
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Update(deg,ryd(q),iyd(q),rfwd(r),ifwd(r));
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

  procedure Speel ( deg : in integer32;
                    idx : in Standard_Integer_VecVecs.VecVec;
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
          Multiply(deg,rp,ip,rx(q),ix(q),rwrk,iwrk);
          Update(deg,ryptr,iyptr,rwrk,iwrk);
          Update(deg,ryd(q),iyd(q),rp,ip);
        else
          Speel(deg,rx,ix,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs);
          q := idk'last-1;
          Multiply(deg,rp,ip,rfwd(q),ifwd(q),rwrk,iwrk);
          Update(deg,ryptr,iyptr,rwrk,iwrk);
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Multiply(deg,rp,ip,rx(r),ix(r),rwrk,iwrk);
            Update(deg,ryd(q),iyd(q),rwrk,iwrk);
            Multiply(deg,rp,ip,rx(q),ix(q),rwrk,iwrk);
            Update(deg,ryd(r),iyd(r),rwrk,iwrk);
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Multiply(deg,rp,ip,rbck(r),ibck(r),rwrk,iwrk);
            Update(deg,ryd(q),iyd(q),rwrk,iwrk);
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Multiply(deg,rp,ip,rcrs(r),icrs(r),rwrk,iwrk);
              Update(deg,ryd(q),iyd(q),rwrk,iwrk);
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Multiply(deg,rp,ip,rfwd(r),ifwd(r),rwrk,iwrk);
            Update(deg,ryd(q),iyd(q),rwrk,iwrk);
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

  procedure Multiply_Factor
              ( deg : in integer32;
                xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
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
      Multiply(deg,rcff,icff,rx(fptr),ix(fptr),racc,iacc);
    else
      rlpw := rpwx(powidx-2);  -- coefficients of higher powers
      ilpw := ipwx(powidx-2);
      Multiply(deg,rcff,icff,rlpw,ilpw,racc,iacc);
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
        Multiply(deg,rwrk,iwrk,rx(fptr),ix(fptr),racc,iacc);
      else
        rlpw := rpwx(powidx-2);  -- coefficients of higher powers
        ilpw := ipwx(powidx-2);
        Multiply(deg,rwrk,iwrk,rlpw,ilpw,racc,iacc);
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

  procedure Multiply_Power
              ( deg,multiplier : in integer32;
                rcff : in Standard_Floating_Vectors.Link_to_Vector; 
                icff : in Standard_Floating_Vectors.Link_to_Vector ) is

    factor : constant double_float := Create(multiplier);

  begin
    for i in rcff'first..deg loop
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

  procedure Speel ( deg : in integer32;
                    xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
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
            Multiply(deg,rpcf,ipcf,rx(pidx),ix(pidx),rwrk,iwrk);
            Update(deg,ryptr,iyptr,rwrk,iwrk);
            Update(deg,ryd(pidx),iyd(pidx),rpcf,ipcf);
          else
            Multiply_Factor(deg,xpk,fck,rx,ix,rpcf,ipcf,rwrk,iwrk,
                            racc,iacc,rpwt,ipwt);
            Multiply(deg,racc,iacc,rx(pidx),ix(pidx),rwrk,iwrk);
            Update(deg,ryptr,iyptr,rwrk,iwrk);
            Multiply_Power(deg,xpk(pidx),racc,iacc);
            Update(deg,ryd(pidx),iyd(pidx),racc,iacc);
          end if;
        else
          Speel(deg,rx,ix,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs);
          pidx := idk'last-1;
          if fck = null then
            Multiply(deg,rpcf,ipcf,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
          else
            Multiply_Factor(deg,xpk,fck,rx,ix,rpcf,ipcf,rwrk,iwrk,
                            racc,iacc,rpwt,ipwt);
            Multiply(deg,racc,iacc,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
          end if;
          Update(deg,ryptr,iyptr,rwrk,iwrk);
          if idk'last = 2 then
            pidx := idk(1); qidx := idk(2);
            if fck = null then
              Multiply(deg,rpcf,ipcf,rx(pidx),ix(pidx),rwrk,iwrk);
              Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
              Multiply(deg,rpcf,ipcf,rx(qidx),ix(qidx),rwrk,iwrk);
              Update(deg,ryd(pidx),iyd(pidx),rwrk,iwrk);
            else -- use the common factor in acc
              Multiply(deg,racc,iacc,rx(pidx),ix(pidx),rwrk,iwrk);
              if xpk(qidx) > 1
               then Multiply_Power(deg,xpk(qidx),rwrk,iwrk);
              end if;
              Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
              Multiply(deg,racc,iacc,rx(qidx),ix(qidx),rwrk,iwrk);
              if xpk(pidx) > 1
               then Multiply_Power(deg,xpk(pidx),rwrk,iwrk);
              end if;
              Update(deg,ryd(pidx),iyd(pidx),rwrk,iwrk);
            end if;
          else -- idk'last > 2 
            pidx := idk'last-2; qidx := idk(1);
            if fck = null then
              Multiply(deg,rpcf,ipcf,rbck(pidx),ibck(pidx),rwrk,iwrk);
              Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
              for j in idk'first+1..idk'last-1 loop
                Multiply(deg,rpcf,ipcf,rcrs(j-1),icrs(j-1),rwrk,iwrk);
                qidx := idk(j);
                Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
              end loop;
              Multiply(deg,rpcf,ipcf,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
              qidx := idk(idk'last);
              Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
            else
              Multiply(deg,racc,iacc,rbck(pidx),ibck(pidx),rwrk,iwrk);
              Multiply_Power(deg,xpk(qidx),rwrk,iwrk);
              Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
              for j in idk'first+1..idk'last-1 loop
                Multiply(deg,racc,iacc,rcrs(j-1),icrs(j-1),rwrk,iwrk);
                qidx := idk(j);
                Multiply_Power(deg,xpk(qidx),rwrk,iwrk);
                Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
              end loop;
              Multiply(deg,racc,iacc,rfwd(pidx),ifwd(pidx),rwrk,iwrk);
              qidx := idk(idk'last);
              Multiply_Power(deg,xpk(qidx),rwrk,iwrk);
              Update(deg,ryd(qidx),iyd(qidx),rwrk,iwrk);
            end if;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

-- EVALUATION OF POWER SERIES COEFFICIENTS :

  procedure EvalCoeff ( c : in Circuit; t : in double_float;
                        rct,ict : out double_float;
                        rcf : out Standard_Floating_Vectors.Vector;
                        icf : out Standard_Floating_Vectors.Vector ) is

    use Standard_Floating_Vectors;

    deg : integer32;
    rlc,ilc : Link_to_Vector;
    zr,zi : double_float;

  begin
    if c.rct = null then     -- also c.ict = null
      rct := 0.0; ict := 0.0;
    else
      deg := c.rct'last;     -- equals c.ict'last
      rct := c.rct(deg);
      ict := c.ict(deg);
      for k in reverse 0..deg-1 loop
        rct := rct*t + c.rct(k);
        ict := ict*t + c.ict(k);
      end loop;
    end if;
    for i in 1..c.nbr loop
      rlc := c.rcf(i); ilc := c.icf(i);
      deg := rlc'last;
      zr := rlc(deg); zi := ilc(deg);
      for k in reverse 0..deg-1 loop
        zr := zr*t + rlc(k);
        zi := zi*t + ilc(k);
      end loop;
      rcf(i) := zr; icf(i) := zi;
    end loop;
  end EvalCoeff;

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

  procedure EvalDiff ( deg : in integer32; c : in Circuit;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec ) is

    use Standard_Floating_Vectors;

  begin
    Speel(deg,c.xps,c.idx,c.fac,c.rcf,c.icf,rx,ix,c.rfwd,c.ifwd,c.rbck,c.ibck,
          c.rcrs,c.icrs,ryd,iyd,c.rwrk,c.iwrk,c.racc,c.iacc,rpwt,ipwt);
    if c.rct /= null and c.ict /= null
     then Update(deg,ryd(ryd'last),iyd(iyd'last),c.rct,c.ict);
    end if;
  end EvalDiff;

  procedure EvalDiff ( c : in Circuits;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                       vy : in Standard_Complex_VecVecs.VecVec;
                       vm : in Standard_Complex_VecMats.VecMat ) is

    vleft : Standard_Complex_Vectors.Link_to_Vector;
    rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
    mleft : Standard_Complex_Matrices.Link_to_Matrix;

  begin
    for i in c'range loop
      EvalDiff(c(i).all,rx,ix,rpwt,ipwt,ryd,iyd);
      rvright := ryd(rx'last+1); ivright := iyd(ix'last+1);
      for j in rvright'range loop  -- the j-th coefficient of vright is
        vleft := vy(j);  -- assigned to the j-th vector of vy at position i
        vleft(i) := Standard_Complex_Numbers.Create(rvright(j),ivright(j));
        rvright(j) := 0.0;         -- reset the value to zero
        ivright(j) := 0.0;
      end loop;
      for j in 1..rx'last loop
        rvright := ryd(j); ivright := iyd(j);
        for k in vm'range loop     -- k-th coefficient in matrix vm(k)
          mleft := vm(k);          -- the row i in vm(k) is the equation
          mleft(i,j)               -- the column j in vm(k) is the variable
            := Standard_Complex_Numbers.Create(rvright(k),ivright(k));
          rvright(k) := 0.0;       -- reset the value to zero
          ivright(k) := 0.0;
        end loop;
      end loop;
    end loop;
  end EvalDiff;

  procedure EvalDiff ( deg : in integer32; c : in Circuits;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                       vy : in Standard_Complex_VecVecs.VecVec;
                       vm : in Standard_Complex_VecMats.VecMat ) is

    vleft : Standard_Complex_Vectors.Link_to_Vector;
    rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
    mleft : Standard_Complex_Matrices.Link_to_Matrix;

  begin
    for i in c'range loop
      EvalDiff(deg,c(i).all,rx,ix,rpwt,ipwt,ryd,iyd);
      rvright := ryd(rx'last+1); ivright := iyd(ix'last+1);
      for j in rvright'first..deg loop  -- the j-th coefficient of vright is
        vleft := vy(j);   -- assigned to the j-th vector of vy at position i
        vleft(i) := Standard_Complex_Numbers.Create(rvright(j),ivright(j));
        rvright(j) := 0.0;         -- reset the value to zero
        ivright(j) := 0.0;
      end loop;
      for j in 1..rx'last loop
        rvright := ryd(j); ivright := iyd(j);
        for k in vm'first..deg loop     -- k-th coefficient in matrix vm(k)
          mleft := vm(k);            -- the row i in vm(k) is the equation
          mleft(i,j)                 -- the column j in vm(k) is the variable
            := Standard_Complex_Numbers.Create(rvright(k),ivright(k));
          rvright(k) := 0.0;       -- reset the value to zero
          ivright(k) := 0.0;
        end loop;
      end loop;
    end loop;
  end EvalDiff;

  procedure Delinearize ( vy,yv : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for k in vy'range loop
      declare
        vyk : constant Standard_Complex_Vectors.Link_to_Vector := vy(k);
        left : Standard_Complex_Vectors.Link_to_Vector;
      begin
        for i in yv'range loop  -- vyk holds k-th coefficient of all series
          left := yv(i);        -- so we assign to coefficients of series i
          left(k) := vyk(i);    -- at position k the i-th value of vyk
        end loop;
      end;
    end loop;
  end Delinearize;

  procedure Delinearize ( deg : in integer32;
                          vy,yv : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for k in vy'first..deg loop
      declare
        vyk : constant Standard_Complex_Vectors.Link_to_Vector := vy(k);
        left : Standard_Complex_Vectors.Link_to_Vector;
      begin
        for i in yv'range loop  -- vyk holds k-th coefficient of all series
          left := yv(i);        -- so we assign to coefficients of series i
          left(k) := vyk(i);    -- at position k the i-th value of vyk
        end loop;
      end;
    end loop;
  end Delinearize;

  procedure EvalDiff ( s : in System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec ) is
  begin
    EvalDiff(s.crc,rx,ix,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm);
    Delinearize(s.vy,s.yv);
  end EvalDiff;

  procedure EvalDiff ( deg : in integer32; s : in System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec ) is
  begin
    EvalDiff(deg,s.crc,rx,ix,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm);
    Delinearize(deg,s.vy,s.yv);
  end EvalDiff;

  procedure EvalDiff ( s : in Link_to_System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec ) is
  begin
    EvalDiff(s.crc,rx,ix,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm);
    Delinearize(s.vy,s.yv);
  end EvalDiff;

  procedure EvalDiff ( deg : in integer32; s : in Link_to_System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec ) is
  begin
    EvalDiff(deg,s.crc,rx,ix,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm);
    Delinearize(deg,s.vy,s.yv);
  end EvalDiff;

-- DEALLOCATORS :

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

  procedure Clear ( s : in out System ) is
  begin
    Clear(s.crc);
    Clear(s.rpwt);
    Clear(s.ipwt);
    Standard_Floating_VecVecs.Clear(s.ryd);
    Standard_Floating_VecVecs.Clear(s.iyd);
    Standard_Complex_VecVecs.Clear(s.vy);
    Standard_Complex_VecVecs.Clear(s.yv);
    Standard_Complex_VecMats.Clear(s.vm);
  end Clear;

  procedure Clear ( s : in out Link_to_System ) is

    procedure free is new unchecked_deallocation(System,Link_to_System);

  begin
    if s /= null then
      Clear(s.all);
      free(s);
    end if;
  end Clear;

  procedure Clear ( s : in out System_Array ) is
  begin
    for k in s'range loop
      Clear(s(k));
    end loop;
  end Clear;

end Standard_Coefficient_Convolutions;
