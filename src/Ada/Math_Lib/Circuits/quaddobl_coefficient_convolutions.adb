with unchecked_deallocation;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Standard_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with Exponent_Indices;
with Standard_Coefficient_Convolutions;

package body QuadDobl_Coefficient_Convolutions is

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

  function Linearized_Allocation
             ( dim,deg : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(0..deg);

  begin
    for k in 0..deg loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Linearized_Allocation;

  function Allocate_Coefficients
             ( nbq,nvr,deg : integer32 )
             return QuadDobl_Complex_VecMats.VecMat is

    res : QuadDobl_Complex_VecMats.VecMat(0..deg);

  begin
    for k in res'range loop
      declare
        mat : QuadDobl_Complex_Matrices.Matrix(1..nbq,1..nvr);
      begin
        for i in 1..nbq loop
          for j in 1..nvr loop
            mat(i,j) := QuadDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        res(k) := new QuadDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Allocate_Coefficients;

  function Create ( c : Circuits; dim,deg : integer32 ) return System is

    neq : constant integer32 := c'last;
    res : System(neq,neq+1,dim,dim+1,deg);
    degdim : constant integer32 := 4*(deg+1)-1;
  
    use Standard_Vector_Splitters;
    use QuadDobl_Vector_Splitters;

  begin
    res.crc := c;
    res.mxe := Exponent_Maxima(c,dim);
    res.rpwt := Standard_Coefficient_Convolutions.Allocate(res.mxe,degdim);
    res.ipwt := Standard_Coefficient_Convolutions.Allocate(res.mxe,degdim);
    res.ryd := Allocate_Floating_Coefficients(dim+1,degdim);
    res.iyd := Allocate_Floating_Coefficients(dim+1,degdim);
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

-- COMPUTING THE POWER TABLE :

  procedure Compute
              ( rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;

    use QuadDobl_Vector_Splitters;

  begin
    for i in xr'range loop
      if mxe(i) > 2 then
        rxpw := rpwt(i); ixpw := ipwt(i);
        Multiply(xr(i),xi(i),xr(i),xi(i),rxpw(1),ixpw(1),u,v,w);
        for k in 2..(mxe(i)-2) loop
          Multiply(rxpw(k-1),ixpw(k-1),xr(i),xi(i),rxpw(k),ixpw(k),u,v,w);
        end loop;
      end if;
    end loop;
  end Compute;

-- REVERSE MODE OF ALGORITHMIC DIFFERENTIATION :

  procedure Speel ( xr,xi : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    use QuadDobl_Vector_Splitters;

  begin
    Multiply(xr(1),xi(1),xr(2),xi(2),rfwd(1),ifwd(1),u,v,w);
    for k in 3..xr'last loop
      Multiply(rfwd(k-2),ifwd(k-2),xr(k),xi(k),rfwd(k-1),ifwd(k-1),u,v,w);
    end loop;
    if xr'last > 2 then
      Multiply(xr(xr'last),xi(xi'last),xr(xr'last-1),xi(xi'last-1),
               rbck(1),ibck(1),u,v,w);
      for k in 2..xr'last-2 loop
        Multiply(rbck(k-1),ibck(k-1),xr(xr'last-k),xi(xi'last-k),
                 rbck(k),ibck(k),u,v,w);
      end loop;
      if xr'last = 3 then
        Multiply(xr(1),xi(1),xr(3),xi(3),rcrs(1),icrs(1),u,v,w);
      else
        Multiply(xr(1),xi(1),rbck(xr'last-3),ibck(xi'last-3),
                 rcrs(1),icrs(1),u,v,w);
        for k in 2..xr'last-3 loop
          Multiply(rfwd(k-1),ifwd(k-1),rbck(xr'last-2-k),ibck(xi'last-2-k),
                   rcrs(k),icrs(k),u,v,w);
        end loop;
        Multiply(rfwd(xr'last-3),ifwd(xi'last-3),xr(xr'last),xi(xi'last),
                 rcrs(xr'last-2),icrs(xi'last-2),u,v,w);
      end if;
    end if;
  end Speel;

  procedure Speel ( xr,xi : in Standard_Floating_VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    p,q,r : integer32;

    use QuadDobl_Vector_Splitters;

  begin
    p := idx(1); q := idx(2);
    Multiply(xr(p),xi(p),xr(q),xi(q),rfwd(1),ifwd(1),u,v,w);
    for k in 3..idx'last loop
      p := k-2; q := idx(k); r := k-1;
      Multiply(rfwd(p),ifwd(p),xr(q),xi(q),rfwd(r),ifwd(r),u,v,w);
    end loop;
    if idx'last > 2 then
      p := idx(idx'last); q := idx(idx'last-1);
      Multiply(xr(p),xi(p),xr(q),xi(q),rbck(1),ibck(1),u,v,w);
      for k in 2..idx'last-2 loop
        p := k-1; q := idx(idx'last-k);
        Multiply(rbck(p),ibck(p),xr(q),xi(q),rbck(k),ibck(k),u,v,w);
      end loop;
      if idx'last = 3 then
        p := idx(1); q := idx(3);
        Multiply(xr(p),xi(p),xr(q),xi(q),rcrs(1),icrs(1),u,v,w);
      else
        p := idx(1); q := idx'last-3;
        Multiply(xr(p),xi(p),rbck(q),ibck(q),rcrs(1),icrs(1),u,v,w);
        for k in 2..idx'last-3 loop
          p := k-1; q := idx'last-2-k;
          Multiply(rfwd(p),ifwd(p),rbck(q),ibck(q),rcrs(k),icrs(k),u,v,w);
        end loop;
        p := idx'last-3; q := idx(idx'last); r := idx'last-2;
        Multiply(rfwd(p),ifwd(p),xr(q),xi(q),rcrs(r),icrs(r),u,v,w);
      end if;
    end if;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    xr,xi : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Integer_Vectors,QuadDobl_Vector_Splitters;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    ryptr : constant Standard_Floating_Vectors.Link_to_Vector
          := ryd(ryd'last);
    iyptr : constant Standard_Floating_Vectors.Link_to_Vector
          := iyd(iyd'last);
    p : Standard_Floating_Vectors.Link_to_Vector;
    q,r : integer32;
    one : constant quad_double := create(1.0);

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          q := idk(1);
          Update(ryptr,iyptr,xr(q),xi(q),u);
          p := ryd(q);
          p(0) := p(0) + hihi_part(one);
         -- p(1) := p(1) + lohi_part(one);
         -- p(2) := p(2) + hilo_part(one);
         -- p(3) := p(3) + lolo_part(one);
        else
          Speel(xr,xi,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs,u,v,w);
          q := idk'last-1;
          Update(ryptr,iyptr,rfwd(q),ifwd(q),u);
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Update(ryd(q),iyd(q),xr(r),xi(r),u);
            Update(ryd(r),iyd(r),xr(q),xi(q),u);
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Update(ryd(q),iyd(q),rbck(r),ibck(r),u);
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Update(ryd(q),iyd(q),rcrs(r),icrs(r),u);
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Update(ryd(q),iyd(q),rfwd(r),ifwd(r),u);
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    xr,xi : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Integer_Vectors,QuadDobl_Vector_Splitters;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    ryptr : constant Standard_Floating_Vectors.Link_to_Vector
          := ryd(ryd'last);
    iyptr : constant Standard_Floating_Vectors.Link_to_Vector
          := iyd(iyd'last);
    rp,ip : Standard_Floating_Vectors.Link_to_Vector;
    q,r : integer32;

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        rp := rcff(k); ip := icff(k);
        if idk'last = 1 then
          q := idk(1);
          Multiply(rp,ip,xr(q),xi(q),rwrk,iwrk,u,v,w);
          Update(ryptr,iyptr,rwrk,iwrk,u);
          Update(ryd(q),iyd(q),rp,ip,u);
        else
          Speel(xr,xi,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs,u,v,w);
          q := idk'last-1;
          Multiply(rp,ip,rfwd(q),ifwd(q),rwrk,iwrk,u,v,w);
          Update(ryptr,iyptr,rwrk,iwrk,u);
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Multiply(rp,ip,xr(r),xi(r),rwrk,iwrk,u,v,w);
            Update(ryd(q),iyd(q),rwrk,iwrk,u);
            Multiply(rp,ip,xr(q),xi(q),rwrk,iwrk,u,v,w);
            Update(ryd(r),iyd(r),rwrk,iwrk,u);
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Multiply(rp,ip,rbck(r),ibck(r),rwrk,iwrk,u,v,w);
            Update(ryd(q),iyd(q),rwrk,iwrk,u);
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Multiply(rp,ip,rcrs(r),icrs(r),rwrk,iwrk,u,v,w);
              Update(ryd(q),iyd(q),rwrk,iwrk,u);
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Multiply(rp,ip,rfwd(r),ifwd(r),rwrk,iwrk,u,v,w);
            Update(ryd(q),iyd(q),rwrk,iwrk,u);
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Multiply_Factor
              ( xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rcff,icff : in Standard_Floating_Vectors.Link_to_Vector;
                rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    rpwx,ipwx : Standard_Floating_VecVecs.Link_to_VecVec;
    rlpw,ilpw : Standard_Floating_Vectors.Link_to_Vector;
    powidx,fptr : integer32;

    use QuadDobl_Vector_Splitters;

  begin
    fptr := facidx(facidx'first);
    rpwx := rpwt(fptr); ipwx := ipwt(fptr);
    powidx := xpk(fptr);   -- power in power table
    if powidx = 2 then
      Multiply(rcff,icff,xr(fptr),xi(fptr),racc,iacc,u,v,w);
    else
      rlpw := rpwx(powidx-2);  -- coefficients of higher powers
      ilpw := ipwx(powidx-2);
      Multiply(rcff,icff,rlpw,ilpw,racc,iacc,u,v,w);
    end if;
    for k in facidx'first+1..facidx'last loop
      for i in rwrk'range loop
        rwrk(i) := racc(i); iwrk(i) := iacc(i);
      end loop;
      fptr := facidx(k);
      rpwx := rpwt(fptr); ipwx := ipwt(fptr);
      powidx := xpk(fptr);   -- power in power table
      if powidx = 2 then
        Multiply(rwrk,iwrk,xr(fptr),xi(fptr),racc,iacc,u,v,w);
      else
        rlpw := rpwx(powidx-2);  -- coefficients of higher powers
        ilpw := ipwx(powidx-2);
        Multiply(rwrk,iwrk,rlpw,ilpw,racc,iacc,u,v,w);
      end if;
    end loop;
  end Multiply_Factor;

  procedure Multiply_Power
              ( multiplier : in integer32;
                rcff : in Standard_Floating_Vectors.Link_to_Vector; 
                icff : in Standard_Floating_Vectors.Link_to_Vector ) is

    factor : constant double_float := Create(multiplier);
    nbr : quad_double;
    dim : constant integer32 := rcff'last;
    deg : constant integer32 := (dim+1)/4-1;
    idx : integer32 := rcff'first;

  begin
    for i in 0..deg loop
      nbr := create(rcff(idx),rcff(idx+1),rcff(idx+2),rcff(idx+3));
      nbr := factor*nbr;
      rcff(idx) := hihi_part(nbr);
      rcff(idx+1) := lohi_part(nbr);
      rcff(idx+2) := hilo_part(nbr);
      rcff(idx+3) := lolo_part(nbr);
      nbr := create(icff(idx),icff(idx+1),icff(idx+2),icff(idx+3));
      nbr := factor*nbr;
      icff(idx) := hihi_part(nbr);
      icff(idx+1) := lohi_part(nbr);
      icff(idx+2) := hilo_part(nbr);
      icff(idx+3) := lolo_part(nbr);
      idx := idx + 4;
    end loop;
  end Multiply_Power;

  procedure Speel
              ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                rcff,icff : in Standard_Floating_VecVecs.VecVec;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Integer_Vectors,QuadDobl_Vector_Splitters;

    idk,xpk,fck : Standard_Integer_Vectors.Link_to_Vector;
    ryptr : constant Standard_Floating_Vectors.Link_to_Vector
          := ryd(ryd'last);
    iyptr : constant Standard_Floating_Vectors.Link_to_Vector
          := iyd(iyd'last);
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
            Multiply(rpcf,ipcf,xr(pidx),xi(pidx),rwrk,iwrk,u,v,w);
            Update(ryptr,iyptr,rwrk,iwrk,u);
            Update(ryd(pidx),iyd(pidx),rpcf,ipcf,u);
          else
            Multiply_Factor(xpk,fck,xr,xi,rpcf,ipcf,rwrk,iwrk,racc,iacc,
                            rpwt,ipwt,u,v,w);
            Multiply(racc,iacc,xr(pidx),xi(pidx),rwrk,iwrk,u,v,w);
            Update(ryptr,iyptr,rwrk,iwrk,u);
            Multiply_Power(xpk(pidx),racc,iacc);
            Update(ryd(pidx),iyd(pidx),racc,iacc,u);
          end if;
        else
          Speel(xr,xi,idk.all,rfwd,ifwd,rbck,ibck,rcrs,icrs,u,v,w);
          pidx := idk'last-1;
          if fck = null then
            Multiply(rpcf,ipcf,rfwd(pidx),ifwd(pidx),rwrk,iwrk,u,v,w);
          else
            Multiply_Factor(xpk,fck,xr,xi,rpcf,ipcf,rwrk,iwrk,racc,iacc,
                            rpwt,ipwt,u,v,w);
            Multiply(racc,iacc,rfwd(pidx),ifwd(pidx),rwrk,iwrk,u,v,w);
          end if;
          Update(ryptr,iyptr,rwrk,iwrk,u);
          if idk'last = 2 then
            pidx := idk(1); qidx := idk(2);
            if fck = null then
              Multiply(rpcf,ipcf,xr(pidx),xi(pidx),rwrk,iwrk,u,v,w);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
              Multiply(rpcf,ipcf,xr(qidx),xi(qidx),rwrk,iwrk,u,v,w);
              Update(ryd(pidx),iyd(pidx),rwrk,iwrk,u);
            else -- use the common factor in acc
              Multiply(racc,iacc,xr(pidx),xi(pidx),rwrk,iwrk,u,v,w);
              if xpk(qidx) > 1
               then Multiply_Power(xpk(qidx),rwrk,iwrk);
              end if;
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
              Multiply(racc,iacc,xr(qidx),xi(qidx),rwrk,iwrk,u,v,w);
              if xpk(pidx) > 1
               then Multiply_Power(xpk(pidx),rwrk,iwrk);
              end if;
              Update(ryd(pidx),iyd(pidx),rwrk,iwrk,u);
            end if;
          else -- idk'last > 2 
            pidx := idk'last-2; qidx := idk(1);
            if fck = null then
              Multiply(rpcf,ipcf,rbck(pidx),ibck(pidx),rwrk,iwrk,u,v,w);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
              for j in idk'first+1..idk'last-1 loop
                Multiply(rpcf,ipcf,rcrs(j-1),icrs(j-1),rwrk,iwrk,u,v,w);
                qidx := idk(j);
                Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
              end loop;
              Multiply(rpcf,ipcf,rfwd(pidx),ifwd(pidx),rwrk,iwrk,u,v,w);
              qidx := idk(idk'last);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
            else
              Multiply(racc,iacc,rbck(pidx),ibck(pidx),rwrk,iwrk,u,v,w);
              Multiply_Power(xpk(qidx),rwrk,iwrk);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
              for j in idk'first+1..idk'last-1 loop
                Multiply(racc,iacc,rcrs(j-1),icrs(j-1),rwrk,iwrk,u,v,w);
                qidx := idk(j);
                Multiply_Power(xpk(qidx),rwrk,iwrk);
                Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
              end loop;
              Multiply(racc,iacc,rfwd(pidx),ifwd(pidx),rwrk,iwrk,u,v,w);
              qidx := idk(idk'last);
              Multiply_Power(xpk(qidx),rwrk,iwrk);
              Update(ryd(qidx),iyd(qidx),rwrk,iwrk,u);
            end if;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

-- EVALUATION AND DIFFERENTIATION ON CIRCUITS :

  procedure EvalDiff
              ( c : in Circuit;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Floating_Vectors,QuadDobl_Vector_Splitters;

  begin
    Speel(c.xps,c.idx,c.fac,c.rcff,c.icff,xr,xi,c.rfwd,c.ifwd,c.rbck,c.ibck,
          c.rcrs,c.icrs,ryd,iyd,c.rwrk,c.iwrk,c.racc,c.iacc,rpwt,ipwt,u,v,w);
    if c.rcst /= null and c.icst /= null
     then Update(ryd(ryd'last),iyd(iyd'last),c.rcst,c.icst,u);
    end if;
  end EvalDiff;

  procedure EvalDiff
              ( c : in Circuits;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is

    vleft : QuadDobl_Complex_Vectors.Link_to_Vector;
    rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
    mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
    rqd,iqd : quad_double;
    idx : integer32;

  begin
    for i in c'range loop
      EvalDiff(c(i).all,xr,xi,rpwt,ipwt,ryd,iyd,u,v,w);
      rvright := ryd(xr'last+1); ivright := iyd(xi'last+1);
      idx := rvright'first;
      for j in vy'range loop  -- the j-th coefficient of vright is
        vleft := vy(j);  -- assigned to the j-th vector of vy at position i
        rqd := Create(rvright(idx),rvright(idx+1),
                      rvright(idx+2),rvright(idx+3));
        iqd := Create(ivright(idx),ivright(idx+1),
                      ivright(idx+2),ivright(idx+3));
        vleft(i) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
        rvright(idx) := 0.0;   ivright(idx) := 0.0;   -- reset to zero
        rvright(idx+1) := 0.0; ivright(idx+1) := 0.0;
        rvright(idx+2) := 0.0; ivright(idx+2) := 0.0;
        rvright(idx+3) := 0.0; ivright(idx+3) := 0.0;
        idx := idx + 4;
      end loop;
      for j in 1..xr'last loop
        rvright := ryd(j); ivright := iyd(j);
        idx := rvright'first;
        for k in vm'range loop     -- k-th coefficient in matrix vm(k)
          mleft := vm(k);          -- the row i in vm(k) is the equation
                                   -- the column j in vm(k) is the variable
          rqd := Create(rvright(idx),rvright(idx+1),
                        rvright(idx+2),rvright(idx+3));
          iqd := Create(ivright(idx),ivright(idx+1),
                        ivright(idx+2),ivright(idx+3));
          mleft(i,j) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
          rvright(idx) := 0.0;   ivright(idx) := 0.0;      -- reset to zero
          rvright(idx+1) := 0.0; ivright(idx+1) := 0.0;
          rvright(idx+2) := 0.0; ivright(idx+2) := 0.0;
          rvright(idx+3) := 0.0; ivright(idx+3) := 0.0;
          idx := idx + 4;
        end loop;
      end loop;
    end loop;
  end EvalDiff;

  procedure Delinearize ( vy,yv : in QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in vy'range loop
      declare
        vyk : constant QuadDobl_Complex_Vectors.Link_to_Vector := vy(k);
        left : QuadDobl_Complex_Vectors.Link_to_Vector;
      begin
        for i in yv'range loop  -- vyk holds k-th coefficient of all series
          left := yv(i);        -- so we assign to coefficients of series i
          left(k) := vyk(i);    -- at position k the i-th value of vyk
        end loop;
      end;
    end loop;
  end Delinearize;

  procedure EvalDiff
              ( s : in System;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    EvalDiff(s.crc,xr,xi,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm,u,v,w);
    Delinearize(s.vy,s.yv);
  end EvalDiff;

  procedure EvalDiff
               ( s : in Link_to_System;
                 xr,xi : in Standard_Floating_VecVecs.VecVec;
                 u,v,w : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    EvalDiff(s.crc,xr,xi,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm,u,v,w);
    Delinearize(s.vy,s.yv);
  end EvalDiff;

-- DEALLOCATORS :

  procedure Clear ( c : in out Circuit ) is
  begin
    Standard_Integer_VecVecs.Clear(c.xps);
    Standard_Integer_VecVecs.Clear(c.idx);
    Standard_Integer_VecVecs.Clear(c.fac);
    Standard_Floating_VecVecs.Clear(c.rcff);
    Standard_Floating_VecVecs.Clear(c.icff);
    Standard_Floating_Vectors.Clear(c.rcst);
    Standard_Floating_Vectors.Clear(c.icst);
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
    Standard_Floating_VecVecVecs.Clear(s.rpwt);
    Standard_Floating_VecVecVecs.Clear(s.ipwt);
    Standard_Floating_VecVecs.Clear(s.ryd);
    Standard_Floating_VecVecs.Clear(s.iyd);
    QuadDobl_Complex_VecVecs.Clear(s.vy);
    QuadDobl_Complex_VecVecs.Clear(s.yv);
    QuadDobl_Complex_VecMats.Clear(s.vm);
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

end QuadDobl_Coefficient_Convolutions;
