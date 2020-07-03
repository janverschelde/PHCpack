with unchecked_deallocation;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Double_Double_Numbers;               use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with Exponent_Indices;
with Standard_Coefficient_Convolutions; -- for power table allocations

package body DoblDobl_Coefficient_Convolutions is

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
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(0..deg);

  begin
    for k in 0..deg loop
      declare
        cff : constant DoblDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Linearized_Allocation;

  function Allocate_Coefficients
             ( nbq,nvr,deg : integer32 )
             return DoblDobl_Complex_VecMats.VecMat is

    res : DoblDobl_Complex_VecMats.VecMat(0..deg);

  begin
    for k in res'range loop
      declare
        mat : DoblDobl_Complex_Matrices.Matrix(1..nbq,1..nvr);
      begin
        for i in 1..nbq loop
          for j in 1..nvr loop
            mat(i,j) := DoblDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        res(k) := new DoblDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Allocate_Coefficients;

  function Create ( c : Circuits; dim,deg : integer32 ) return System is

    neq : constant integer32 := c'last;
    res : System(neq,neq+1,dim,dim+1,deg);

    use Standard_Vector_Splitters;
    use DoblDobl_Vector_Splitters;

  begin
    res.crc := c;
    res.mxe := Exponent_Maxima(c,dim);
    res.rhpwt := Standard_Coefficient_Convolutions.Allocate(res.mxe,deg);
    res.ihpwt := Standard_Coefficient_Convolutions.Allocate(res.mxe,deg);
    res.rlpwt := Standard_Coefficient_Convolutions.Allocate(res.mxe,deg);
    res.ilpwt := Standard_Coefficient_Convolutions.Allocate(res.mxe,deg);
    res.rhyd := Allocate_Floating_Coefficients(dim+1,deg);
    res.ihyd := Allocate_Floating_Coefficients(dim+1,deg);
    res.rlyd := Allocate_Floating_Coefficients(dim+1,deg);
    res.ilyd := Allocate_Floating_Coefficients(dim+1,deg);
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
              ( rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec ) is

    rhxpw,ihxpw,rlxpw,ilxpw : Standard_Floating_VecVecs.Link_to_VecVec;

    use DoblDobl_Vector_Splitters;

  begin
    for i in rhx'range loop
      if mxe(i) > 2 then
        rhxpw := rhpwt(i); ihxpw := ihpwt(i);
        rlxpw := rlpwt(i); ilxpw := ilpwt(i);
        Multiply(rhx(i),ihx(i),rlx(i),ilx(i),rhx(i),ihx(i),rlx(i),ilx(i),
                 rhxpw(1),ihxpw(1),rlxpw(1),ilxpw(1));
        for k in 2..(mxe(i)-2) loop
          Multiply(rhxpw(k-1),ihxpw(k-1),rlxpw(k-1),ilxpw(k-1),
                   rhx(i),ihx(i),rlx(i),ilx(i),
                   rhxpw(k),ihxpw(k),rlxpw(k),ilxpw(k));
        end loop;
      end if;
    end loop;
  end Compute;

-- REVERSE MODE OF ALGORITHMIC DIFFERENTIATION :

  procedure Speel
              ( rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec ) is

    use DoblDobl_Vector_Splitters;

  begin
    Multiply(rhx(1),ihx(1),rlx(1),ilx(1),rhx(2),ihx(2),rlx(2),ilx(2),
             rhfwd(1),ihfwd(1),rlfwd(1),ilfwd(1));
    for k in 3..rhx'last loop
      Multiply(rhfwd(k-2),ihfwd(k-2),rlfwd(k-2),ilfwd(k-2),
               rhx(k),ihx(k),rlx(k),ilx(k),
               rhfwd(k-1),ihfwd(k-1),rlfwd(k-1),ilfwd(k-1));
    end loop;
    if rhx'last > 2 then
      Multiply(rhx(rhx'last),ihx(ihx'last),rlx(rhx'last),ilx(ihx'last),
               rhx(rhx'last-1),ihx(ihx'last-1),rlx(rhx'last-1),ilx(ihx'last-1),
               rhbck(1),ihbck(1),rlbck(1),ilbck(1));
      for k in 2..rhx'last-2 loop
        Multiply(rhbck(k-1),ihbck(k-1),rlbck(k-1),ilbck(k-1),
                 rhx(rhx'last-k),ihx(ihx'last-k),
                 rlx(rlx'last-k),ilx(ilx'last-k),
                 rhbck(k),ihbck(k),rlbck(k),ilbck(k));
      end loop;
      if rhx'last = 3 then
        Multiply(rhx(1),ihx(1),rlx(1),ilx(1),rhx(3),ihx(3),rlx(3),ilx(3),
                 rhcrs(1),ihcrs(1),rlcrs(1),ilcrs(1));
      else
        Multiply(rhx(1),ihx(1),rlx(1),ilx(1),
                 rhbck(rhx'last-3),ihbck(ihx'last-3),
                 rlbck(rlx'last-3),ilbck(ilx'last-3),
                 rhcrs(1),ihcrs(1),rlcrs(1),ilcrs(1));
        for k in 2..rhx'last-3 loop
          Multiply(rhfwd(k-1),ihfwd(k-1),rlfwd(k-1),ilfwd(k-1),
                   rhbck(rhx'last-2-k),ihbck(ihx'last-2-k),
                   rlbck(rlx'last-2-k),ilbck(ilx'last-2-k),
                   rhcrs(k),ihcrs(k),rlcrs(k),ilcrs(k));
        end loop;
        Multiply(rhfwd(rhx'last-3),ihfwd(ihx'last-3),
                 rlfwd(rlx'last-3),ilfwd(ilx'last-3),
                 rhx(rhx'last),ihx(ihx'last),rlx(rlx'last),ilx(ilx'last),
                 rhcrs(rhx'last-2),ihcrs(ihx'last-2),
                 rlcrs(rlx'last-2),ilcrs(ilx'last-2));
      end if;
    end if;
  end Speel;

  procedure Speel
              ( rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                idx : in Standard_Integer_Vectors.Vector;
                rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec ) is

    p,q,r : integer32;

    use DoblDobl_Vector_Splitters;

  begin
    p := idx(1); q := idx(2);
    Multiply(rhx(p),ihx(p),rlx(p),ilx(p),rhx(q),ihx(q),rlx(q),ilx(q),
             rhfwd(1),ihfwd(1),rlfwd(1),ilfwd(1));
    for k in 3..idx'last loop
      p := k-2; q := idx(k); r := k-1;
      Multiply(rhfwd(p),ihfwd(p),rlfwd(p),ilfwd(p),
               rhx(q),ihx(q),rlx(q),ilx(q),
               rhfwd(r),ihfwd(r),rlfwd(r),ilfwd(r));
    end loop;
    if idx'last > 2 then
      p := idx(idx'last); q := idx(idx'last-1);
      Multiply(rhx(p),ihx(p),rlx(p),ilx(p),rhx(q),ihx(q),rlx(q),ilx(q),
               rhbck(1),ihbck(1),rlbck(1),ilbck(1));
      for k in 2..idx'last-2 loop
        p := k-1; q := idx(idx'last-k);
        Multiply(rhbck(p),ihbck(p),rlbck(p),ilbck(p),
                 rhx(q),ihx(q),rlx(q),ilx(q),
                 rhbck(k),ihbck(k),rlbck(k),ilbck(k));
      end loop;
      if idx'last = 3 then
        p := idx(1); q := idx(3);
        Multiply(rhx(p),ihx(p),rlx(p),ilx(p),
                 rhx(q),ihx(q),rlx(q),ilx(q),
                 rhcrs(1),ihcrs(1),rlcrs(1),ilcrs(1));
      else
        p := idx(1); q := idx'last-3;
        Multiply(rhx(p),ihx(p),rlx(p),ilx(p),
                 rhbck(q),ihbck(q),rlbck(q),ilbck(q),
                 rhcrs(1),ihcrs(1),rlcrs(1),ilcrs(1));
        for k in 2..idx'last-3 loop
          p := k-1; q := idx'last-2-k;
          Multiply(rhfwd(p),ihfwd(p),rlfwd(p),ilfwd(p),
                   rhbck(q),ihbck(q),rlbck(q),ilbck(q),
                   rhcrs(k),ihcrs(k),rlcrs(k),ilcrs(k));
        end loop;
        p := idx'last-3; q := idx(idx'last); r := idx'last-2;
        Multiply(rhfwd(p),ihfwd(p),rlfwd(p),ilfwd(p),
                 rhx(q),ihx(q),rlx(q),ilx(q),
                 rhcrs(r),ihcrs(r),rlcrs(r),ilcrs(r));
      end if;
    end if;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                    rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                    rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                    rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                    rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                    rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                    rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec;
                    rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                    rlyd,ilyd : in Standard_Floating_VecVecs.VecVec ) is

    use Standard_Integer_Vectors,DoblDobl_Vector_Splitters;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    rhyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := rhyd(rhyd'last);
    ihyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := ihyd(ihyd'last);
    rlyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := rlyd(rlyd'last);
    ilyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := ilyd(ilyd'last);
    p : Standard_Floating_Vectors.Link_to_Vector;
    q,r : integer32;
    one : constant double_double := create(1.0);

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          q := idk(1);
          Update(rhyptr,ihyptr,rlyptr,ilyptr,rhx(q),ihx(q),rlx(q),ilx(q));
          p := rhyd(q); p(0) := p(0) + hi_part(one);
          p := rlyd(q); p(0) := p(0) + lo_part(one);
        else
          Speel(rhx,ihx,rlx,ilx,idk.all,rhfwd,ihfwd,rlfwd,ilfwd,
                rhbck,ihbck,rlbck,ilbck,rhcrs,ihcrs,rlcrs,ilcrs);
          q := idk'last-1;
          Update(rhyptr,ihyptr,rlyptr,ilyptr,
                 rhfwd(q),ihfwd(q),rlfwd(q),ilfwd(q));
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),
                   rhx(r),ihx(r),rlx(r),ilx(r));
            Update(rhyd(r),ihyd(r),rlyd(r),ilyd(r),
                   rhx(q),ihx(q),rlx(q),ilx(q));
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),
                   rhbck(r),ihbck(r),rlbck(r),ilbck(r));
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),
                     rhcrs(r),ihcrs(r),rlcrs(r),ilcrs(r));
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),
                   rhfwd(r),ihfwd(r),rlfwd(r),ilfwd(r));
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rhcff,ihcff : in Standard_Floating_VecVecs.VecVec;
                    rlcff,ilcff : in Standard_Floating_VecVecs.VecVec;
                    rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                    rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                    rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                    rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                    rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                    rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                    rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec;
                    rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                    rlyd,ilyd : in Standard_Floating_VecVecs.VecVec;
                    rhwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    rlwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    ihwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    ilwrk : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Integer_Vectors,DoblDobl_Vector_Splitters;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    rhyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := rhyd(rhyd'last);
    ihyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := ihyd(ihyd'last);
    rlyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := rlyd(rlyd'last);
    ilyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := ilyd(ilyd'last);
    rhp,ihp,rlp,ilp : Standard_Floating_Vectors.Link_to_Vector;
    q,r : integer32;

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        rhp := rhcff(k); ihp := ihcff(k);
        rlp := rlcff(k); ilp := ilcff(k);
        if idk'last = 1 then
          q := idk(1);
          Multiply(rhp,ihp,rlp,ilp,rhx(q),ihx(q),rlx(q),ilx(q),
                   rhwrk,ihwrk,rlwrk,ilwrk);
          Update(rhyptr,ihyptr,rlyptr,ilyptr,rhwrk,ihwrk,rlwrk,ilwrk);
          Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),rhp,ihp,rlp,ilp);
        else
          Speel(rhx,ihx,rlx,ilx,idk.all,rhfwd,ihfwd,rlfwd,ilfwd,
                rhbck,ihbck,rlbck,ilbck,rhcrs,ihcrs,rlcrs,ilcrs);
          q := idk'last-1;
          Multiply(rhp,ihp,rlp,ilp,rhfwd(q),ihfwd(q),rlfwd(q),ilfwd(q),
                   rhwrk,ihwrk,rlwrk,ilwrk);
          Update(rhyptr,ihyptr,rlyptr,ilyptr,rhwrk,ihwrk,rlwrk,ilwrk);
          if idk'last = 2 then
            q := idk(2); r := idk(1);
            Multiply(rhp,ihp,rlp,ilp,rhx(r),ihx(r),rlx(r),ilx(r),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),rhwrk,ihwrk,rlwrk,ilwrk);
            Multiply(rhp,ihp,rlp,ilp,rhx(q),ihx(q),rlx(q),ilx(q),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            Update(rhyd(r),ihyd(r),rlyd(r),ilyd(r),rhwrk,ihwrk,rlwrk,ilwrk);
          else -- idk'last > 2 
            q := idk(1); r := idk'last-2;
            Multiply(rhp,ihp,rlp,ilp,rhbck(r),ihbck(r),rlbck(r),ilbck(r),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),rhwrk,ihwrk,rlwrk,ilwrk);
            for j in idk'first+1..idk'last-1 loop
              q := idk(j); r := j-1;
              Multiply(rhp,ihp,rlp,ilp,rhcrs(r),ihcrs(r),rlcrs(r),ilcrs(r),
                       rhwrk,ihwrk,rlwrk,ilwrk);
              Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),rhwrk,ihwrk,rlwrk,ilwrk);
            end loop;
            q := idk(idk'last); r := idk'last-2;
            Multiply(rhp,ihp,rlp,ilp,rhfwd(r),ihfwd(r),rlfwd(r),ilfwd(r),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            Update(rhyd(q),ihyd(q),rlyd(q),ilyd(q),rhwrk,ihwrk,rlwrk,ilwrk);
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Multiply_Factor
              ( xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
                rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhcff,ihcff : in Standard_Floating_Vectors.Link_to_Vector;
                rlcff,ilcff : in Standard_Floating_Vectors.Link_to_Vector;
                rhwrk,ihwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rlwrk,ilwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rhacc,ihacc : in Standard_Floating_Vectors.Link_to_Vector;
                rlacc,ilacc : in Standard_Floating_Vectors.Link_to_Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec ) is

    rhpwx,ihpwx,rlpwx,ilpwx : Standard_Floating_VecVecs.Link_to_VecVec;
    rhlpw,ihlpw,rllpw,illpw : Standard_Floating_Vectors.Link_to_Vector;
    powidx,fptr : integer32;

    use DoblDobl_Vector_Splitters;

  begin
    fptr := facidx(facidx'first);
    rhpwx := rhpwt(fptr); rlpwx := rlpwt(fptr);
    ihpwx := ihpwt(fptr); ilpwx := ilpwt(fptr);
    powidx := xpk(fptr);   -- power in power table
    if powidx = 2 then
      Multiply(rhcff,ihcff,rlcff,ilcff,rhx(fptr),ihx(fptr),
               rlx(fptr),ilx(fptr),rhacc,ihacc,rlacc,ilacc);
    else
      rhlpw := rhpwx(powidx-2);  -- coefficients of higher powers
      rllpw := rlpwx(powidx-2); ihlpw := ihpwx(powidx-2);
      illpw := ilpwx(powidx-2);
      Multiply(rhcff,ihcff,rlcff,ilcff,rhlpw,ihlpw,rllpw,illpw,
               rhacc,ihacc,rlacc,ilacc);
    end if;
    for k in facidx'first+1..facidx'last loop
      for i in rhwrk'range loop
        rhwrk(i) := rhacc(i); ihwrk(i) := ihacc(i);
        rlwrk(i) := rlacc(i); ilwrk(i) := ilacc(i);
      end loop;
      fptr := facidx(k);
      rhpwx := rhpwt(fptr); rlpwx := rlpwt(fptr);
      ihpwx := ihpwt(fptr); ilpwx := ilpwt(fptr);
      powidx := xpk(fptr);   -- power in power table
      if powidx = 2 then
        Multiply(rhwrk,ihwrk,rlwrk,ilwrk,rhx(fptr),ihx(fptr),
                 rlx(fptr),ilx(fptr),rhacc,ihacc,rlacc,ilacc);
      else
        rhlpw := rhpwx(powidx-2);  -- coefficients of higher powers
        ihlpw := ihpwx(powidx-2); rllpw := rlpwx(powidx-2);
        illpw := ilpwx(powidx-2);
        Multiply(rhwrk,ihwrk,rlwrk,ilwrk,rhlpw,ihlpw,rllpw,illpw,
                 rhacc,ihacc,rlacc,ilacc);
      end if;
    end loop;
  end Multiply_Factor;

  procedure Multiply_Power
              ( multiplier : in integer32;
                rhcff : in Standard_Floating_Vectors.Link_to_Vector; 
                ihcff : in Standard_Floating_Vectors.Link_to_Vector; 
                rlcff : in Standard_Floating_Vectors.Link_to_Vector; 
                ilcff : in Standard_Floating_Vectors.Link_to_Vector ) is

    factor : constant double_float := Create(multiplier);
    dd : double_double;

  begin
    for i in rhcff'range loop
      dd := create(rhcff(i),rlcff(i)); dd := factor*dd;
      rhcff(i) := hi_part(dd); rlcff(i) := lo_part(dd);
      dd := create(ihcff(i),ilcff(i)); dd := factor*dd;
      ihcff(i) := hi_part(dd); ilcff(i) := lo_part(dd);
     -- rhcff(i) := factor*rhcff(i); ihcff(i) := factor*ihcff(i);
     -- rlcff(i) := factor*rlcff(i); ilcff(i) := factor*ilcff(i);
    end loop;
  end Multiply_Power;

  procedure Speel
              ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                rhcff,ihcff : in Standard_Floating_VecVecs.VecVec;
                rlcff,ilcff : in Standard_Floating_VecVecs.VecVec;
                rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec;
                rhyd,ihyd,rlyd,ilyd : in Standard_Floating_VecVecs.VecVec;
                rhwrk,ihwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rlwrk,ilwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rhacc,ihacc : in Standard_Floating_Vectors.Link_to_Vector;
                rlacc,ilacc : in Standard_Floating_Vectors.Link_to_Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec ) is

    use Standard_Integer_Vectors,DoblDobl_Vector_Splitters;

    idk,xpk,fck : Standard_Integer_Vectors.Link_to_Vector;
    rhyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := rhyd(rhyd'last);
    ihyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := ihyd(ihyd'last);
    rlyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := rlyd(rlyd'last);
    ilyptr : constant Standard_Floating_Vectors.Link_to_Vector
           := ilyd(ilyd'last);
    rhpcf,ihpcf,rlpcf,ilpcf : Standard_Floating_Vectors.Link_to_Vector;
    pidx,qidx : integer32;

  begin
    for k in idx'range loop
      idk := idx(k);           -- the k-th exponent index 
      if idk /= null then
        xpk := xps(k);         -- the k-th exponent vector
        fck := fac(k);         -- the k-th factor index
        rhpcf := rhcff(k); ihpcf := ihcff(k);
        rlpcf := rlcff(k); ilpcf := ilcff(k);
        if idk'last = 1 then
          pidx := idk(1);
          if fck = null then
            Multiply(rhpcf,ihpcf,rlpcf,ilpcf,rhx(pidx),ihx(pidx),
                     rlx(pidx),ilx(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
            Update(rhyptr,ihyptr,rlyptr,ilyptr,rhwrk,ihwrk,rlwrk,ilwrk);
            Update(rhyd(pidx),ihyd(pidx),rlyd(pidx),ilyd(pidx),
                   rhpcf,ihpcf,ilpcf,ilpcf);
          else
            Multiply_Factor(xpk,fck,rhx,ihx,rlx,ilx,rhpcf,ihpcf,rlpcf,ilpcf,
                            rhwrk,ihwrk,rlwrk,ilwrk,rhacc,ihacc,rlacc,ilacc,
                            rhpwt,ihpwt,rlpwt,ilpwt);
            Multiply(rhacc,ihacc,rlacc,ilacc,rhx(pidx),ihx(pidx),
                     rlx(pidx),ilx(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
            Update(rhyptr,ihyptr,rlyptr,ilyptr,rhwrk,ihwrk,rlwrk,ilwrk);
            Multiply_Power(xpk(pidx),rhacc,ihacc,rlacc,ilacc);
            Update(rhyd(pidx),ihyd(pidx),rlyd(pidx),ilyd(pidx),
                   rhacc,ihacc,rlacc,ilacc);
          end if;
        else
          Speel(rhx,ihx,rlx,ilx,idk.all,rhfwd,ihfwd,rlfwd,ilfwd,
                rhbck,ihbck,rlbck,ilbck,rhcrs,ihcrs,rlcrs,ilcrs);
          pidx := idk'last-1;
          if fck = null then
            Multiply(rhpcf,ihpcf,rlpcf,ilpcf,rhfwd(pidx),ihfwd(pidx),
                     rlfwd(pidx),ilfwd(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
          else
            Multiply_Factor(xpk,fck,rhx,ihx,rlx,ilx,rhpcf,ihpcf,rlpcf,ilpcf,
                            rhwrk,ihwrk,rlwrk,ilwrk,rhacc,ihacc,rlacc,ilacc,
                            rhpwt,ihpwt,rlpwt,ilpwt);
            Multiply(rhacc,ihacc,rlacc,ilacc,rhfwd(pidx),ihfwd(pidx),
                     rlfwd(pidx),ilfwd(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
          end if;
          Update(rhyptr,ihyptr,rlyptr,ilyptr,rhwrk,ihwrk,rlwrk,ilwrk);
          if idk'last = 2 then
            pidx := idk(1); qidx := idk(2);
            if fck = null then
              Multiply(rhpcf,ihpcf,rlpcf,ilpcf,rhx(pidx),ihx(pidx),
                       rlx(pidx),ilx(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
              Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
              Multiply(rhpcf,ihpcf,rlpcf,ilpcf,rhx(qidx),ihx(qidx),
                       rlx(qidx),ilx(qidx),rhwrk,ihwrk,rlwrk,ilwrk);
              Update(rhyd(pidx),ihyd(pidx),rlyd(pidx),ilyd(pidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            else -- use the common factor in acc
              Multiply(rhacc,ihacc,rlacc,ilacc,rhx(pidx),ihx(pidx),
                       rlx(pidx),ilx(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
              if xpk(qidx) > 1
               then Multiply_Power(xpk(qidx),rhwrk,ihwrk,rlwrk,ilwrk);
              end if;
              Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
              Multiply(rhacc,ihacc,rlacc,ilacc,rhx(qidx),ihx(qidx),
                       rlx(qidx),ilx(qidx),rhwrk,ihwrk,rlwrk,ilwrk);
              if xpk(pidx) > 1
               then Multiply_Power(xpk(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
              end if;
              Update(rhyd(pidx),ihyd(pidx),rlyd(pidx),ilyd(pidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            end if;
          else -- idk'last > 2 
            pidx := idk'last-2; qidx := idk(1);
            if fck = null then
              Multiply(rhpcf,ihpcf,rlpcf,ilpcf,rhbck(pidx),ihbck(pidx),
                       rlbck(pidx),ilbck(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
              Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
              for j in idk'first+1..idk'last-1 loop
                Multiply(rhpcf,ihpcf,rlpcf,ilpcf,rhcrs(j-1),ihcrs(j-1),
                         rlcrs(j-1),ilcrs(j-1),rhwrk,ihwrk,rlwrk,ilwrk);
                qidx := idk(j);
                Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                       rhwrk,ihwrk,rlwrk,ilwrk);
              end loop;
              Multiply(rhpcf,ihpcf,rlpcf,ilpcf,rhfwd(pidx),ihfwd(pidx),
                       rlfwd(pidx),ilfwd(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
              qidx := idk(idk'last);
              Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            else
              Multiply(rhacc,ihacc,rlacc,ilacc,rhbck(pidx),ihbck(pidx),
                       rlbck(pidx),ilbck(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
              Multiply_Power(xpk(qidx),rhwrk,ihwrk,rlwrk,ilwrk);
              Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
              for j in idk'first+1..idk'last-1 loop
                Multiply(rhacc,ihacc,rlacc,ilacc,rhcrs(j-1),ihcrs(j-1),
                         rlcrs(j-1),ilcrs(j-1),rhwrk,ihwrk,rlwrk,ilwrk);
                qidx := idk(j);
                Multiply_Power(xpk(qidx),rhwrk,ihwrk,rlwrk,ilwrk);
                Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                       rhwrk,ihwrk,rlwrk,ilwrk);
              end loop;
              Multiply(rhacc,ihacc,rlacc,ilacc,rhfwd(pidx),ihfwd(pidx),
                       rlfwd(pidx),ilfwd(pidx),rhwrk,ihwrk,rlwrk,ilwrk);
              qidx := idk(idk'last);
              Multiply_Power(xpk(qidx),rhwrk,ihwrk,rlwrk,ilwrk);
              Update(rhyd(qidx),ihyd(qidx),rlyd(qidx),ilyd(qidx),
                     rhwrk,ihwrk,rlwrk,ilwrk);
            end if;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

-- EVALUATION AND DIFFERENTIATION ON CIRCUITS :

  procedure EvalDiff
              ( c : in Circuit;
                rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                rlyd,ilyd : in Standard_Floating_VecVecs.VecVec ) is

    use Standard_Floating_Vectors,DoblDobl_Vector_Splitters;

  begin
    Speel(c.xps,c.idx,c.fac,c.rhcf,c.ihcf,c.rlcf,c.ilcf,rhx,ihx,rlx,ilx,
          c.rhfwd,c.ihfwd,c.rlfwd,c.ilfwd,c.rhbck,c.ihbck,c.rlbck,c.ilbck,
          c.rhcrs,c.ihcrs,c.rlcrs,c.ilcrs,rhyd,ihyd,rlyd,ilyd,
          c.rhwrk,c.ihwrk,c.rlwrk,c.ilwrk,c.rhacc,c.ihacc,c.rlacc,c.ilacc,
          rhpwt,ihpwt,rlpwt,ilpwt);
    if c.rhct /= null and c.ihct /= null and
       c.rlct /= null and c.ilct /= null then
      Update(rhyd(rhyd'last),ihyd(ihyd'last),rlyd(rlyd'last),ilyd(ilyd'last),
             c.rhct,c.ihct,c.rlct,c.ilct);
    end if;
  end EvalDiff;

  procedure EvalDiff
              ( c : in Circuits;
                rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                rlyd,ilyd : in Standard_Floating_VecVecs.VecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat ) is

    vleft : DoblDobl_Complex_Vectors.Link_to_Vector;
    rhvright,ihvright : Standard_Floating_Vectors.Link_to_Vector;
    rlvright,ilvright : Standard_Floating_Vectors.Link_to_Vector;
    mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
    rdd,idd : double_double;

  begin
    for i in c'range loop
      EvalDiff(c(i).all,rhx,ihx,rlx,ilx,rhpwt,ihpwt,rlpwt,ilpwt,
               rhyd,ihyd,rlyd,ilyd);
      rhvright := rhyd(rhx'last+1); ihvright := ihyd(ihx'last+1);
      rlvright := rlyd(rlx'last+1); ilvright := ilyd(ilx'last+1);
      for j in rhvright'range loop  -- the j-th coefficient of vright is
        vleft := vy(j);  -- assigned to the j-th vector of vy at position i
        rdd := Create(rhvright(j),rlvright(j));
        idd := Create(ihvright(j),ilvright(j));
        vleft(i) := DoblDobl_Complex_Numbers.Create(rdd,idd);
        rhvright(j) := 0.0;         -- reset the value to zero
        ihvright(j) := 0.0; rlvright(j) := 0.0; ilvright(j) := 0.0;
      end loop;
      for j in 1..rhx'last loop
        rhvright := rhyd(j); ihvright := ihyd(j);
        rlvright := rlyd(j); ilvright := ilyd(j);
        for k in vm'range loop     -- k-th coefficient in matrix vm(k)
          mleft := vm(k);          -- the row i in vm(k) is the equation
                                   -- the column j in vm(k) is the variable
          rdd := Create(rhvright(k),rlvright(k));
          idd := Create(ihvright(k),ilvright(k));
          mleft(i,j) := DoblDobl_Complex_Numbers.Create(rdd,idd);
          rhvright(k) := 0.0;       -- reset the value to zero
          ihvright(k) := 0.0; rlvright(k) := 0.0; ilvright(k) := 0.0;
        end loop;
      end loop;
    end loop;
  end EvalDiff;

  procedure Delinearize ( vy,yv : in DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in vy'range loop
      declare
        vyk : constant DoblDobl_Complex_Vectors.Link_to_Vector := vy(k);
        left : DoblDobl_Complex_Vectors.Link_to_Vector;
      begin
        for i in yv'range loop  -- vyk holds k-th coefficient of all series
          left := yv(i);        -- so we assign to coefficients of series i
          left(k) := vyk(i);    -- at position k the i-th value of vyk
        end loop;
      end;
    end loop;
  end Delinearize;

  procedure EvalDiff ( s : in System;
                       rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                       rlx,ilx : in Standard_Floating_VecVecs.VecVec ) is
  begin
    EvalDiff(s.crc,rhx,ihx,rlx,ilx,s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,
             s.rhyd,s.ihyd,s.rlyd,s.ilyd,s.vy,s.vm);
    Delinearize(s.vy,s.yv);
  end EvalDiff;

  procedure EvalDiff ( s : in Link_to_System;
                       rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                       rlx,ilx : in Standard_Floating_VecVecs.VecVec ) is
  begin
    EvalDiff(s.crc,rhx,ihx,rlx,ilx,s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,
             s.rhyd,s.ihyd,s.rlyd,s.ilyd,s.vy,s.vm);
    Delinearize(s.vy,s.yv);
  end EvalDiff;

-- DEALLOCATORS :

  procedure Clear ( c : in out Circuit ) is
  begin
    Standard_Integer_VecVecs.Clear(c.xps);
    Standard_Integer_VecVecs.Clear(c.idx);
    Standard_Integer_VecVecs.Clear(c.fac);
    Standard_Floating_VecVecs.Clear(c.rhcf);
    Standard_Floating_VecVecs.Clear(c.ihcf);
    Standard_Floating_VecVecs.Clear(c.rlcf);
    Standard_Floating_VecVecs.Clear(c.ilcf);
    Standard_Floating_Vectors.Clear(c.rhct);
    Standard_Floating_Vectors.Clear(c.ihct);
    Standard_Floating_Vectors.Clear(c.rlct);
    Standard_Floating_Vectors.Clear(c.ilct);
    Standard_Floating_VecVecs.Clear(c.rhfwd);
    Standard_Floating_VecVecs.Clear(c.ihfwd);
    Standard_Floating_VecVecs.Clear(c.rlfwd);
    Standard_Floating_VecVecs.Clear(c.ilfwd);
    Standard_Floating_VecVecs.Clear(c.rhbck);
    Standard_Floating_VecVecs.Clear(c.ihbck);
    Standard_Floating_VecVecs.Clear(c.rlbck);
    Standard_Floating_VecVecs.Clear(c.ilbck);
    Standard_Floating_VecVecs.Clear(c.rhcrs);
    Standard_Floating_VecVecs.Clear(c.ihcrs);
    Standard_Floating_VecVecs.Clear(c.rlcrs);
    Standard_Floating_VecVecs.Clear(c.ilcrs);
    Standard_Floating_Vectors.Clear(c.rhwrk);
    Standard_Floating_Vectors.Clear(c.ihwrk);
    Standard_Floating_Vectors.Clear(c.rlwrk);
    Standard_Floating_Vectors.Clear(c.ilwrk);
    Standard_Floating_Vectors.Clear(c.rhacc);
    Standard_Floating_Vectors.Clear(c.ihacc);
    Standard_Floating_Vectors.Clear(c.rlacc);
    Standard_Floating_Vectors.Clear(c.ilacc);
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
    Standard_Floating_VecVecVecs.Clear(s.rhpwt);
    Standard_Floating_VecVecVecs.Clear(s.ihpwt);
    Standard_Floating_VecVecVecs.Clear(s.rlpwt);
    Standard_Floating_VecVecVecs.Clear(s.ilpwt);
    Standard_Floating_VecVecs.Clear(s.rhyd);
    Standard_Floating_VecVecs.Clear(s.ihyd);
    Standard_Floating_VecVecs.Clear(s.rlyd);
    Standard_Floating_VecVecs.Clear(s.ilyd);
    DoblDobl_Complex_VecVecs.Clear(s.vy);
    DoblDobl_Complex_VecVecs.Clear(s.yv);
    DoblDobl_Complex_VecMats.Clear(s.vm);
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

end DoblDobl_Coefficient_Convolutions;
