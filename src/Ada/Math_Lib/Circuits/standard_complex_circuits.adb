with unchecked_deallocation;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Exponent_Indices;

package body Standard_Complex_Circuits is

-- CONSTRUCTORS :

  function Allocate ( nbr,dim : integer32 ) return Circuit is

    res : Circuit(nbr);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    forward : constant Standard_Complex_Vectors.Vector(1..dim-1)
            := (1..dim-1 => zero);
    backward : constant Standard_Complex_Vectors.Vector(1..dim-2)
             := (1..dim-2 => zero);
    cross : constant Standard_Complex_Vectors.Vector(1..dim-2)
          := (1..dim-2 => zero);

  begin
    res.dim := dim;
    res.fwd := new Standard_Complex_Vectors.Vector'(forward);
    res.bck := new Standard_Complex_Vectors.Vector'(backward);
    res.crs := new Standard_Complex_Vectors.Vector'(cross);
    return res;
  end Allocate;

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

  function Create ( c : Circuits; dim : integer32 ) return System is

    res : System(c'last,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number 
         := Standard_Complex_Numbers.Create(0.0);

  begin
    res.crc := c;
    res.mxe := Exponent_Maxima(c,dim);
    res.pwt := Allocate(res.mxe);
    res.yd := new Standard_Complex_Vectors.Vector'(0..dim => zero);
    return res;
  end Create;

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUITS :

  procedure EvalDiff
              ( s : in out System;
                x : in Standard_Complex_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,x,s.pwt);
    EvalDiff(s.crc,x,s.yd,s.pwt,s.fx,s.jm);
  end EvalDiff;

  procedure EvalDiff
              ( s : in Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,x,s.pwt);
    EvalDiff(s.crc,x,s.yd,s.pwt,s.fx,s.jm);
  end EvalDiff;

  procedure EvalDiff
              ( c : in Circuits;
                x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                pwt : in Standard_Complex_VecVecs.VecVec;
                fx : out Standard_Complex_Vectors.Vector;
                jm : out Standard_Complex_Matrices.Matrix ) is
  begin
    for i in c'range loop
      Speel(c(i).all,x,yd,pwt);
      fx(i) := yd(0);
      for j in jm'range(2) loop
        jm(i,j) := yd(j);
        yd(j) := Standard_Complex_Numbers.Create(0.0);
      end loop;
    end loop;
  end EvalDiff;

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUIT :

  procedure Speel ( c : in Circuit;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                    h : out Standard_Complex_Matrices.Matrix ) is
  begin
    Speel(c.xps,c.cff,c.cst,x,yd,c.fwd,c.bck,c.crs,h);
  end Speel;

  procedure Speel ( c : in Circuit;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector ) is
  begin
    Speel(c.xps,c.cff,c.cst,x,yd,c.fwd,c.bck,c.crs);
  end Speel;

  procedure Speel ( c : in Circuit;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                    pwt : in Standard_Complex_VecVecs.VecVec ) is
  begin
    Speel(c.xps,c.idx,c.fac,c.cff,c.cst,x,yd,c.fwd,c.bck,c.crs,pwt);
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    cff : in Standard_Complex_Vectors.Vector;
                    cst : in Standard_Complex_Numbers.Complex_Number;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                    fwd : in Standard_Complex_Vectors.Link_to_Vector;
                    bck : in Standard_Complex_Vectors.Link_to_Vector;
                    crs : in Standard_Complex_Vectors.Link_to_Vector ) is

    idk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    kcff : Standard_Complex_Numbers.Complex_Number;

    use Standard_Complex_Numbers,Standard_Integer_Vectors;

  begin
    yd(0) := cst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          idx1 := idk(1); kcff := cff(k);
          yd(0) := yd(0) + kcff*x(idx1);
          yd(idx1) := yd(idx1) + kcff;
        else
          Forward_Backward_Cross(idk.all,x,fwd,bck,crs);
          kcff := cff(k); yd(0) := yd(0) + kcff*fwd(idk'last-1); 
          if idk'last = 2 then
            idx1 := idk(1); idx2 := idk(2);
            yd(idx2) := yd(idx2) + kcff*x(idx1);
            yd(idx1) := yd(idx1) + kcff*x(idx2);
          else -- idk'last > 2
            idx1 := idk(1);
            yd(idx1) := yd(idx1) + kcff*bck(idk'last-2);
            for j in idk'first+1..idk'last-1 loop
              idx1 := idk(j);
              yd(idx1) := yd(idx1) + kcff*crs(j-1);
            end loop;
            idx1 := idk(idk'last);
            yd(idx1) := yd(idx1) + kcff*fwd(idk'last-2);
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Indexed_Speel
              ( idx : in Standard_Integer_Vectors.Vector;
                cff : in Standard_Complex_Numbers.Complex_Number;
                x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                fwd : in Standard_Complex_Vectors.Link_to_Vector;
                bck : in Standard_Complex_Vectors.Link_to_Vector;
                crs : in Standard_Complex_Vectors.Link_to_Vector;
                h : in out Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;

    sz : constant integer32 := idx'last;
    idx1,idx2,idx3,idx4,idx5,idx6 : integer32;
    acc : Complex_Number;

  begin
    Forward_Backward_Cross(idx,x,fwd,bck,crs);
    yd(0) := yd(0) + cff*fwd(sz-1); 
    if sz = 2 then
      idx1 := idx(1); idx2 := idx(2);
      yd(idx2) := yd(idx2) + cff*x(idx1);
      yd(idx1) := yd(idx1) + cff*x(idx2);
      h(idx1,idx2) := h(idx1,idx2) + cff;
    else -- sz > 2
      idx1 := idx(1);
      yd(idx1) := yd(idx1) + cff*bck(sz-2);
      for j in idx'first+1..sz-1 loop
        idx1 := idx(j); yd(idx1) := yd(idx1) + cff*crs(j-1);
      end loop;
      idx1 := idx(sz); yd(idx1) := yd(idx1) + cff*fwd(sz-2);
      if sz = 3 then
        idx1 := idx(1); idx2 := idx(2); idx3 := idx(3);
        h(idx1,idx2) := h(idx1,idx2) + cff*x(idx3);
        h(idx1,idx3) := h(idx1,idx3) + cff*x(idx2);
        h(idx2,idx3) := h(idx2,idx3) + cff*x(idx1);
      else -- sz > 3
        idx5 := idx(sz-1); idx6 := idx(sz); idx3 := sz-3;
       -- last element is copy of fwd(sz-3), multiplied with cff
        h(idx5,idx6) := h(idx5,idx6) + cff*fwd(idx3);
       -- first element is copy of bck(sz-3), multiplied with cff
        idx1 := idx(1); idx2 := idx(2);
        h(idx1,idx2) := h(idx1,idx2) + cff*bck(idx3);
        if sz = 4 then -- special case for all rows
          idx3 := idx(3); idx4 := idx(4);
          acc := cff*x(idx2);
          h(idx1,idx3) := h(idx1,idx3) + acc*x(idx6);
          h(idx1,idx4) := h(idx1,idx4) + acc*x(idx5);
          acc := cff*x(idx1);
          h(idx2,idx3) := h(idx2,idx3) + acc*x(idx6);
          h(idx2,idx4) := h(idx2,idx4) + acc*x(idx5);
        else -- sz > 4
         -- first row is special, starts with x(idx(2)) after diagonal
          idx3 := idx(3); idx4 := sz-4;
          acc := cff*x(idx2);
          h(idx1,idx3) := h(idx1,idx3) + acc*bck(idx4);
          for k in 4..sz-2 loop
            idx4 := idx(k-1); idx5 := idx(k); idx6 := sz-k-1;
            acc := acc*x(idx4);
            h(idx1,idx5) := h(idx1,idx5) + acc*bck(idx6);
          end loop;
          idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
          acc := acc*x(idx4);
          h(idx1,idx5) := h(idx1,idx5) + acc*x(idx6);
          h(idx1,idx6) := h(idx1,idx6) + acc*x(idx5);
         -- second row is special, starts with x(idx(1)) after diagonal
          acc := cff*x(idx1); idx4 := sz-4;
          h(idx2,idx3) := h(idx2,idx3) + acc*bck(idx4);
          for k in 4..sz-2 loop
            idx4 := idx(k-1); idx5 := idx(k); idx6 := sz-k-1;
            acc := acc*x(idx4);
            h(idx2,idx5) := h(idx2,idx5) + acc*bck(idx6);
          end loop;
          idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
          acc := acc*x(idx4);
          h(idx2,idx5) := h(idx2,idx5) + acc*x(idx6);
          h(idx2,idx6) := h(idx2,idx6) + acc*x(idx5);
         -- the row with index sz-2 has a general formula
          idx3 := sz-4; acc := cff*fwd(idx3);
          h(idx4,idx5) := h(idx4,idx5) + acc*x(idx6);
          h(idx4,idx6) := h(idx4,idx6) + acc*x(idx5);
          for rw in 3..sz-3 loop  -- row rw starts with fwd(rw-2)
            idx1 := idx(rw); idx2 := idx(rw+1);
            acc := cff*fwd(rw-2);
            h(idx1,idx2) := h(idx1,idx2) + acc*bck(sz-rw-2);
            for k in rw+2..sz-2 loop
              idx2 := idx(k-1); idx3 := idx(k); idx4 := sz-k-1;
              acc := acc*x(idx2);
              h(idx1,idx3) := h(idx1,idx3) + acc*bck(idx4);
            end loop;
            idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
            acc := acc*x(idx4);
            h(idx1,idx5) := h(idx1,idx5) + acc*x(idx6);
            h(idx1,idx6) := h(idx1,idx6) + acc*x(idx5);
          end loop;
        end if;
      end if;
    end if;
  end Indexed_Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    cff : in Standard_Complex_Vectors.Vector;
                    cst : in Standard_Complex_Numbers.Complex_Number;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                    fwd : in Standard_Complex_Vectors.Link_to_Vector;
                    bck : in Standard_Complex_Vectors.Link_to_Vector;
                    crs : in Standard_Complex_Vectors.Link_to_Vector;
                    h : out Standard_Complex_Matrices.Matrix ) is

    dim : constant integer32 := x'last;
    idk : Standard_Integer_Vectors.Link_to_Vector;
    idx1 : integer32;
    kcff : Standard_Complex_Numbers.Complex_Number;

    use Standard_Complex_Numbers,Standard_Integer_Vectors;

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        h(i,j) := Create(0.0);
      end loop;
    end loop;
    yd(0) := cst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        kcff := cff(k);
        if idk'last = 1 then
          idx1 := idk(1);
          yd(0) := yd(0) + kcff*x(idx1);
          yd(idx1) := yd(idx1) + kcff;
        else
          Indexed_Speel(idk.all,kcff,x,yd,fwd,bck,crs,h);
        end if;
      end if;
    end loop;
    for i in 2..dim loop
      for j in 1..(i-1) loop
        h(i,j) := h(j,i);
      end loop;
    end loop;
  end Speel;

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    cff : in Standard_Complex_Vectors.Vector;
                    cst : in Standard_Complex_Numbers.Complex_Number;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                    fwd : in Standard_Complex_Vectors.Link_to_Vector;
                    bck : in Standard_Complex_Vectors.Link_to_Vector;
                    crs : in Standard_Complex_Vectors.Link_to_Vector;
                    pwt : in Standard_Complex_VecVecs.VecVec ) is

    idk,fck,xpk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    kcff,acc,wrk : Standard_Complex_Numbers.Complex_Number;
    factor : double_float;

    use Standard_Complex_Numbers,Standard_Integer_Vectors;

  begin
    yd(0) := cst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        kcff := cff(k);
        fck := fac(k);
        idx1 := idk(1);
        if fck = null then
          if idk'last = 1 then
            yd(0) := yd(0) + kcff*x(idx1);
            yd(idx1) := yd(idx1) + kcff;
          else
            Forward_Backward_Cross(idk.all,x,fwd,bck,crs);
            yd(0) := yd(0) + kcff*fwd(idk'last-1); 
            if idk'last = 2 then
              idx2 := idk(2);
              yd(idx2) := yd(idx2) + kcff*x(idx1);
              yd(idx1) := yd(idx1) + kcff*x(idx2);
            else -- idk'last > 2
              yd(idx1) := yd(idx1) + kcff*bck(idk'last-2);
              for j in idk'first+1..idk'last-1 loop
                idx2 := idk(j);
                yd(idx2) := yd(idx2) + kcff*crs(j-1);
              end loop;
              idx2 := idk(idk'last);
              yd(idx2) := yd(idx2) + kcff*fwd(idk'last-2);
            end if;
          end if;
        else
          xpk := xps(k);
          if idk'last = 1 then
            Multiply_Factor(xpk,fck,x,kcff,pwt,acc);
            yd(0) := yd(0) + acc*x(idx1);
            factor := Create(xpk(idx1));
            yd(idx1) := yd(idx1) + factor*acc;
          else
            Forward_Backward_Cross(idk.all,x,fwd,bck,crs);
            Multiply_Factor(xpk,fck,x,kcff,pwt,acc);
            yd(0) := yd(0) + acc*fwd(idk'last-1); 
            if idk'last = 2 then
              idx2 := idk(2);
              wrk := acc*x(idx1); factor := Create(xpk(idx2));
              wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
              wrk := acc*x(idx2); factor := Create(xpk(idx1));
              wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
            else -- idk'last > 2
              wrk := acc*bck(idk'last-2);
              factor := Create(xpk(idx1)); wrk := factor*wrk;
              yd(idx1) := yd(idx1) + wrk;
              for j in idk'first+1..idk'last-1 loop
                wrk := acc*crs(j-1); idx2 := idk(j);
                factor := Create(xpk(idx2)); wrk := factor*wrk;
                yd(idx2) := yd(idx2) + wrk;
              end loop;
              wrk := acc*fwd(idk'last-2); idx2 := idk(idk'last);
              factor := Create(xpk(idx2)); wrk := factor*wrk;
              yd(idx2) := yd(idx2) + wrk;
            end if;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

-- AUXILIARY PROCEDURES :

  procedure Forward ( x : in Standard_Complex_Vectors.Link_to_Vector;
                      f : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
  end Forward;

  procedure Forward_Backward
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
    b(b'first) := x(x'last)*x(x'last-1);
    for k in 2..x'last-2 loop
      b(k) := b(k-1)*x(x'last-k);
    end loop;
  end Forward_Backward;

  procedure Fused_Forward_Backward
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    b(b'first) := x(x'last)*x(x'last-1);
    for k in 2..x'last-2 loop
      f(k) := f(k-1)*x(k+1);
      b(k) := b(k-1)*x(x'last-k);
    end loop;
    if f'last > 1
     then f(f'last) := f(f'last-1)*x(x'last);
    end if;
  end Fused_Forward_Backward;

  procedure Forward_Backward_Cross
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
    if x'last > 2 then
      b(b'first) := x(x'last)*x(x'last-1);
      for k in 2..x'last-2 loop
        b(k) := b(k-1)*x(x'last-k);
      end loop;
      if x'last = 3 then
        c(1) := x(1)*x(3);
      else
        c(1) := x(1)*b(x'last-3);
        for k in 2..x'last-3 loop
          c(k) := f(k-1)*b(x'last-2-k);
        end loop;
        c(x'last-2) := x(x'last)*f(x'last-3);
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(1) := x(idx(1))*x(idx(2));
    for k in 2..idx'last-1 loop
      f(k) := f(k-1)*x(idx(k+1));
    end loop;
    if idx'last > 2 then
      b(b'first) := x(idx(idx'last))*x(idx(idx'last-1));
      for k in 2..idx'last-2 loop
        b(k) := b(k-1)*x(idx(idx'last-k));
      end loop;
      if idx'last = 3 then
        c(1) := x(idx(1))*x(idx(3));
      else
        c(1) := x(idx(1))*b(idx'last-3);
        for k in 2..idx'last-3 loop
          c(k) := f(k-1)*b(idx'last-2-k);
        end loop;
        c(idx'last-2) := x(idx(idx'last))*f(idx'last-3);
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Fused_Forward_Backward_Cross
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

    firstend,lastend,plusidx,minidx : integer32;

  begin
    if x'last >= 8 then
      if x'last mod 2 = 0 then
        lastend := x'last-4;
        firstend := lastend/2;
        f(f'first) := x(x'first)*x(x'first+1);
        b(b'first) := x(x'last)*x(x'last-1);
        for k in 2..firstend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
        end loop;
        minidx := firstend+1; plusidx := minidx+1;
        for k in firstend+1..lastend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
          c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          plusidx := plusidx + 1;
          c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          minidx := minidx - 1;
        end loop;
      else
        lastend := x'last-4;
        firstend := (x'last-3)/2;
        f(f'first) := x(x'first)*x(x'first+1);
        b(b'first) := x(x'last)*x(x'last-1);
        for k in 2..firstend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
        end loop;
        plusidx := firstend+1;
        c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
        minidx := plusidx;
        for k in firstend+1..lastend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
          plusidx := plusidx + 1;
          c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          minidx := minidx - 1;
          c(minidx) := f(minidx-1)*b(x'last-2-minidx);
        end loop;
      end if;
      plusidx := lastend+1;
      f(plusidx) := f(plusidx-1)*x(plusidx+1);
      b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
      plusidx := plusidx+1;
      f(plusidx) := f(plusidx-1)*x(plusidx+1);
      b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
      f(f'last) := f(f'last-1)*x(x'last);
      c(1) := x(1)*b(x'last-3);
      c(x'last-2) := x(x'last)*f(x'last-3);
    else
      f(f'first) := x(x'first)*x(x'first+1);
      b(b'first) := x(x'last)*x(x'last-1);
      for k in 2..x'last-2 loop
        f(k) := f(k-1)*x(k+1);
        b(k) := b(k-1)*x(x'last-k);
      end loop;
      if f'last > 1
       then f(f'last) := f(f'last-1)*x(x'last);
      end if;
      if x'last > 3 then
        c(1) := x(1)*b(x'last-3);
        for k in 2..x'last-3 loop
          c(k) := f(k-1)*b(x'last-2-k);
        end loop;
        c(x'last-2) := x(x'last)*f(x'last-3);
      elsif x'last = 3 then
        c(1) := x(1)*x(3);
      end if;
    end if;
  end Fused_Forward_Backward_Cross;

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(mxe'range);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    for k in mxe'range loop
      if mxe(k) > 1 then
        res(k) := new Standard_Complex_Vectors.Vector'(1..mxe(k)-1 => zero);
      end if;
    end loop;
    return res;
  end Allocate;

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                pwt : in Standard_Complex_VecVecs.VecVec ) is

    lnk : Standard_Complex_Vectors.Link_to_Vector;

    use Standard_Complex_Numbers;

  begin
    for k in x'range loop
      if mxe(k) > 1 then
        lnk := pwt(k);
        lnk(1) := x(k)*x(k);
        for i in 2..mxe(k)-1 loop
          lnk(i) := lnk(i-1)*x(k);
        end loop;
      end if;
    end loop;
  end Power_Table;

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                cff : in Standard_Complex_Numbers.Complex_Number;
                pwt : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Numbers.Complex_Number ) is

    pwx : Standard_Complex_Vectors.Link_to_Vector;
    idx,powidx : integer32;

    use Standard_Complex_Numbers;

  begin
    idx := fac(fac'first); powidx := xps(idx);
    if powidx = 2 then
      res := cff*x(idx);
    else
      pwx := pwt(idx);
      res := cff*pwx(powidx-2);
    end if;
    for k in fac'first+1..fac'last loop
      idx := fac(k); powidx := xps(idx);
      if powidx = 2 then
        res := res*x(idx);
      else
        pwx := pwt(idx);
        res := res*pwx(powidx-2);
      end if;
    end loop;
  end Multiply_Factor;

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit ) is
  begin
    Standard_Complex_Vectors.Clear(c.fwd);
    Standard_Complex_Vectors.Clear(c.bck);
    Standard_Complex_Vectors.Clear(c.crs);
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

  procedure Clear ( s : in out System ) is
  begin
    Clear(s.crc);
    Standard_Complex_Vectors.Clear(s.yd);
    Standard_Complex_VecVecs.Clear(s.pwt);
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

end Standard_Complex_Circuits;
