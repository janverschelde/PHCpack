with Quad_Double_Numbers;                use Quad_Double_Numbers;

package body QuadDobl_Vector_Splitters is

  procedure Split ( x : in Complex_Number;
                    rehihi,imhihi,relohi,imlohi : out double_float;
                    rehilo,imhilo,relolo,imlolo : out double_float ) is

    nbr : quad_double;

  begin
    nbr := REAL_PART(x);
    rehihi := hihi_part(nbr); relohi := lohi_part(nbr);
    rehilo := hilo_part(nbr); relolo := lolo_part(nbr);
    nbr := IMAG_PART(x);
    imhihi := hihi_part(nbr); imlohi := lohi_part(nbr);
    imhilo := hilo_part(nbr); imlolo := lolo_part(nbr);
  end Split;

  procedure Merge ( x : out Complex_Number;
                    rehihi,imhihi,relohi,imlohi : in double_float;
                    rehilo,imhilo,relolo,imlolo : in double_float ) is

    real_part : constant quad_double
              := Quad_Double_Numbers.create(rehihi,relohi,rehilo,relolo);
    imag_part : constant quad_double
              := Quad_Double_Numbers.create(imhihi,imlohi,imhilo,imlolo);

  begin
    x := QuadDobl_Complex_Numbers.Create(real_part,imag_part);
  end Merge;

  procedure Split ( x : in QuadDobl_Complex_Vectors.Vector;
                    xrhh,xihh : out Standard_Floating_Vectors.Vector;
                    xrlh,xilh : out Standard_Floating_Vectors.Vector;
                    xrhl,xihl : out Standard_Floating_Vectors.Vector;
                    xrll,xill : out Standard_Floating_Vectors.Vector ) is
  begin
    for k in x'range loop
      Split(x(k),xrhh(k),xihh(k),xrlh(k),xilh(k),
                 xrhl(k),xihl(k),xrll(k),xill(k));
    end loop;
  end Split;

  procedure Merge ( x : out QuadDobl_Complex_Vectors.Vector;
                    xrhh,xihh : in Standard_Floating_Vectors.Vector;
                    xrlh,xilh : in Standard_Floating_Vectors.Vector;
                    xrhl,xihl : in Standard_Floating_Vectors.Vector;
                    xrll,xill : in Standard_Floating_Vectors.Vector ) is
  begin
    for k in x'range loop
      Merge(x(k),xrhh(k),xihh(k),xrlh(k),xilh(k),
                 xrhl(k),xihl(k),xrll(k),xill(k));
    end loop;
  end Merge;

-- 2-VECTOR REPRESENTATIONS :

  procedure Two_Split
              ( x : in QuadDobl_Complex_Vectors.Vector;
                xr,xi : out Standard_Floating_Vectors.Vector ) is

    xrehihi,ximhihi,xrelohi,ximlohi : double_float;
    xrehilo,ximhilo,xrelolo,ximlolo : double_float;
    idx : integer32 := xr'first;

  begin
    for k in x'range loop
      Split(x(k),xrehihi,ximhihi,xrelohi,ximlohi,
                 xrehilo,ximhilo,xrelolo,ximlolo);
      xr(idx)   := xrehihi;  xi(idx)   := ximhihi;
      xr(idx+1) := xrelohi;  xi(idx+1) := ximlohi;
      xr(idx+2) := xrehilo;  xi(idx+2) := ximhilo;
      xr(idx+3) := xrelolo;  xi(idx+3) := ximlolo;
      idx := idx + 4;
    end loop;
  end Two_Split;

  procedure Two_Merge
              ( x : out QuadDobl_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Vector ) is

    idx : integer32 := xr'first;

  begin
    for k in x'range loop
      Merge(x(k),xr(idx),xi(idx),xr(idx+1),xi(idx+1),
            xr(idx+2),xi(idx+2),xr(idx+3),xi(idx+3));
      idx := idx + 4;
    end loop;
  end Two_Merge;

-- SPLITTERS AND MERGERS WITH AND WITHOUT MEMORY ALLOCATIONS :

  procedure Split_Complex
              ( x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                xr,xi : out Standard_Floating_Vectors.Link_to_Vector ) is

    use QuadDobl_Complex_Vectors;

  begin
    if x /= null then
      if x'first = 1 then
        declare
          dim : constant integer32 := x'last;
          vr,vi : Standard_Floating_Vectors.Vector(1..4*dim);
        begin
          Two_Split(x.all,vr,vi);
          xr := new Standard_Floating_Vectors.Vector'(vr);
          xi := new Standard_Floating_Vectors.Vector'(vi);
        end;
      else -- x'first = 0
        declare
          deg : constant integer32 := x'last;
          dim : constant integer32 := 4*(deg+1)-1;
          vr,vi : Standard_Floating_Vectors.Vector(0..dim);
        begin
          Two_Split(x.all,vr,vi);
          xr := new Standard_Floating_Vectors.Vector'(vr);
          xi := new Standard_Floating_Vectors.Vector'(vi);
        end;
      end if;
    end if;
  end Split_Complex;

  procedure Split_Complex
              ( x : in QuadDobl_Complex_VecVecs.VecVec;
                xr,xi : out Standard_Floating_VecVecs.VecVec ) is

    use QuadDobl_Complex_Vectors;

  begin
    for k in x'range loop
      if x(k) /= null
       then Split_Complex(x(k),xr(k),xi(k));
      end if;
    end loop;
  end Split_Complex;

  procedure Split_Complex
              ( x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                xr,xi : out Standard_Floating_VecVecs.Link_to_VecVec ) is

    use QuadDobl_Complex_Vectors,QuadDobl_Complex_VecVecs;

  begin
    if x /= null then
      declare
        vr,vi : Standard_Floating_VecVecs.VecVec(x'range);
      begin
        for k in x'range loop
          if x(k) /= null
           then Split_Complex(x(k),vr(k),vi(k));
          end if;
        end loop;
        xr := new Standard_Floating_VecVecs.VecVec'(vr);
        xi := new Standard_Floating_VecVecs.VecVec'(vi);
      end;
    end if;
  end Split_Complex;

-- MEMORY ALLOCATORS :

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return QuadDobl_Complex_Vectors.Link_to_Vector is

    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    cff : constant QuadDobl_Complex_Vectors.Vector(0..deg) := (0..deg => zero);
    res : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector'(cff);

  begin
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..dim);

  begin
    for k in res'range loop
      res(k) := Allocate_Complex_Coefficients(deg);
    end loop;
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return QuadDobl_Complex_VecVecs.Link_to_VecVec is

    res : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    cff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
        := Allocate_Complex_Coefficients(dim,deg);

  begin
    res := new QuadDobl_Complex_VecVecs.VecVec'(cff);
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in res'range loop
      declare
        v : constant QuadDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new QuadDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return QuadDobl_Complex_VecVecs.Link_to_VecVec is

    res : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    vv : QuadDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in vv'range loop
      declare
        v : constant QuadDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        vv(i) := new QuadDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    res := new QuadDobl_Complex_VecVecs.VecVec'(vv);
    return res;
  end Allocate;

-- SPLITTERS AND MERGERS ON ALLOCATED VECTORS :

  procedure Complex_Parts
              ( x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector ) is

    idx : integer32 := xr'first;
    nbr : quad_double;

  begin
    for k in x'range loop
      nbr := QuadDobl_Complex_Numbers.REAL_PART(x(k));
      xr(idx) := hihi_part(nbr);
      xr(idx+1) := lohi_part(nbr);
      xr(idx+2) := hilo_part(nbr);
      xr(idx+3) := lolo_part(nbr);
      nbr := QuadDobl_Complex_Numbers.IMAG_PART(x(k));
      xi(idx) := hihi_part(nbr);
      xi(idx+1) := lohi_part(nbr);
      xi(idx+2) := hilo_part(nbr);
      xi(idx+3) := lolo_part(nbr);
      idx := idx + 4;
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( x : in QuadDobl_Complex_VecVecs.VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(x(k),xr(k),xi(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(x(k),xr(k),xi(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( deg : in integer32;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector ) is

    idx : integer32 := xr'first;
    nbr : quad_double;

  begin
    for k in x'first..deg loop
      nbr := QuadDobl_Complex_Numbers.REAL_PART(x(k));
      xr(idx) := hihi_part(nbr);
      xr(idx+1) := lohi_part(nbr);
      xr(idx+2) := hilo_part(nbr);
      xr(idx+3) := lolo_part(nbr);
      nbr := QuadDobl_Complex_Numbers.IMAG_PART(x(k));
      xi(idx) := hihi_part(nbr);
      xi(idx+1) := lohi_part(nbr);
      xi(idx+2) := hilo_part(nbr);
      xi(idx+3) := lolo_part(nbr);
      idx := idx + 4;
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( deg : in integer32;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(deg,x(k),xr(k),xi(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( deg : in integer32;
                x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(deg,x(k),xr(k),xi(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Merge
              ( xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := xr'first;
    rqd,iqd : quad_double;

  begin
    for k in cvx'range loop
      rqd := Create(xr(idx),xr(idx+1),xr(idx+2),xr(idx+3));
      iqd := Create(xi(idx),xi(idx+1),xi(idx+2),xi(idx+3));
      cvx(k) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
      idx := idx + 4;
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.Link_to_VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(xr(k),xi(k),cvx(k));
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(xr(k),xi(k),cvx(k));
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( deg : in integer32; 
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := xr'first;
    rqd,iqd : quad_double;

  begin
    for k in cvx'first..deg loop
      rqd := Create(xr(idx),xr(idx+1),xr(idx+2),xr(idx+3));
      iqd := Create(xi(idx),xi(idx+1),xi(idx+2),xi(idx+3));
      cvx(k) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
      idx := idx + 4;
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( deg : in integer32;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.Link_to_VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(deg,xr(k),xi(k),cvx(k));
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( deg : in integer32;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(deg,xr(k),xi(k),cvx(k));
    end loop;
  end Complex_Merge;

-- VECTOR COMPUTATIONS :

  procedure Add ( x,y : in Standard_Floating_Vectors.Vector;
                  z : out Standard_Floating_Vectors.Vector ) is

    i,j,k : integer32;
    s,t : double_float;
    u,v : double_float; -- double-length accumulator
    c0,c1,c2,c3,s0,s1,s2,s3 : double_float; -- for the renorm4
    za,zb : boolean;

  begin
    z(0) := 0.0; z(1) := 0.0; z(2) := 0.0; z(3) := 0.0;
    i := 0; j := 0; k := 0;
    if abs(x(i)) > abs(y(j))
     then u := x(i); i := i+1;
     else u := y(j); j := j+1;
    end if;
    if abs(x(i)) > abs(y(j))
     then v := x(i); i := i+1;  
     else v := y(j); j := j+1;
    end if;
   -- Double_Double_Basics.quick_two_sum(u,v,u,v);
    s := u + v; t := v - (s - u); u := s; v := t; 
    while k < 4 loop
      if (i >= 4 and j >= 4) then
        z(k) := u;
        if k < 3
         then k := k+1; z(k) := v;
        end if;
        exit;
      end if;
      if i >= 4 then
        t := y(j); j := j+1;
      elsif j >= 4  then
        t := x(i); i := i+1;
      elsif abs(x(i)) > abs(y(j)) then
        t := x(i); i := i+1;
      else 
        t := y(j); j := j+1;
      end if;
     -- Quad_Double_Renormalizations.quick_three_accum(u,v,s,t);
     -- Double_Double_Basics.two_sum(v,t,s,v);
      s := v + t; c0 := s - v; v := (v - (s - c0)) + (t - c0);
     -- Double_Double_Basics.two_sum(u,s,s,u);
      c0 := s; s := u + s; c1 := s - u; u := (u - (s - c1)) + (c0 - c1);
      za := (u /= 0.0);
      zb := (v /= 0.0);
      if za and zb then
        null;
      else
        if not zb
         then v := u; u := s;
         else u := s;
        end if;
        s := 0.0;
      end if;
      if s /= 0.0
       then z(k) := s; k := k+1;
      end if;
    end loop;
    for k in i..3 loop                    -- add the rest
      z(3) := z(3) + x(k);
    end loop;
    for k in j..3 loop
      z(3) := z(3) + y(k);
    end loop;
   -- Quad_Double_Renormalizations.renorm4(z(0),z(1),z(2),z(3));
    c0 := z(0); c1 := z(1); c2 := z(2); c3 := z(3);
   -- Double_Double_Basics.quick_two_sum(c2,c3,s0,c3);
    s0 := c2 + c3; c3 := c3 - (s0 - c2);
   -- Double_Double_Basics.quick_two_sum(c1,s0,s0,c2);
    s := c1 + s0; c2 := s0 - (s - c1); s0 := s;
   -- Double_Double_Basics.quick_two_sum(c0,s0,c0,c1);
    s := c0 + s0; c1 := s0 - (s - c0); c0 := s;
    s0 := c0; s1 := c1;
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,c2,s1,s2);
      s := s1 + c2; s2 := c2 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,c3,s2,s3);
        s := s2 + c3; s3 := c3 - (s - s2); s2 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(s0,c2,s0,s1);
      s := s0 + c2; s1 := c2 - (s - s0); s0 := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s0,c3,s0,s1);
        s := s0 + c3; s1 := c3 - (s - s0); s0 := s;
      end if;
    end if;
   -- c0 := s0; c1 := s1; c2 := s2; c3 := s3;
    z(0) := s0; z(1) := s1; z(2) := s2; z(3) := s3;
  end Add;

  procedure Add ( offset : in integer32;
                  x,y,z : in Standard_Floating_Vectors.Link_to_Vector ) is

    i,j,k : integer32;
    s,t : double_float;
    u,v : double_float; -- double-length accumulator
    c0,c1,c2,c3,s0,s1,s2,s3 : double_float; -- for the renorm4
    za,zb : boolean;

  begin
    z(offset) := 0.0; z(offset+1) := 0.0;
    z(offset+2) := 0.0; z(offset+3) := 0.0;
    i := offset; j := offset; k := offset;
    if abs(x(i)) > abs(y(j))
     then u := x(i); i := i+1;
     else u := y(j); j := j+1;
    end if;
    if abs(x(i)) > abs(y(j))
     then v := x(i); i := i+1;  
     else v := y(j); j := j+1;
    end if;
   -- Double_Double_Basics.quick_two_sum(u,v,u,v);
    s := u + v; t := v - (s - u); u := s; v := t; 
    while k < offset+4 loop
      if (i >= offset+4 and j >= offset+4) then
        z(k) := u;
        if k < offset+3
         then k := k+1; z(k) := v;
        end if;
        exit;
      end if;
      if i >= offset+4 then
        t := y(j); j := j+1;
      elsif j >= offset+4  then
        t := x(i); i := i+1;
      elsif abs(x(i)) > abs(y(j)) then
        t := x(i); i := i+1;
      else 
        t := y(j); j := j+1;
      end if;
     -- Quad_Double_Renormalizations.quick_three_accum(u,v,s,t);
     -- Double_Double_Basics.two_sum(v,t,s,v);
      s := v + t; c0 := s - v; v := (v - (s - c0)) + (t - c0);
     -- Double_Double_Basics.two_sum(u,s,s,u);
      c0 := s; s := u + s; c1 := s - u; u := (u - (s - c1)) + (c0 - c1);
      za := (u /= 0.0);
      zb := (v /= 0.0);
      if za and zb then
        null;
      else
        if not zb
         then v := u; u := s;
         else u := s;
        end if;
        s := 0.0;
      end if;
      if s /= 0.0
       then z(k) := s; k := k+1;
      end if;
    end loop;
    for k in i..offset+3 loop                    -- add the rest
      z(3) := z(3) + x(k);
    end loop;
    for k in j..offset+3 loop
      z(3) := z(3) + y(k);
    end loop;
   -- Quad_Double_Renormalizations.renorm4(z(0),z(1),z(2),z(3));
    c0 := z(offset); c1 := z(offset+1);
    c2 := z(offset+2); c3 := z(offset+3);
   -- Double_Double_Basics.quick_two_sum(c2,c3,s0,c3);
    s0 := c2 + c3; c3 := c3 - (s0 - c2);
   -- Double_Double_Basics.quick_two_sum(c1,s0,s0,c2);
    s := c1 + s0; c2 := s0 - (s - c1); s0 := s;
   -- Double_Double_Basics.quick_two_sum(c0,s0,c0,c1);
    s := c0 + s0; c1 := s0 - (s - c0); c0 := s;
    s0 := c0; s1 := c1;
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,c2,s1,s2);
      s := s1 + c2; s2 := c2 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,c3,s2,s3);
        s := s2 + c3; s3 := c3 - (s - s2); s2 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(s0,c2,s0,s1);
      s := s0 + c2; s1 := c2 - (s - s0); s0 := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s0,c3,s0,s1);
        s := s0 + c3; s1 := c3 - (s - s0); s0 := s;
      end if;
    end if;
   -- c0 := s0; c1 := s1; c2 := s2; c3 := s3;
    z(offset) := s0; z(offset+1) := s1;
    z(offset+2) := s2; z(offset+3) := s3;
  end Add;

  procedure Add ( x,y,z : in Standard_Floating_Vectors.Link_to_Vector ) is

    i,j,k : integer32;
    s,t : double_float;
    u,v : double_float; -- double-length accumulator
    c0,c1,c2,c3,s0,s1,s2,s3 : double_float; -- for the renorm4
    za,zb : boolean;

  begin
    z(0) := 0.0; z(1) := 0.0; z(2) := 0.0; z(3) := 0.0;
    i := 0; j := 0; k := 0;
    if abs(x(i)) > abs(y(j))
     then u := x(i); i := i+1;
     else u := y(j); j := j+1;
    end if;
    if abs(x(i)) > abs(y(j))
     then v := x(i); i := i+1;  
     else v := y(j); j := j+1;
    end if;
   -- Double_Double_Basics.quick_two_sum(u,v,u,v);
    s := u + v; t := v - (s - u); u := s; v := t; 
    while k < 4 loop
      if (i >= 4 and j >= 4) then
        z(k) := u;
        if k < 3
         then k := k+1; z(k) := v;
        end if;
        exit;
      end if;
      if i >= 4 then
        t := y(j); j := j+1;
      elsif j >= 4  then
        t := x(i); i := i+1;
      elsif abs(x(i)) > abs(y(j)) then
        t := x(i); i := i+1;
      else 
        t := y(j); j := j+1;
      end if;
     -- Quad_Double_Renormalizations.quick_three_accum(u,v,s,t);
     -- Double_Double_Basics.two_sum(v,t,s,v);
      s := v + t; c0 := s - v; v := (v - (s - c0)) + (t - c0);
     -- Double_Double_Basics.two_sum(u,s,s,u);
      c0 := s; s := u + s; c1 := s - u; u := (u - (s - c1)) + (c0 - c1);
      za := (u /= 0.0);
      zb := (v /= 0.0);
      if za and zb then
        null;
      else
        if not zb
         then v := u; u := s;
         else u := s;
        end if;
        s := 0.0;
      end if;
      if s /= 0.0
       then z(k) := s; k := k+1;
      end if;
    end loop;
    for k in i..3 loop                    -- add the rest
      z(3) := z(3) + x(k);
    end loop;
    for k in j..3 loop
      z(3) := z(3) + y(k);
    end loop;
   -- Quad_Double_Renormalizations.renorm4(z(0),z(1),z(2),z(3));
    c0 := z(0); c1 := z(1); c2 := z(2); c3 := z(3);
   -- Double_Double_Basics.quick_two_sum(c2,c3,s0,c3);
    s0 := c2 + c3; c3 := c3 - (s0 - c2);
   -- Double_Double_Basics.quick_two_sum(c1,s0,s0,c2);
    s := c1 + s0; c2 := s0 - (s - c1); s0 := s;
   -- Double_Double_Basics.quick_two_sum(c0,s0,c0,c1);
    s := c0 + s0; c1 := s0 - (s - c0); c0 := s;
    s0 := c0; s1 := c1;
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,c2,s1,s2);
      s := s1 + c2; s2 := c2 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,c3,s2,s3);
        s := s2 + c3; s3 := c3 - (s - s2); s2 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(s0,c2,s0,s1);
      s := s0 + c2; s1 := c2 - (s - s0); s0 := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s0,c3,s0,s1);
        s := s0 + c3; s1 := c3 - (s - s0); s0 := s;
      end if;
    end if;
   -- c0 := s0; c1 := s1; c2 := s2; c3 := s3;
    z(0) := s0; z(1) := s1; z(2) := s2; z(3) := s3;
  end Add;

  procedure Add ( zrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  zihh : in Standard_Floating_Vectors.Link_to_Vector;
                  zrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  zilh : in Standard_Floating_Vectors.Link_to_Vector;
                  zrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  zihl : in Standard_Floating_Vectors.Link_to_Vector;
                  zrll : in Standard_Floating_Vectors.Link_to_Vector;
                  zill : in Standard_Floating_Vectors.Link_to_Vector;
                  xrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  xihh : in Standard_Floating_Vectors.Link_to_Vector;
                  xrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  xilh : in Standard_Floating_Vectors.Link_to_Vector;
                  xrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  xihl : in Standard_Floating_Vectors.Link_to_Vector;
                  xrll : in Standard_Floating_Vectors.Link_to_Vector;
                  xill : in Standard_Floating_Vectors.Link_to_Vector;
                  yrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  yihh : in Standard_Floating_Vectors.Link_to_Vector;
                  yrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  yilh : in Standard_Floating_Vectors.Link_to_Vector;
                  yrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  yihl : in Standard_Floating_Vectors.Link_to_Vector;
                  yrll : in Standard_Floating_Vectors.Link_to_Vector;
                  yill : in Standard_Floating_Vectors.Link_to_Vector;
                  x,y,z : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in zrhh'range loop
      x(0) := xrhh(k); x(1) := xrlh(k); x(2) := xrhl(k); x(3) := xrll(k);
      y(0) := yrhh(k); y(1) := yrlh(k); y(2) := yrhl(k); y(3) := yrll(k);
      Add(x,y,z);
      zrhh(k) := z(0); zrlh(k) := z(1); zrhl(k) := z(2); zrll(k) := z(3); 
      x(0) := xihh(k); x(1) := xilh(k); x(2) := xihl(k); x(3) := xill(k);
      y(0) := yihh(k); y(1) := yilh(k); y(2) := yihl(k); y(3) := yill(k);
      Add(x,y,z);
      zihh(k) := z(0); zilh(k) := z(1); zihl(k) := z(2); zill(k) := z(3); 
    end loop;
  end Add;

  procedure Two_Add 
              ( zr : in Standard_Floating_Vectors.Link_to_Vector;
                zi : in Standard_Floating_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                yr : in Standard_Floating_Vectors.Link_to_Vector;
                yi : in Standard_Floating_Vectors.Link_to_Vector ) is

    dim : constant integer32 := zr'last/4;
    idx : integer32 := zr'first; -- allow for zero start index

  begin
    for k in zr'first..dim loop
      Add(idx,xr,yr,zr);
      Add(idx,xi,yi,zi);
      idx := idx + 4;
    end loop;
  end Two_Add;

  procedure Update ( offset : in integer32;
                     z,y,x : in Standard_Floating_Vectors.Link_to_Vector ) is

    i,j,k : integer32;
    s,t : double_float;
    u,v : double_float; -- double-length accumulator
    c0,c1,c2,c3,s0,s1,s2,s3 : double_float; -- for the renorm4
    za,zb : boolean;

  begin
    x(0) := z(offset); x(1) := z(offset+1);
    x(2) := z(offset+2); x(3) := z(offset+3);
    i := 0; j := offset; k := offset;
    if abs(x(i)) > abs(y(j))
     then u := x(i); i := i+1;
     else u := y(j); j := j+1;
    end if;
    if abs(x(i)) > abs(y(j))
     then v := x(i); i := i+1;  
     else v := y(j); j := j+1;
    end if;
   -- Double_Double_Basics.quick_two_sum(u,v,u,v);
    s := u + v; t := v - (s - u); u := s; v := t; 
    while k < offset+4 loop
      if (i >= 4 and j >= offset+4) then
        z(k) := u;
        if k < offset+3
         then k := k+1; z(k) := v;
        end if;
        exit;
      end if;
      if i >= 4 then
        t := y(j); j := j+1;
      elsif j >= offset+4  then
        t := x(i); i := i+1;
      elsif abs(x(i)) > abs(y(j)) then
        t := x(i); i := i+1;
      else 
        t := y(j); j := j+1;
      end if;
     -- Quad_Double_Renormalizations.quick_three_accum(u,v,s,t);
     -- Double_Double_Basics.two_sum(v,t,s,v);
      s := v + t; c0 := s - v; v := (v - (s - c0)) + (t - c0);
     -- Double_Double_Basics.two_sum(u,s,s,u);
      c0 := s; s := u + s; c1 := s - u; u := (u - (s - c1)) + (c0 - c1);
      za := (u /= 0.0);
      zb := (v /= 0.0);
      if za and zb then
        null;
      else
        if not zb
         then v := u; u := s;
         else u := s;
        end if;
        s := 0.0;
      end if;
      if s /= 0.0
       then z(k) := s; k := k+1;
      end if;
    end loop;
    for k in i..3 loop                    -- add the rest
      z(offset+3) := z(offset+3) + x(k);
    end loop;
    for k in j..offset+3 loop
      z(3) := z(3) + y(k);
    end loop;
   -- Quad_Double_Renormalizations.renorm4(z(0),z(1),z(2),z(3));
    c0 := z(offset); c1 := z(offset+1);
    c2 := z(offset+2); c3 := z(offset+3);
   -- Double_Double_Basics.quick_two_sum(c2,c3,s0,c3);
    s0 := c2 + c3; c3 := c3 - (s0 - c2);
   -- Double_Double_Basics.quick_two_sum(c1,s0,s0,c2);
    s := c1 + s0; c2 := s0 - (s - c1); s0 := s;
   -- Double_Double_Basics.quick_two_sum(c0,s0,c0,c1);
    s := c0 + s0; c1 := s0 - (s - c0); c0 := s;
    s0 := c0; s1 := c1;
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,c2,s1,s2);
      s := s1 + c2; s2 := c2 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,c3,s2,s3);
        s := s2 + c3; s3 := c3 - (s - s2); s2 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(s0,c2,s0,s1);
      s := s0 + c2; s1 := c2 - (s - s0); s0 := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,c3,s1,s2);
        s := s1 + c3; s2 := c3 - (s - s1); s1 := s;
      else
       -- Double_Double_Basics.quick_two_sum(s0,c3,s0,s1);
        s := s0 + c3; s1 := c3 - (s - s0); s0 := s;
      end if;
    end if;
   -- c0 := s0; c1 := s1; c2 := s2; c3 := s3;
    z(offset) := s0; z(offset+1) := s1;
    z(offset+2) := s2; z(offset+3) := s3;
  end Update;

  procedure Update
              ( zr : in Standard_Floating_Vectors.Link_to_Vector;
                zi : in Standard_Floating_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                wrk : in Standard_Floating_Vectors.Link_to_Vector ) is

    dim : constant integer32 := zr'last/4;
    idx : integer32 := zr'first; -- allow for zero start index

  begin
    for k in zr'first..dim loop
      Update(idx,zr,xr,wrk);
      Update(idx,zi,xi,wrk);
      idx := idx + 4;
    end loop;
  end Update;

  procedure Update_Product
              ( zrhh,zihh,zrlh,zilh : in out double_float;
                zrhl,zihl,zrll,zill : in out double_float;
                xrhh,xihh,xrlh,xilh : in double_float;
                xrhl,xihl,xrll,xill : in double_float;
                yrhh,yihh,yrlh,yilh : in double_float;
                yrhl,yihl,yrll,yill : in double_float;
                x,y,z : in Standard_Floating_Vectors.Link_to_Vector ) is

    QD_SPLITTER : constant double_float := 134217729.0; -- 2^27 + 1
    QD_SPLIT_THRESH : constant double_float := 6.69692879491417e+299; -- 2^996

    p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,q5,bb,s : double_float;
    p6,p7,p8,p9,q6,q7,q8,q9,r0,r1,t0,t1,s0,s1,s2,s3 : double_float;
    xrhh_hi,xrhh_lo,xihh_hi,xihh_lo : double_float;
    xrhl_hi,xrhl_lo,xihl_hi,xihl_lo : double_float;
    xrlh_hi,xrlh_lo,xilh_hi,xilh_lo : double_float;
    xrll_hi,xrll_lo,xill_hi,xill_lo : double_float;
    yrhh_hi,yrhh_lo,yihh_hi,yihh_lo : double_float;
    yrhl_hi,yrhl_lo,yihl_hi,yihl_lo : double_float;
    yrlh_hi,yrlh_lo,yilh_hi,yilh_lo : double_float;
    yrll_hi,yrll_lo,yill_hi,yill_lo : double_float;

  begin
   -- zre = xre*yre - xim*yim
   -- (1) compute xre*yre
   -- Double_Double_Basics.two_prod(xrhh,yrhh,p0,q0);
    p0 := xrhh*yrhh;
   -- Double_Double_Basics.split(xrhh,xrhh_hi,xrhh_lo);
    if ( xrhh > QD_SPLIT_THRESH or xrhh < -QD_SPLIT_THRESH ) then
      bb := xrhh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xrhh_hi := s - (s - bb);
      xrhh_lo := xrhh - xrhh_hi;
      xrhh_hi := xrhh_hi*268435456.0;  -- 2^28
      xrhh_lo := 268435456.0;          -- 2^28
    else
      s := QD_SPLITTER * xrhh;
      xrhh_hi := s - (s - xrhh);
      xrhh_lo := xrhh - xrhh_hi;
    end if;
   -- Double_Double_Basics.split(yrhh,yrhh_hi,yrhh_lo);
    if ( yrhh > QD_SPLIT_THRESH or yrhh < -QD_SPLIT_THRESH ) then
      bb := yrhh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yrhh_hi := s - (s - bb);
      yrhh_lo := yrhh - yrhh_hi;
      yrhh_hi := yrhh_hi*268435456.0;  -- 2^28
      yrhh_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yrhh;
      yrhh_hi := s - (s - yrhh);
      yrhh_lo := yrhh - yrhh_hi;
    end if;
    q0 := ((xrhh_hi*yrhh_hi - p0) + xrhh_hi*yrhh_lo + xrhh_lo*yrhh_hi)
          + xrhh_lo*yrhh_lo;
   -- Double_Double_Basics.two_prod(xrhh,yrlh,p1,q1);
    p1 := xrhh*yrlh;
   -- Double_Double_Basics.split(xrhh,xrhh_hi,xrhh_lo); --> done
   -- Double_Double_Basics.split(yrlh,yrlh_hi,yrlh_lo);
    if ( yrlh > QD_SPLIT_THRESH or yrlh < -QD_SPLIT_THRESH ) then
      bb := yrlh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yrlh_hi := s - (s - bb);
      yrlh_lo := yrlh - yrlh_hi;
      yrlh_hi := yrlh_hi*268435456.0;  -- 2^28
      yrlh_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yrlh;
      yrlh_hi := s - (s - yrlh);
      yrlh_lo := yrlh - yrlh_hi;
    end if;
    q1 := ((xrhh_hi*yrlh_hi - p1) + xrhh_hi*yrlh_lo + xrhh_lo*yrlh_hi)
          + xrhh_lo*yrlh_lo;
   -- Double_Double_Basics.two_prod(xrlh,yrhh,p2,q2);
    p2 := xrlh*yrhh;
   -- Double_Double_Basics.split(xrlh,xrlh_hi,xrlh_lo);
    if ( xrlh > QD_SPLIT_THRESH or xrlh < -QD_SPLIT_THRESH ) then
      bb := xrlh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xrlh_hi := s - (s - bb);
      xrlh_lo := xrlh - xrlh_hi;
      xrlh_hi := xrlh_hi*268435456.0;  -- 2^28
      xrlh_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * xrlh;
      xrlh_hi := s - (s - xrlh);
      xrlh_lo := xrlh - xrlh_hi;
    end if;
   -- Double_Double_Basics.split(yrhh,yrhh_hi,yrhh_lo); --> done
    q2 := ((xrlh_hi*yrhh_hi - p2) + xrlh_hi*yrhh_lo + xrlh_lo*yrhh_hi)
          + xrlh_lo*yrhh_lo;
   -- Double_Double_Basics.two_prod(xrhh,yrhl,p3,q3);
    p3 := xrhh*yrhl;
   -- Double_Double_Basics.split(xrhh,xrhh_hi,xrhh_lo); --> done
   -- Double_Double_Basics.split(yrhl,yrhl_hi,yrhl_lo);
    if ( yrhl > QD_SPLIT_THRESH or yrhl < -QD_SPLIT_THRESH ) then
      bb := yrhl*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yrhl_hi := s - (s - bb);
      yrhl_lo := yrhl - yrhl_hi;
      yrhl_hi := yrhl_hi*268435456.0;  -- 2^28
      yrhl_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yrhl;
      yrhl_hi := s - (s - yrhl);
      yrhl_lo := yrhl - yrhl_hi;
    end if;
    q3 := ((xrhh_hi*yrhl_hi - p3) + xrhh_hi*yrhl_lo + xrhh_lo*yrhl_hi)
          + xrhh_lo*yrhl_lo;
   -- Double_Double_Basics.two_prod(xrlh,yrlh,p4,q4);
    p4 := xrlh*yrlh;
   -- Double_Double_Basics.split(xrlh,xrlh_hi,xrlh_lo); --> done
   -- Double_Double_Basics.split(yrlh,yrlh_hi,yrlh_lo); --> done
    q4 := ((xrlh_hi*yrlh_hi - p4) + xrlh_hi*yrlh_lo + xrlh_lo*yrlh_hi)
          + xrlh_lo*yrlh_lo;
   -- Double_Double_Basics.two_prod(xrhl,yrhh,p5,q5);
    p5 := xrhl*yrhh;
   -- Double_Double_Basics.split(xrhl,xrhl_hi,xrhl_lo);
    if ( xrhl > QD_SPLIT_THRESH or xrhl < -QD_SPLIT_THRESH ) then
      bb := xrhl*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xrhl_hi := s - (s - bb);
      xrhl_lo := xrhl - xrhl_hi;
      xrhl_hi := xrhl_hi*268435456.0;  -- 2^28
      xrhl_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * xrhl;
      xrhl_hi := s - (s - xrhl);
      xrhl_lo := xrhl - xrhl_hi;
    end if;
   -- Double_Double_Basics.split(yrhh,yrhh_hi,yrhh_lo); --> done
    q5 := ((xrhl_hi*yrhh_hi - p5) + xrhl_hi*yrhh_lo + xrhl_lo*yrhh_hi)
          + xrhl_lo*yrhh_lo;
   -- Quad_Double_Renormalizations.three_sum(p1,p2,q0);
   -- Double_Double_Basics.two_sum(p1,p2,s0,s1);
    s0 := p1 + p2; bb := s0 - p1; s1 := (p1 - (s0 - bb)) + (p2 - bb);
   -- Double_Double_Basics.two_sum(q0,s0,p1,s2);
    p1 := q0 + s0; bb := p1 - q0; s2 := (q0 - (p1 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p2,q0);
    p2 := s1 + s2; bb := p2 - s1; q0 := (s1 - (p2 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p2,q1,q2);
   -- Double_Double_Basics.two_sum(p2,q1,s0,s1);
    s0 := p2 + q1; bb := s0 - p2; s1 := (p2 - (s0 - bb)) + (q1 - bb);
   -- Double_Double_Basics.two_sum(q2,s0,p2,s2);
    p2 := q2 + s0; bb := p2 - q2; s2 := (q2 - (p2 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,q1,q2);
    q1 := s1 + s2; bb := q1 - s1; q2 := (s1 - (q1 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p3,p4,p5);
   -- Double_Double_Basics.two_sum(p3,p4,s0,s1);
    s0 := p3 + p4; bb := s0 - p3; s1 := (p3 - (s0 - bb)) + (p4 - bb);
   -- Double_Double_Basics.two_sum(p5,s0,p3,s2);
    p3 := p5 + s0; bb := p3 - p5; s2 := (p5 - (p3 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p4,p5);
    p4 := s1 + s2; bb := p4 - s1; p5 := (s1 - (p4 - bb)) + (s2 - bb);
   -- Double_Double_Basics.two_sum(p2,p3,s0,t0);
    s0 := p2 + p3; bb := s0 - p2; t0 := (p2 - (s0 - bb)) + (p3 - bb);
   -- Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s1 := q1 + p4; bb := s1 - q1; t1 := (q1 - (s1 - bb)) + (p4 - bb);
    s2 := q2 + p5;
   -- Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s := s1 + t0; bb := s - s1; t0 := (s1 - (s - bb)) + (t0 - bb); s1 := s;
    s2 := s2 + (t0 + t1);
   -- Double_Double_Basics.two_prod(xrhh,yrll,p6,q6);
    p6 := xrhh*yrll;
   -- Double_Double_Basics.split(xrhh,xrhh_hi,xrhh_lo); --> done
   -- Double_Double_Basics.split(yrll,yrll_hi,yrll_lo);
    if ( yrll > QD_SPLIT_THRESH or yrll < -QD_SPLIT_THRESH ) then
      bb := yrll*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yrll_hi := s - (s - bb);
      yrll_lo := yrll - yrll_hi;
      yrll_hi := yrll_hi*268435456.0;  -- 2^28
      yrll_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yrll;
      yrll_hi := s - (s - yrll);
      yrll_lo := yrll - yrll_hi;
    end if;
    q6 := ((xrhh_hi*yrll_hi - p6) + xrhh_hi*yrll_lo + xrhh_lo*yrll_hi)
          + xrhh_lo*yrll_lo;
   -- Double_Double_Basics.two_prod(xrlh,yrhl,p7,q7);
    p7 := xrlh*yrhl;
   -- Double_Double_Basics.split(xrlh,xrlh_hi,xrlh_lo); --> done
   -- Double_Double_Basics.split(yrhl,yrhl_hi,yrhl_lo); --> done
    q7 := ((xrlh_hi*yrhl_hi - p7) + xrlh_hi*yrhl_lo + xrlh_lo*yrhl_hi)
          + xrlh_lo*yrhl_lo;
   -- Double_Double_Basics.two_prod(xrhl,yrlh,p8,q8);
    p8 := xrhl*yrlh;
   -- Double_Double_Basics.split(xrhl,xrhl_hi,xrhl_lo); --> done
   -- Double_Double_Basics.split(yrlh,yrlh_hi,yrlh_lo); --> done
    q8 := ((xrhl_hi*yrlh_hi - p8) + xrhl_hi*yrlh_lo + xrhl_lo*yrlh_hi)
          + xrhl_lo*yrlh_lo;
   -- Double_Double_Basics.two_prod(xrll,yrhh,p9,q9);
    p9 := xrll*yrhh;
   -- Double_Double_Basics.split(xrll,xrll_hi,xrll_lo);
    if ( xrll > QD_SPLIT_THRESH or xrll < -QD_SPLIT_THRESH ) then
      bb := xrll*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xrll_hi := s - (s - bb);
      xrll_lo := xrll - xrll_hi;
      xrll_hi := xrll_hi*268435456.0;  -- 2^28
      xrll_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * xrll;
      xrll_hi := s - (s - xrll);
      xrll_lo := xrll - xrll_hi;
    end if;
   -- Double_Double_Basics.split(yrhh,yrhh_hi,yrhh_lo); --> done
    q9 := ((xrll_hi*yrhh_hi - p9) + xrll_hi*yrhh_lo + xrll_lo*yrhh_hi)
          + xrll_lo*yrhh_lo;
   -- Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    s := q0 + q3; bb := s - q0; q3 := (q0 - (s - bb)) + (q3 - bb); q0 := s;
   -- Double_Double_Basics.two_sum(q4,q5,q4,q5);
    s := q4 + q5; bb := s - q4; q5 := (q4 - (s - bb)) + (q5 - bb); q4 := s;
   -- Double_Double_Basics.two_sum(p6,p7,p6,p7);
    s := p6 + p7; bb := s - p6; p7 := (p6 - (s - bb)) + (p7 - bb); p6 := s;
   -- Double_Double_Basics.two_sum(p8,p9,p8,p9);
    s := p8 + p9; bb := s - p8; p9 := (p8 - (s - bb)) + (p9 - bb); p8 := s;
   -- Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t0 := q0 + q4; bb := t0 - q0; t1 := (q0 - (t0 - bb)) + (q4 - bb);
    t1 := t1 + (q3 + q5);
   -- Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r0 := p6 + p8; bb := r0 - p6; r1 := (p6 - (r0 - bb)) + (p8 - bb);
    r1 := r1 + (p7 + p9); 
   -- Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q3 := t0 + r0; bb := q3 - t0; q4 := (t0 - (q3 - bb)) + (r0 - bb);
    q4 := q4 + (t1 + r1); 
   -- Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t0 := q3 + s1; bb := t0 - q3; t1 := (q3 - (t0 - bb)) + (s1 - bb);
    t1 := t1 + q4;
    t1 := t1 + xrlh * yrll + xrhl * yrhl + xrll * yrlh
        + q6 + q7 + q8 + q9 + s2;
   -- Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
   -- s0 = bb, c0 = p0, c1 = p1, c2 = s0, c3 = t0, c4 = t1
   -- Double_Double_Basics.quick_two_sum(t0,t1,bb,t1);
    bb := t0 + t1; t1 := t1 - (bb - t0);
   -- Double_Double_Basics.quick_two_sum(s0,bb,bb,t0);
    s := s0 + bb; t0 := bb - (s - s0); bb := s;
   -- Double_Double_Basics.quick_two_sum(p1,bb,bb,s0);
    s := p1 + bb; s0 := bb - (s - p1); bb := s;
   -- Double_Double_Basics.quick_two_sum(p0,bb,p0,p1);
    s := p0 + bb; p1 := bb - (s - p0); p0 := s; -- bb := p0; s1 := p1;
   -- Double_Double_Basics.quick_two_sum(p0,p1,bb,s1);
    bb := p0 + p1; s1 := p1 - (bb - p0);
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,s0,s1,s2);
      s := s1 + s0; s2 := s0 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,t0,s2,s3);
        s := s2 + t0; s3 := t0 - (s - s2); s2 := s;
        if s3 /= 0.0
         then s3 := s3 + t1;
         else s2 := s2 + t1;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(bb,s0,bb,s1);
      s := bb + s0; s1 := s0 - (s - bb); bb := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(bb,t0,bb,s1);
        s := bb + t0; s1 := t0 - (s - bb); bb := s;
        if s1 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        else
         -- Double_Double_Basics.quick_two_sum(bb,t1,bb,s1);
          s := bb + t1; s1 := t1 - (s - bb); bb := s;
        end if;
      end if;
    end if;
    p0 := bb; p1 := s1; s0 := s2; t0 := s3;
   -- zhihi := p0; zlohi := p1; zhilo := s0; zlolo := t0;
   -- (2) update relevant part of zr variables
    x(0) := zrhh; x(1) := zrlh; x(2) := zrhl; x(3) := zrll;
    y(0) := p0;   y(1) := p1;   y(2) := s0;   y(3) := t0;
    Add(x,y,z);
    zrhh := z(0); zrlh := z(1); zrhl := z(2); zrll := z(3);
   -- (3) compute xim*yim
   -- Double_Double_Basics.two_prod(xihh,yihh,p0,q0);
    p0 := xihh*yihh;
   -- Double_Double_Basics.split(xihh,xihh_hi,xihh_lo);
    if ( xihh > QD_SPLIT_THRESH or xihh < -QD_SPLIT_THRESH ) then
      bb := xihh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xihh_hi := s - (s - bb);
      xihh_lo := xihh - xihh_hi;
      xihh_hi := xihh_hi*268435456.0;  -- 2^28
      xihh_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * xihh;
      xihh_hi := s - (s - xihh);
      xihh_lo := xihh - xihh_hi;
    end if;
   -- Double_Double_Basics.split(yihh,yihh_hi,yihh_lo);
    if ( yihh > QD_SPLIT_THRESH or yihh < -QD_SPLIT_THRESH ) then
      bb := yihh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yihh_hi := s - (s - bb);
      yihh_lo := yihh - yihh_hi;
      yihh_hi := yihh_hi*268435456.0;  -- 2^28
      yihh_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yihh;
      yihh_hi := s - (s - yihh);
      yihh_lo := yihh - yihh_hi;
    end if;
    q0 := ((xihh_hi*yihh_hi - p0) + xihh_hi*yihh_lo + xihh_lo*yihh_hi)
          + xihh_lo*yihh_lo;
   -- Double_Double_Basics.two_prod(xihh,yilh,p1,q1);
    p1 := xihh*yilh;
   -- Double_Double_Basics.split(xihh,xihh_hi,xihh_lo); --> done
   -- Double_Double_Basics.split(yilh,yilh_hi,yilh_lo);
    if ( yilh > QD_SPLIT_THRESH or yilh < -QD_SPLIT_THRESH ) then
      bb := yilh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yilh_hi := s - (s - bb);
      yilh_lo := yilh - yilh_hi;
      yilh_hi := yilh_hi*268435456.0;  -- 2^28
      yilh_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yilh;
      yilh_hi := s - (s - yilh);
      yilh_lo := yilh - yilh_hi;
    end if;
    q1 := ((xihh_hi*yilh_hi - p1) + xihh_hi*yilh_lo + xihh_lo*yilh_hi)
          + xihh_lo*yilh_lo;
   -- Double_Double_Basics.two_prod(xilh,yihh,p2,q2);
    p2 := xilh*yihh;
   -- Double_Double_Basics.split(xilh,xilh_hi,xilh_lo);
    if ( xilh > QD_SPLIT_THRESH or xilh < -QD_SPLIT_THRESH ) then
      bb := xilh*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xilh_hi := s - (s - bb);
      xilh_lo := xilh - xilh_hi;
      xilh_hi := xilh_hi*268435456.0;  -- 2^28
      xilh_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * xilh;
      xilh_hi := s - (s - xilh);
      xilh_lo := xilh - xilh_hi;
    end if;
   -- Double_Double_Basics.split(yihh,yihh_hi,yihh_lo); --> done
    q2 := ((xilh_hi*yihh_hi - p2) + xilh_hi*yihh_lo + xilh_lo*yihh_hi)
          + xilh_lo*yihh_lo;
   -- Double_Double_Basics.two_prod(xihh,yihl,p3,q3);
    p3 := xihh*yihl;
   -- Double_Double_Basics.split(xihh,xihh_hi,xihh_lo); --> done
   -- Double_Double_Basics.split(yihl,yihl_hi,yihl_lo);
    if ( yihl > QD_SPLIT_THRESH or yihl < -QD_SPLIT_THRESH ) then
      bb := yihl*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yihl_hi := s - (s - bb);
      yihl_lo := yilh - yihl_hi;
      yihl_hi := yihl_hi*268435456.0;  -- 2^28
      yihl_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yihl;
      yihl_hi := s - (s - yihl);
      yihl_lo := yihl - yihl_hi;
    end if;
    q3 := ((xihh_hi*yihl_hi - p3) + xihh_hi*yihl_lo + xihh_lo*yihl_hi)
          + xihh_lo*yihl_lo;
   -- Double_Double_Basics.two_prod(xilh,yilh,p4,q4);
    p4 := xilh*yilh;
   -- Double_Double_Basics.split(xilh,xilh_hi,xilh_lo); --> done
   -- Double_Double_Basics.split(yilh,yilh_hi,yilh_lo); --> done
    q4 := ((xilh_hi*yilh_hi - p4) + xilh_hi*yilh_lo + xilh_lo*yilh_hi)
          + xilh_lo*yilh_lo;
   -- Double_Double_Basics.two_prod(xihl,yihh,p5,q5);
    p5 := xihl*yihh;
   -- Double_Double_Basics.split(xihl,xihl_hi,xihl_lo);
    if ( xihl > QD_SPLIT_THRESH or xihl < -QD_SPLIT_THRESH ) then
      bb := xihl*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xihl_hi := s - (s - bb);
      xihl_lo := xihl - xihl_hi;
      xihl_hi := xihl_hi*268435456.0;  -- 2^28
      xihl_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * xihl;
      xihl_hi := s - (s - xihl);
      xihl_lo := xihl - xihl_hi;
    end if;
   -- Double_Double_Basics.split(yihh,yihh_hi,yihh_lo); --> done
    q5 := ((xihl_hi*yihh_hi - p5) + xihl_hi*yihh_lo + xihl_lo*yihh_hi)
          + xihl_lo*yihh_lo;
   -- Quad_Double_Renormalizations.three_sum(p1,p2,q0);
   -- Double_Double_Basics.two_sum(p1,p2,s0,s1);
    s0 := p1 + p2; bb := s0 - p1; s1 := (p1 - (s0 - bb)) + (p2 - bb);
   -- Double_Double_Basics.two_sum(q0,s0,p1,s2);
    p1 := q0 + s0; bb := p1 - q0; s2 := (q0 - (p1 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p2,q0);
    p2 := s1 + s2; bb := p2 - s1; q0 := (s1 - (p2 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p2,q1,q2);
   -- Double_Double_Basics.two_sum(p2,q1,s0,s1);
    s0 := p2 + q1; bb := s0 - p2; s1 := (p2 - (s0 - bb)) + (q1 - bb);
   -- Double_Double_Basics.two_sum(q2,s0,p2,s2);
    p2 := q2 + s0; bb := p2 - q2; s2 := (q2 - (p2 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,q1,q2);
    q1 := s1 + s2; bb := q1 - s1; q2 := (s1 - (q1 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p3,p4,p5);
   -- Double_Double_Basics.two_sum(p3,p4,s0,s1);
    s0 := p3 + p4; bb := s0 - p3; s1 := (p3 - (s0 - bb)) + (p4 - bb);
   -- Double_Double_Basics.two_sum(p5,s0,p3,s2);
    p3 := p5 + s0; bb := p3 - p5; s2 := (p5 - (p3 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p4,p5);
    p4 := s1 + s2; bb := p4 - s1; p5 := (s1 - (p4 - bb)) + (s2 - bb);
   -- Double_Double_Basics.two_sum(p2,p3,s0,t0);
    s0 := p2 + p3; bb := s0 - p2; t0 := (p2 - (s0 - bb)) + (p3 - bb);
   -- Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s1 := q1 + p4; bb := s1 - q1; t1 := (q1 - (s1 - bb)) + (p4 - bb);
    s2 := q2 + p5;
   -- Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s := s1 + t0; bb := s - s1; t0 := (s1 - (s - bb)) + (t0 - bb); s1 := s;
    s2 := s2 + (t0 + t1);
   -- Double_Double_Basics.two_prod(xihh,yill,p6,q6);
    p6 := xihh*yill;
   -- Double_Double_Basics.split(xihh,xihh_hi,xihh_lo); --> done
   -- Double_Double_Basics.split(yill,yill_hi,yill_lo);
    if ( yill > QD_SPLIT_THRESH or yill < -QD_SPLIT_THRESH ) then
      bb := yill*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      yill_hi := s - (s - bb);
      yill_lo := yill - yill_hi;
      yill_hi := yill_hi*268435456.0;  -- 2^28
      yill_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * yill;
      yill_hi := s - (s - yill);
      yill_lo := yill - yill_hi;
    end if;
    q6 := ((xihh_hi*yill_hi - p6) + xihh_hi*yill_lo + xihh_lo*yill_hi)
          + xihh_lo*yill_lo;
   -- Double_Double_Basics.two_prod(xilh,yihl,p7,q7);
    p7 := xilh*yihl;
   -- Double_Double_Basics.split(xilh,xilh_hi,xilh_lo); --> done
   -- Double_Double_Basics.split(yihl,yihl_hi,yihl_lo); --> done
    q7 := ((xilh_hi*yihl_hi - p7) + xilh_hi*yihl_lo + xilh_lo*yihl_hi)
          + xilh_lo*yihl_lo;
   -- Double_Double_Basics.two_prod(xihl,yilh,p8,q8);
    p8 := xihl*yilh;
   -- Double_Double_Basics.split(xihl,xihl_hi,xihl_lo); --> done
   -- Double_Double_Basics.split(yilh,yilh_hi,yilh_lo); --> done
    q8 := ((xihl_hi*yilh_hi - p8) + xihl_hi*yilh_lo + xihl_lo*yilh_hi)
          + xihl_lo*yilh_lo;
   -- Double_Double_Basics.two_prod(xill,yihh,p9,q9);
    p9 := xill*yihh;
   -- Double_Double_Basics.split(xill,xill_hi,xill_lo);
    if ( xill > QD_SPLIT_THRESH or xill < -QD_SPLIT_THRESH ) then
      bb := xill*3.7252902984619140625E-09;  -- 2^-28
      s := QD_SPLITTER * bb;
      xill_hi := s - (s - bb);
      xill_lo := xill - xill_hi;
      xill_hi := xill_hi*268435456.0;  -- 2^28
      xill_lo := 268435456.0;     -- 2^28
    else
      s := QD_SPLITTER * xill;
      xill_hi := s - (s - xill);
      xill_lo := xill - xill_hi;
    end if;
   -- Double_Double_Basics.split(yihh,yihh_hi,yihh_lo); --> done
    q9 := ((xill_hi*yihh_hi - p9) + xill_hi*yihh_lo + xill_lo*yihh_hi)
          + xill_lo*yihh_lo;
   -- Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    s := q0 + q3; bb := s - q0; q3 := (q0 - (s - bb)) + (q3 - bb); q0 := s;
   -- Double_Double_Basics.two_sum(q4,q5,q4,q5);
    s := q4 + q5; bb := s - q4; q5 := (q4 - (s - bb)) + (q5 - bb); q4 := s;
   -- Double_Double_Basics.two_sum(p6,p7,p6,p7);
    s := p6 + p7; bb := s - p6; p7 := (p6 - (s - bb)) + (p7 - bb); p6 := s;
   -- Double_Double_Basics.two_sum(p8,p9,p8,p9);
    s := p8 + p9; bb := s - p8; p9 := (p8 - (s - bb)) + (p9 - bb); p8 := s;
   -- Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t0 := q0 + q4; bb := t0 - q0; t1 := (q0 - (t0 - bb)) + (q4 - bb);
    t1 := t1 + (q3 + q5);
   -- Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r0 := p6 + p8; bb := r0 - p6; r1 := (p6 - (r0 - bb)) + (p8 - bb);
    r1 := r1 + (p7 + p9); 
   -- Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q3 := t0 + r0; bb := q3 - t0; q4 := (t0 - (q3 - bb)) + (r0 - bb);
    q4 := q4 + (t1 + r1); 
   -- Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t0 := q3 + s1; bb := t0 - q3; t1 := (q3 - (t0 - bb)) + (s1 - bb);
    t1 := t1 + q4;
    t1 := t1 + xilh * yill + xihl * yihl + xill * yilh
        + q6 + q7 + q8 + q9 + s2;
   -- Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
   -- s0 = bb, c0 = p0, c1 = p1, c2 = s0, c3 = t0, c4 = t1
   -- Double_Double_Basics.quick_two_sum(t0,t1,bb,t1);
    bb := t0 + t1; t1 := t1 - (bb - t0);
   -- Double_Double_Basics.quick_two_sum(s0,bb,bb,t0);
    s := s0 + bb; t0 := bb - (s - s0); bb := s;
   -- Double_Double_Basics.quick_two_sum(p1,bb,bb,s0);
    s := p1 + bb; s0 := bb - (s - p1); bb := s;
   -- Double_Double_Basics.quick_two_sum(p0,bb,p0,p1);
    s := p0 + bb; p1 := bb - (s - p0); p0 := s; -- bb := p0; s1 := p1;
   -- Double_Double_Basics.quick_two_sum(p0,p1,bb,s1);
    bb := p0 + p1; s1 := p1 - (bb - p0);
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,s0,s1,s2);
      s := s1 + s0; s2 := s0 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,t0,s2,s3);
        s := s2 + t0; s3 := t0 - (s - s2); s2 := s;
        if s3 /= 0.0
         then s3 := s3 + t1;
         else s2 := s2 + t1;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(bb,s0,bb,s1);
      s := bb + s0; s1 := s0 - (s - bb); bb := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(bb,t0,bb,s1);
        s := bb + t0; s1 := t0 - (s - bb); bb := s;
        if s1 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        else
         -- Double_Double_Basics.quick_two_sum(bb,t1,bb,s1);
          s := bb + t1; s1 := t1 - (s - bb); bb := s;
        end if;
      end if;
    end if;
    p0 := bb; p1 := s1; s0 := s2; t0 := s3;
   -- zhihi := p0; zlohi := p1; zhilo := s0; zlolo := t0;
   -- (4) update relevant part of zr variables
    x(0) := zrhh; x(1) := zrlh; x(2) := zrhl; x(3) := zrll;
    y(0) := -p0;  y(1) := -p1;  y(2) := -s0;  y(3) := -t0;
    Add(x,y,z);
    zrhh := z(0); zrlh := z(1); zrhl := z(2); zrll := z(3);
   -- zim = xre*yim + xim * yre
   -- (5) compute xre*yim
   -- Double_Double_Basics.two_prod(xrhh,yihh,p0,q0);
    p0 := xrhh*yihh;
    q0 := ((xrhh_hi*yihh_hi - p0) + xrhh_hi*yihh_lo + xrhh_lo*yihh_hi)
          + xrhh_lo*yihh_lo;
   -- Double_Double_Basics.two_prod(xrhh,yilh,p1,q1);
    p1 := xrhh*yilh;
    q1 := ((xrhh_hi*yilh_hi - p1) + xrhh_hi*yilh_lo + xrhh_lo*yilh_hi)
          + xrhh_lo*yilh_lo;
   -- Double_Double_Basics.two_prod(xrlh,yihh,p2,q2);
    p2 := xrlh*yihh;
    q2 := ((xrlh_hi*yihh_hi - p2) + xrlh_hi*yihh_lo + xrlh_lo*yihh_hi)
          + xrlh_lo*yihh_lo;
   -- Double_Double_Basics.two_prod(xrhh,yihl,p3,q3);
    p3 := xrhh*yihl;
    q3 := ((xrhh_hi*yihl_hi - p3) + xrhh_hi*yihl_lo + xrhh_lo*yihl_hi)
          + xrhh_lo*yihl_lo;
   -- Double_Double_Basics.two_prod(xrlh,yilh,p4,q4);
    p4 := xrlh*yilh;
    q4 := ((xrlh_hi*yilh_hi - p4) + xrlh_hi*yilh_lo + xrlh_lo*yilh_hi)
          + xrlh_lo*yilh_lo;
   -- Double_Double_Basics.two_prod(xrhl,yihh,p5,q5);
    p5 := xrhl*yihh;
    q5 := ((xrhl_hi*yihh_hi - p5) + xrhl_hi*yihh_lo + xrhl_lo*yihh_hi)
          + xrhl_lo*yihh_lo;
   -- Quad_Double_Renormalizations.three_sum(p1,p2,q0);
   -- Double_Double_Basics.two_sum(p1,p2,s0,s1);
    s0 := p1 + p2; bb := s0 - p1; s1 := (p1 - (s0 - bb)) + (p2 - bb);
   -- Double_Double_Basics.two_sum(q0,s0,p1,s2);
    p1 := q0 + s0; bb := p1 - q0; s2 := (q0 - (p1 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p2,q0);
    p2 := s1 + s2; bb := p2 - s1; q0 := (s1 - (p2 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p2,q1,q2);
   -- Double_Double_Basics.two_sum(p2,q1,s0,s1);
    s0 := p2 + q1; bb := s0 - p2; s1 := (p2 - (s0 - bb)) + (q1 - bb);
   -- Double_Double_Basics.two_sum(q2,s0,p2,s2);
    p2 := q2 + s0; bb := p2 - q2; s2 := (q2 - (p2 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,q1,q2);
    q1 := s1 + s2; bb := q1 - s1; q2 := (s1 - (q1 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p3,p4,p5);
   -- Double_Double_Basics.two_sum(p3,p4,s0,s1);
    s0 := p3 + p4; bb := s0 - p3; s1 := (p3 - (s0 - bb)) + (p4 - bb);
   -- Double_Double_Basics.two_sum(p5,s0,p3,s2);
    p3 := p5 + s0; bb := p3 - p5; s2 := (p5 - (p3 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p4,p5);
    p4 := s1 + s2; bb := p4 - s1; p5 := (s1 - (p4 - bb)) + (s2 - bb);
   -- Double_Double_Basics.two_sum(p2,p3,s0,t0);
    s0 := p2 + p3; bb := s0 - p2; t0 := (p2 - (s0 - bb)) + (p3 - bb);
   -- Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s1 := q1 + p4; bb := s1 - q1; t1 := (q1 - (s1 - bb)) + (p4 - bb);
    s2 := q2 + p5;
   -- Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s := s1 + t0; bb := s - s1; t0 := (s1 - (s - bb)) + (t0 - bb); s1 := s;
    s2 := s2 + (t0 + t1);
   -- Double_Double_Basics.two_prod(xrhh,yill,p6,q6);
    p6 := xrhh*yill;
    q6 := ((xrhh_hi*yill_hi - p6) + xrhh_hi*yill_lo + xrhh_lo*yill_hi)
          + xrhh_lo*yill_lo;
   -- Double_Double_Basics.two_prod(xrlh,yihl,p7,q7);
    p7 := xrlh*yihl;
    q7 := ((xrlh_hi*yihl_hi - p7) + xrlh_hi*yihl_lo + xrlh_lo*yihl_hi)
          + xrlh_lo*yihl_lo;
   -- Double_Double_Basics.two_prod(xrhl,yilh,p8,q8);
    p8 := xrhl*yilh;
    q8 := ((xrhl_hi*yilh_hi - p8) + xrhl_hi*yilh_lo + xrhl_lo*yilh_hi)
          + xrhl_lo*yilh_lo;
   -- Double_Double_Basics.two_prod(xrll,yihh,p9,q9);
    p9 := xrll*yihh;
    q9 := ((xrll_hi*yihh_hi - p9) + xrll_hi*yihh_lo + xrll_lo*yihh_hi)
          + xrll_lo*yihh_lo;
   -- Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    s := q0 + q3; bb := s - q0; q3 := (q0 - (s - bb)) + (q3 - bb); q0 := s;
   -- Double_Double_Basics.two_sum(q4,q5,q4,q5);
    s := q4 + q5; bb := s - q4; q5 := (q4 - (s - bb)) + (q5 - bb); q4 := s;
   -- Double_Double_Basics.two_sum(p6,p7,p6,p7);
    s := p6 + p7; bb := s - p6; p7 := (p6 - (s - bb)) + (p7 - bb); p6 := s;
   -- Double_Double_Basics.two_sum(p8,p9,p8,p9);
    s := p8 + p9; bb := s - p8; p9 := (p8 - (s - bb)) + (p9 - bb); p8 := s;
   -- Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t0 := q0 + q4; bb := t0 - q0; t1 := (q0 - (t0 - bb)) + (q4 - bb);
    t1 := t1 + (q3 + q5);
   -- Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r0 := p6 + p8; bb := r0 - p6; r1 := (p6 - (r0 - bb)) + (p8 - bb);
    r1 := r1 + (p7 + p9); 
   -- Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q3 := t0 + r0; bb := q3 - t0; q4 := (t0 - (q3 - bb)) + (r0 - bb);
    q4 := q4 + (t1 + r1); 
   -- Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t0 := q3 + s1; bb := t0 - q3; t1 := (q3 - (t0 - bb)) + (s1 - bb);
    t1 := t1 + q4;
    t1 := t1 + xrlh * yill + xrhl * yihl + xrll * yilh
        + q6 + q7 + q8 + q9 + s2;
   -- Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
   -- s0 = bb, c0 = p0, c1 = p1, c2 = s0, c3 = t0, c4 = t1
   -- Double_Double_Basics.quick_two_sum(t0,t1,bb,t1);
    bb := t0 + t1; t1 := t1 - (bb - t0);
   -- Double_Double_Basics.quick_two_sum(s0,bb,bb,t0);
    s := s0 + bb; t0 := bb - (s - s0); bb := s;
   -- Double_Double_Basics.quick_two_sum(p1,bb,bb,s0);
    s := p1 + bb; s0 := bb - (s - p1); bb := s;
   -- Double_Double_Basics.quick_two_sum(p0,bb,p0,p1);
    s := p0 + bb; p1 := bb - (s - p0); p0 := s; -- bb := p0; s1 := p1;
   -- Double_Double_Basics.quick_two_sum(p0,p1,bb,s1);
    bb := p0 + p1; s1 := p1 - (bb - p0);
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,s0,s1,s2);
      s := s1 + s0; s2 := s0 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,t0,s2,s3);
        s := s2 + t0; s3 := t0 - (s - s2); s2 := s;
        if s3 /= 0.0
         then s3 := s3 + t1;
         else s2 := s2 + t1;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(bb,s0,bb,s1);
      s := bb + s0; s1 := s0 - (s - bb); bb := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(bb,t0,bb,s1);
        s := bb + t0; s1 := t0 - (s - bb); bb := s;
        if s1 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        else
         -- Double_Double_Basics.quick_two_sum(bb,t1,bb,s1);
          s := bb + t1; s1 := t1 - (s - bb); bb := s;
        end if;
      end if;
    end if;
    p0 := bb; p1 := s1; s0 := s2; t0 := s3;
   -- zhihi := p0; zlohi := p1; zhilo := s0; zlolo := t0;
   -- (6) update relevant part of zi variables
    x(0) := zihh; x(1) := zilh; x(2) := zihl; x(3) := zill;
    y(0) := p0;   y(1) := p1;   y(2) := s0;   y(3) := t0;
    Add(x,y,z);
    zihh := z(0); zilh := z(1); zihl := z(2); zill := z(3);
   -- (7) compute xim*yre
   -- Double_Double_Basics.two_prod(xihh,yrhh,p0,q0);
    p0 := xihh*yrhh;
    q0 := ((xihh_hi*yrhh_hi - p0) + xihh_hi*yrhh_lo + xihh_lo*yrhh_hi)
          + xihh_lo*yrhh_lo;
   -- Double_Double_Basics.two_prod(xihh,yrlh,p1,q1);
    p1 := xihh*yrlh;
    q1 := ((xihh_hi*yrlh_hi - p1) + xihh_hi*yrlh_lo + xihh_lo*yrlh_hi)
          + xihh_lo*yrlh_lo;
   -- Double_Double_Basics.two_prod(xilh,yrhh,p2,q2);
    p2 := xilh*yrhh;
    q2 := ((xilh_hi*yrhh_hi - p2) + xilh_hi*yrhh_lo + xilh_lo*yrhh_hi)
          + xilh_lo*yrhh_lo;
   -- Double_Double_Basics.two_prod(xihh,yrhl,p3,q3);
    p3 := xihh*yrhl;
    q3 := ((xihh_hi*yrhl_hi - p3) + xihh_hi*yrhl_lo + xihh_lo*yrhl_hi)
          + xihh_lo*yrhl_lo;
   -- Double_Double_Basics.two_prod(xilh,yrlh,p4,q4);
    p4 := xilh*yrlh;
    q4 := ((xilh_hi*yrlh_hi - p4) + xilh_hi*yrlh_lo + xilh_lo*yrlh_hi)
          + xilh_lo*yrlh_lo;
   -- Double_Double_Basics.two_prod(xihl,yrhh,p5,q5);
    p5 := xihl*yrhh;
    q5 := ((xihl_hi*yrhh_hi - p5) + xihl_hi*yrhh_lo + xihl_lo*yrhh_hi)
          + xihl_lo*yrhh_lo;
   -- Quad_Double_Renormalizations.three_sum(p1,p2,q0);
   -- Double_Double_Basics.two_sum(p1,p2,s0,s1);
    s0 := p1 + p2; bb := s0 - p1; s1 := (p1 - (s0 - bb)) + (p2 - bb);
   -- Double_Double_Basics.two_sum(q0,s0,p1,s2);
    p1 := q0 + s0; bb := p1 - q0; s2 := (q0 - (p1 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p2,q0);
    p2 := s1 + s2; bb := p2 - s1; q0 := (s1 - (p2 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p2,q1,q2);
   -- Double_Double_Basics.two_sum(p2,q1,s0,s1);
    s0 := p2 + q1; bb := s0 - p2; s1 := (p2 - (s0 - bb)) + (q1 - bb);
   -- Double_Double_Basics.two_sum(q2,s0,p2,s2);
    p2 := q2 + s0; bb := p2 - q2; s2 := (q2 - (p2 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,q1,q2);
    q1 := s1 + s2; bb := q1 - s1; q2 := (s1 - (q1 - bb)) + (s2 - bb);
   -- Quad_Double_Renormalizations.three_sum(p3,p4,p5);
   -- Double_Double_Basics.two_sum(p3,p4,s0,s1);
    s0 := p3 + p4; bb := s0 - p3; s1 := (p3 - (s0 - bb)) + (p4 - bb);
   -- Double_Double_Basics.two_sum(p5,s0,p3,s2);
    p3 := p5 + s0; bb := p3 - p5; s2 := (p5 - (p3 - bb)) + (s0 - bb);
   -- Double_Double_Basics.two_sum(s1,s2,p4,p5);
    p4 := s1 + s2; bb := p4 - s1; p5 := (s1 - (p4 - bb)) + (s2 - bb);
   -- Double_Double_Basics.two_sum(p2,p3,s0,t0);
    s0 := p2 + p3; bb := s0 - p2; t0 := (p2 - (s0 - bb)) + (p3 - bb);
   -- Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s1 := q1 + p4; bb := s1 - q1; t1 := (q1 - (s1 - bb)) + (p4 - bb);
    s2 := q2 + p5;
   -- Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s := s1 + t0; bb := s - s1; t0 := (s1 - (s - bb)) + (t0 - bb); s1 := s;
    s2 := s2 + (t0 + t1);
   -- Double_Double_Basics.two_prod(xihh,yrll,p6,q6);
    p6 := xihh*yrll;
    q6 := ((xihh_hi*yrll_hi - p6) + xihh_hi*yrll_lo + xihh_lo*yrll_hi)
          + xihh_lo*yrll_lo;
   -- Double_Double_Basics.two_prod(xilh,yrhl,p7,q7);
    p7 := xilh*yrhl;
    q7 := ((xilh_hi*yrhl_hi - p7) + xilh_hi*yrhl_lo + xilh_lo*yrhl_hi)
          + xilh_lo*yrhl_lo;
   -- Double_Double_Basics.two_prod(xihl,yrlh,p8,q8);
    p8 := xihl*yrlh;
    q8 := ((xihl_hi*yrlh_hi - p8) + xihl_hi*yrlh_lo + xihl_lo*yrlh_hi)
          + xihl_lo*yrlh_lo;
   -- Double_Double_Basics.two_prod(xill,yrhh,p9,q9);
    p9 := xill*yrhh;
    q9 := ((xill_hi*yrhh_hi - p9) + xill_hi*yrhh_lo + xill_lo*yrhh_hi)
          + xill_lo*yrhh_lo;
   -- Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    s := q0 + q3; bb := s - q0; q3 := (q0 - (s - bb)) + (q3 - bb); q0 := s;
   -- Double_Double_Basics.two_sum(q4,q5,q4,q5);
    s := q4 + q5; bb := s - q4; q5 := (q4 - (s - bb)) + (q5 - bb); q4 := s;
   -- Double_Double_Basics.two_sum(p6,p7,p6,p7);
    s := p6 + p7; bb := s - p6; p7 := (p6 - (s - bb)) + (p7 - bb); p6 := s;
   -- Double_Double_Basics.two_sum(p8,p9,p8,p9);
    s := p8 + p9; bb := s - p8; p9 := (p8 - (s - bb)) + (p9 - bb); p8 := s;
   -- Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t0 := q0 + q4; bb := t0 - q0; t1 := (q0 - (t0 - bb)) + (q4 - bb);
    t1 := t1 + (q3 + q5);
   -- Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r0 := p6 + p8; bb := r0 - p6; r1 := (p6 - (r0 - bb)) + (p8 - bb);
    r1 := r1 + (p7 + p9); 
   -- Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q3 := t0 + r0; bb := q3 - t0; q4 := (t0 - (q3 - bb)) + (r0 - bb);
    q4 := q4 + (t1 + r1); 
   -- Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t0 := q3 + s1; bb := t0 - q3; t1 := (q3 - (t0 - bb)) + (s1 - bb);
    t1 := t1 + q4;
    t1 := t1 + xilh * yrll + xihl * yrhl + xill * yrlh
        + q6 + q7 + q8 + q9 + s2;
   -- Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
   -- s0 = bb, c0 = p0, c1 = p1, c2 = s0, c3 = t0, c4 = t1
   -- Double_Double_Basics.quick_two_sum(t0,t1,bb,t1);
    bb := t0 + t1; t1 := t1 - (bb - t0);
   -- Double_Double_Basics.quick_two_sum(s0,bb,bb,t0);
    s := s0 + bb; t0 := bb - (s - s0); bb := s;
   -- Double_Double_Basics.quick_two_sum(p1,bb,bb,s0);
    s := p1 + bb; s0 := bb - (s - p1); bb := s;
   -- Double_Double_Basics.quick_two_sum(p0,bb,p0,p1);
    s := p0 + bb; p1 := bb - (s - p0); p0 := s; -- bb := p0; s1 := p1;
   -- Double_Double_Basics.quick_two_sum(p0,p1,bb,s1);
    bb := p0 + p1; s1 := p1 - (bb - p0);
    if s1 /= 0.0 then
     -- Double_Double_Basics.quick_two_sum(s1,s0,s1,s2);
      s := s1 + s0; s2 := s0 - (s - s1); s1 := s;
      if s2 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s2,t0,s2,s3);
        s := s2 + t0; s3 := t0 - (s - s2); s2 := s;
        if s3 /= 0.0
         then s3 := s3 + t1;
         else s2 := s2 + t1;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      end if;
    else
     -- Double_Double_Basics.quick_two_sum(bb,s0,bb,s1);
      s := bb + s0; s1 := s0 - (s - bb); bb := s;
      if s1 /= 0.0 then
       -- Double_Double_Basics.quick_two_sum(s1,t0,s1,s2);
        s := s1 + t0; s2 := t0 - (s - s1); s1 := s;
        if s2 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s2,t1,s2,s3);
          s := s2 + t1; s3 := t1 - (s - s2); s2 := s;
        else
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        end if;
      else
       -- Double_Double_Basics.quick_two_sum(bb,t0,bb,s1);
        s := bb + t0; s1 := t0 - (s - bb); bb := s;
        if s1 /= 0.0 then
         -- Double_Double_Basics.quick_two_sum(s1,t1,s1,s2);
          s := s1 + t1; s2 := t1 - (s - s1); s1 := s;
        else
         -- Double_Double_Basics.quick_two_sum(bb,t1,bb,s1);
          s := bb + t1; s1 := t1 - (s - bb); bb := s;
        end if;
      end if;
    end if;
    p0 := bb; p1 := s1; s0 := s2; t0 := s3;
   -- zhihi := p0; zlohi := p1; zhilo := s0; zlolo := t0;
   -- (6) update relevant part of zi variables
    x(0) := zihh; x(1) := zilh; x(2) := zihl; x(3) := zill;
    y(0) := p0;   y(1) := p1;   y(2) := s0;   y(3) := t0;
    Add(x,y,z);
    zihh := z(0); zilh := z(1); zihl := z(2); zill := z(3);
  end Update_Product;

  procedure Inner_Product
              ( zrhh,zihh,zrlh,zilh : out double_float;
                zrhl,zihl,zrll,zill : out double_float;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                x,y,z : in Standard_Floating_Vectors.Link_to_Vector ) is

    dim : constant integer32 := xr'last/4;
    idx : integer32 := xr'first;

  begin
    zrhh := 0.0; zihh := 0.0; zrlh := 0.0; zilh := 0.0;
    zrhl := 0.0; zihl := 0.0; zrll := 0.0; zill := 0.0;
    for k in 1..dim loop
      Update_Product(zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill,
                     xr(idx),xi(idx),xr(idx+1),xi(idx+1),
                     xr(idx+2),xi(idx+2),xr(idx+3),xi(idx+3),
                     yr(idx),yi(idx),yr(idx+1),yi(idx+1),
                     yr(idx+2),yi(idx+2),yr(idx+3),yi(idx+3),x,y,z);
      idx := idx + 4;
    end loop;
  end Inner_Product;

  procedure Multiply
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                yr : in Standard_Floating_Vectors.Link_to_Vector;
                yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr : in Standard_Floating_Vectors.Link_to_Vector;
                zi : in Standard_Floating_Vectors.Link_to_Vector;
                x,y,z : in Standard_Floating_Vectors.Link_to_Vector ) is

    deg : constant integer32 := (xr'last+1)/4 - 1;
    zk,xk,yk : integer32;

  begin
   -- product(0) := first(0)*second(0);
    zr(0) := 0.0; zr(1) := 0.0; zr(2) := 0.0; zr(3) := 0.0;
    zi(0) := 0.0; zi(1) := 0.0; zi(2) := 0.0; zi(3) := 0.0;
    Update_Product(zr(0),zi(0),zr(1),zi(1),zr(2),zi(2),zr(3),zi(3),
                   xr(0),xi(0),xr(1),xi(1),xr(2),xi(2),xr(3),xi(3),
                   yr(0),yi(0),yr(1),yi(1),yr(2),yi(2),yr(3),yi(3),x,y,z);
    zk := 4;
    for k in 1..deg loop
     -- product(k) := first(0)*second(k);
      zr(zk) := 0.0; zr(zk+1) := 0.0; zr(zk+2) := 0.0; zr(zk+3) := 0.0;
      zi(zk) := 0.0; zi(zk+1) := 0.0; zi(zk+2) := 0.0; zi(zk+3) := 0.0;
      Update_Product(zr(zk),zi(zk),zr(zk+1),zi(zk+1),
                     zr(zk+2),zi(zk+2),zr(zk+3),zi(zk+3),
                     xr(0),xi(0),xr(1),xi(1),xr(2),xi(2),xr(3),xi(3),
                     yr(zk),yi(zk),yr(zk+1),yi(zk+1),
                     yr(zk+2),yi(zk+2),yr(zk+3),yi(zk+3),x,y,z);
      xk := 4; yk := zk-4;
      for i in 1..k loop
       -- product(k) := product(k) + first(i)*second(k-i);
        Update_Product(zr(zk),zi(zk),zr(zk+1),zi(zk+1),
                       zr(zk+2),zi(zk+2),zr(zk+3),zi(zk+3),
                       xr(xk),xi(xk),xr(xk+1),xi(xk+1),
                       xr(xk+2),xi(xk+2),xr(xk+3),xi(xk+3),
                       yr(yk),yi(yk),yr(yk+1),yi(yk+1),
                       yr(yk+2),yi(yk+2),yr(yk+3),yi(yk+3),x,y,z);
        xk := xk + 4;
        yk := yk - 4;
      end loop;
      zk := zk + 4;
    end loop;
  end Multiply;

end QuadDobl_Vector_Splitters;
