with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Basics;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Renormalizations;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Floating_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;

procedure ts_perfqdvc is

-- DESCRIPTION :
--   Development of better performing computations with vectors of
--   complex numbers in quad double precision.

  procedure Split ( x : in Complex_Number;
                    rehihi,imhihi,relohi,imlohi : out double_float;
                    rehilo,imhilo,relolo,imlolo : out double_float ) is

  -- DESCRIPTION :
  --   Splits the complex number in quad double precision into 8 parts.

  -- ON ENTRY :
  --   x        a quad double complex number.

  -- ON RETURN :
  --   rehihi   highest double of the real part of x,
  --            or high_part(high_part(real_part(x)));
  --   imhihi   highest double of the imaginary part of x;
  --            or high_part(high_part(imag_part(x)));
  --   relohi   second highest double of the real part of x;
  --            or low_part(high_part(real_part(x)));
  --   imlohi   second highest double of the imaginary part of x;
  --            or low_part(high_part(imag_part(x)));
  --   rehilo   second lowest double of the real part of x;
  --            or high_part(low_part(real_part(x)));
  --   imhilo   second lowest double of the imaginary part of x;
  --            or high_part(low_part(imag_part(x)));
  --   relolo   lowest double of the real part of x;
  --            or low_part(low_part(imag_part(x)));
  --   imlolo   lowest double of the imaginary part of x;
  --            or low_part(low_part(real_part(x))).

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

  -- DESCRIPTION :
  --   Merges the 8 doubles into a quad double precision complex number.

  -- ON ENTRY :
  --   rehihi   highest double of the real part for x;
  --   imhihi   highest double of the imaginary part for x;
  --   relohi   second highest double of the real part for x;
  --   imlohi   second highest double of the imaginary part for x;
  --   rehilo   second lowest double of the real part for x;
  --   imhilo   second lowest double of the imaginary part for x;
  --   relolo   lowest double of the real part for x;
  --   imlolo   lowest double of the imaginary part for x.

  -- ON RETURN :
  --   x        a quad double complex number, with
  --            rehihi = high_part(high_part(real_part(x))),
  --            imhihi = high_part(high_part(imag_part(x))),
  --            relohi = low_part(high_part(real_part(x))),
  --            imlohi = low_part(high_part(imag_part(x))),
  --            rehilo = high_part(low_part(real_part(x))),
  --            imhilo = high_part(low_part(imag_part(x))),
  --            relolo = low_part(low_part(real_part(x))),
  --            imlolo = low_part(low_part(imag_part(x))).

    real_part : constant quad_double
              := Quad_Double_Numbers.create(rehihi,relohi,rehilo,relolo);
    imag_part : constant quad_double
              := Quad_Double_Numbers.create(imhihi,imlohi,imhilo,imlolo);

  begin
    x := QuadDobl_Complex_Numbers.Create(real_part,imag_part);
  end Merge;

  procedure Two_Split
              ( x : in QuadDobl_Complex_Vectors.Vector;
                xr,xi : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Splits the vector of quad double complex numbers into
  --   two vectors of doubles, one with the real and the other
  --   with the imaginary parts.

  -- REQUIRED : x'first = xr'first = xi'first
  --   and xr'last = xi'last = 4*x'last.

  -- ON ENTRY :
  --   x        a vector of quad double complex numbers.

  -- ON RETURN :
  --   xr       the real parts of the complex numbers in x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the real part of x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the real part of x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the real part of x(k),
  --            x(4*(k-1) + 4) = lowest double of the real part of x(k);
  --   xi       the imaginary parts of the complex numbers in x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the imag part of x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the imag part of x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the imag part of x(k),
  --            x(4*(k-1) + 4) = lowest double of the imag part of x(k).

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

  procedure Split ( x : in QuadDobl_Complex_Vectors.Vector;
                    xrhh,xihh : out Standard_Floating_Vectors.Vector;
                    xrlh,xilh : out Standard_Floating_Vectors.Vector;
                    xrhl,xihl : out Standard_Floating_Vectors.Vector;
                    xrll,xill : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Splits the vector x of quad double precision complex numbers
  --   into 8 vectors with real and imaginary parts,
  --   and for each part into four doubles.

  -- REQUIRED :
  --   x'range = xrhh'range = xihh'range = xrlh'range = xilh'range
  --           = xrhl'range = xihl'range = xrll'range = xill'range.

  -- ON ENTRY :
  --   x        a vector of quad double precision complex numbers.

  -- ON RETURN :
  --   xrhh     highest doubles of the real parts of x;
  --   xihh     highest doubles of the imaginary parts of x;
  --   xrlh     second highest doubles of the real parts of x;
  --   xilh     second highest doubles of the imaginary parts of x;
  --   xrhl     second lowest doubles of the real parts of x;
  --   xihl     second lowest doubles of the imaginary parts of x;
  --   xrll     lowest doubles of the real parts of x;
  --   xill     lowest doubles of the imaginary parts of x.

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

  -- DESCRIPTION :
  --   Merges the 8 vectors of doubles into one vector of
  --   quad double precision complex numbers.

  -- REQUIRED :
  --   x'range = xrhh'range = xihh'range = xrlh'range = xilh'range
  --           = xrhl'range = xihl'range = xrll'range = xill'range.

  -- ON ENTRY :
  --   xrhh     highest doubles of the real parts for x;
  --   xihh     highest doubles of the imaginary parts for x;
  --   xrlh     second highest doubles of the real parts for x;
  --   xilh     second highest doubles of the imaginary parts for x;
  --   xrhl     second lowest doubles of the real parts for x;
  --   xihl     second lowest doubles of the imaginary parts for x;
  --   xrll     lowest doubles of the real parts for x;
  --   xill     lowest doubles of the imaginary parts for x.

  -- ON RETURN :
  --   x        a vector of quad double precision complex numbers, with
  --            xrhh = high_part(high_part(real_part(x))),
  --            xihh = high_part(high_part(imag_part(x))),
  --            xrlh = low_part(high_part(real_part(x))),
  --            xilh = low_part(high_part(imag_part(x))),
  --            xrhl = high_part(low_part(real_part(x))),
  --            xihl = high_part(low_part(imag_part(x))),
  --            xrll = low_part(low_part(real_part(x))),
  --            xill = low_part(low_part(imag_part(x))).

  begin
    for k in x'range loop
      Merge(x(k),xrhh(k),xihh(k),xrlh(k),xilh(k),
                 xrhl(k),xihl(k),xrll(k),xill(k));
    end loop;
  end Merge;

  procedure Two_Merge
              ( x : out QuadDobl_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Merges the 2 vectors of doubles into one vector of
  --   quad double precision complex numbers.

  -- REQUIRED : x'first = xr'first = xi'first
  --   and xr'last = xi'last = 4*x'last.

  -- ON ENTRY :
  --   xr       the real parts of the complex numbers for x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the real part for x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the real part for x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the real part for x(k),
  --            x(4*(k-1) + 4) = lowest double of the real part for x(k);
  --   xi       the imaginary parts of the complex numbers for x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the imag part for x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the imag part for x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the imag part for x(k),
  --            x(4*(k-1) + 4) = lowest double of the imag part for x(k).

  -- ON RETURN :
  --   x        a vector of quad double complex numbers,
  --            with real parts taken xr and imaginary parts from xi.

    idx : integer32 := xr'last;

  begin
    for k in x'range loop
      Merge(x(k),xr(idx),xi(idx),xr(idx+1),xi(idx+1),
            xr(idx+2),xi(idx+2),xr(idx+3),xi(idx+3));
      idx := idx + 4;
    end loop;
  end Two_Merge;

  procedure Add ( x,y : in Standard_Floating_Vectors.Vector;
                  z : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Adds two quad double numbers, x and y,
  --   and assigns the sum to the quad double number z.
  --   The vector representation for quad doubles is assumed.

  -- REQUIRED : x'range = y'range = z'range = 0..3.

  -- ON ENTRY :
  --   x        a quad double number, with
  --            x(0) the highest double of x,
  --            x(1) the second highest double of x,
  --            x(2) the second lowest double of x,
  --            x(3) the lowest double of x;
  --   y        a quad double number, with
  --            y(0) the highest double of y,
  --            y(1) the second highest double of y,
  --            y(2) the second lowest double of y,
  --            y(3) the lowest double of y.

  -- ON RETURN :
  --   z        the sum of x + y, with
  --            z(0) the highest double of x + y,
  --            z(1) the second highest double of x + y,
  --            z(2) the second lowest double of x + y,
  --            z(3) the lowest double of x + y.

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

  -- DESCRIPTION :
  --   Adds two quad double numbers, x and y,
  --   and assigns the sum to the quad double number z.
  --   The vector representation for quad doubles is assumed.

  -- REQUIRED : x'range = y'range = z'range = 0..3.

  -- ON ENTRY :
  --   x        a quad double number, with
  --            x(offset) the highest double of x,
  --            x(offset+1) the second highest double of x,
  --            x(offset+2) the second lowest double of x,
  --            x(offset+3) the lowest double of x;
  --   y        a quad double number, with
  --            y(offset) the highest double of y,
  --            y(offset+1) the second highest double of y,
  --            y(offset+2) the second lowest double of y,
  --            y(offset+3) the lowest double of y.

  -- ON RETURN :
  --   z        the sum of x + y, with
  --            z(offset) the highest double of x + y,
  --            z(offset+1) the second highest double of x + y,
  --            z(offset+2) the second lowest double of x + y,
  --            z(offset+3) the lowest double of x + y.

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
                  x,y,z : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Adds two vectors of quad double complex numbers,
  --   in the 8-vector representation of doubles.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrhh     highest double of the real part of x;
  --   xihh     highest double of the imaginary part of x;
  --   xrlh     second highest double of the real part of x;
  --   xilh     second highest double of the imaginary part of x;
  --   xrhl     second lowest double of the real part of x;
  --   xihl     second lowest double of the imaginary part of x;
  --   xrll     lowest double of the real part of x;
  --   xill     lowest double of the imaginary part of x;
  --   yrhh     highest double of the real part of y;
  --   yihh     highest double of the imaginary part of y;
  --   yrlh     second highest double of the real part of y;
  --   yilh     second highest double of the imaginary part of y;
  --   yrhl     second lowest double of the real part of y;
  --   yihl     second lowest double of the imaginary part of y;
  --   yrll     lowest double of the real part of y;
  --   yill     lowest double of the imaginary part of y;
  --   x        a 4-vector work space of range 0..3;
  --   y        a 4-vector work space of range 0..3;
  --   z        a 4-vector work space of range 0..3.

  -- ON RETURN :
  --   zrhh     highest double of the real part of the sum;
  --   zihh     highest double of the imaginary part of the sum;
  --   zrlh     second highest double of the real part of the sum;
  --   zilh     second highest double of the imaginary part of the sum;
  --   zrhl     second lowest double of the real part of the sum;
  --   zihl     second lowest double of the imaginary part of the sum;
  --   zrll     lowest double of the real part of the sum;
  --   zill     lowest double of the imaginary part of the sum.

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

  -- DESCRIPTION :
  --   Adds two vectors of quad double complex numbers,
  --   in the 2-vector representation of doubles.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xr       real parts of the vector x;
  --   xi       imaginary parts of the vector x;
  --   yr       real parts of the vector y;
  --   yi       imaginary parts of the vector y.

  -- ON RETURN :
  --   zr       real parts of the sum;
  --   zi       imaginary parts of the sum.

    dim : constant integer32 := zr'last/4;
    idx : integer32 := zr'first; -- allow for zero start index

  begin
    for k in zr'first..dim loop
      Add(idx,xr,yr,zr);
      Add(idx,xi,yi,zi);
      idx := idx + 4;
    end loop;
  end Two_Add;

  function Inner_Product
             ( x,y : QuadDobl_Complex_Vectors.Vector )
             return Complex_Number is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors x and y,
  --   without taking complex conjugates.

  -- REQUIRED : x'range = y'range.

    res : Complex_Number := Create(integer(0));

  begin
    for k in x'range loop
      res := res + x(k)*y(k);
    end loop;
    return res;   
  end Inner_Product;

  procedure Update_Product
              ( zrhh,zihh,zrlh,zilh : in out double_float;
                zrhl,zihl,zrll,zill : in out double_float;
                xrhh,xihh,xrlh,xilh : in double_float;
                xrhl,xihl,xrll,xill : in double_float;
                yrhh,yihh,yrlh,yilh : in double_float;
                yrhl,yihl,yrll,yill : in double_float;
                x,y,z : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Adds the product of two quad double complex numbers to z,
  --   represented by 8 doubles.

  -- ON ENTRY :
  --   zrhh     current highest double in the real part of the inner product;
  --   zihh     current highest double in the imaginary part of the product;
  --   zrlh     current second highest double in the real part of the product;
  --   zilh     current second highest double in the imaginary part;
  --   zrlh     current second lowest double in the real part of the product;
  --   zilh     current second lowest double in the imaginary part;
  --   zrll     current lowest double in the real part of the inner product;
  --   zill     current lowest double in the imaginary part of the product;
  --   xrhh     highest double in the real part of the 1st number
  --   xihh     highest double in the imaginary part of the 1st number
  --   xrlh     second highest double in the real part of the 1st number
  --   xilh     second highest double in the imaginary part of the 1st number;
  --   xrlh     second lowest double in the real part of the 1st number
  --   xilh     second lowest double in the imaginary part of the 1st number;
  --   xrll     lowest double in the real part of the 1st number;
  --   xill     lowest double in the imaginary part of the 1st number;
  --   yrhh     highest double in the real part of the 2nd number;
  --   yihh     highest double in the imaginary part of the 2nd number;
  --   yrlh     second highest double in the real part of the 2nd number;
  --   yilh     second highest double in the imaginary part of the 2nd number;
  --   yrlh     second lowest double in the real part of the 2nd number
  --   yilh     second lowest double in the imaginary part of the 2nd number;
  --   yrll     lowest double in the real part of the 2nd number;
  --   yill     lowest double in the imaginary part of the 2nd number;
  --   x,y,z    vectors of range 0..3 as work space.

  -- ON RETURN :
  --   zrhh     updated highest double in the real part of the inner product;
  --   zihh     updated highest double in the imaginary part of the product;
  --   zrlh     updated second highest double in the real part of the product;
  --   zilh     updated second highest double in the imaginary part;
  --   zrlh     updated second lowest double in the real part of the product;
  --   zilh     updated second lowest double in the imaginary part;
  --   zrll     updated lowest double in the real part of the inner product;
  --   zill     updated lowest double in the imaginary part of the product.

   -- zhihi,zlohi,zhilo,zlolo : double_float;
    p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,q5 : double_float;
    p6,p7,p8,p9,q6,q7,q8,q9,r0,r1,t0,t1,s0,s1,s2 : double_float;

  begin
   -- zre = xre*yre - xim*yim
   -- (1) compute xre*yre
    Double_Double_Basics.two_prod(xrhh,yrhh,p0,q0);
    Double_Double_Basics.two_prod(xrhh,yrlh,p1,q1);
    Double_Double_Basics.two_prod(xrlh,yrhh,p2,q2);
    Double_Double_Basics.two_prod(xrhh,yrhl,p3,q3);
    Double_Double_Basics.two_prod(xrlh,yrlh,p4,q4);
    Double_Double_Basics.two_prod(xrhl,yrhh,p5,q5);
    Quad_Double_Renormalizations.three_sum(p1,p2,q0);
    Quad_Double_Renormalizations.three_sum(p2,q1,q2);
    Quad_Double_Renormalizations.three_sum(p3,p4,p5);
    Double_Double_Basics.two_sum(p2,p3,s0,t0);
    Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s2 := q2 + p5;
    Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s2 := s2 + (t0 + t1);
    Double_Double_Basics.two_prod(xrhh,yrll,p6,q6);
    Double_Double_Basics.two_prod(xrlh,yrhl,p7,q7);
    Double_Double_Basics.two_prod(xrhl,yrlh,p8,q8);
    Double_Double_Basics.two_prod(xrll,yrhh,p9,q9);
    Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    Double_Double_Basics.two_sum(q4,q5,q4,q5);
    Double_Double_Basics.two_sum(p6,p7,p6,p7);
    Double_Double_Basics.two_sum(p8,p9,p8,p9);
    Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t1 := t1 + (q3 + q5);
    Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r1 := r1 + (p7 + p9); 
    Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q4 := q4 + (t1 + r1); 
    Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t1 := t1 + q4;
    t1 := t1 + xrlh * yrll + xrhl * yrhl + xrll * yrlh
        + q6 + q7 + q8 + q9 + s2;
    Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
   -- zhihi := p0; zlohi := p1; zhilo := s0; zlolo := t0;
   -- (2) update relevant part of zr variables
    x(0) := zrhh; x(1) := zrlh; x(2) := zrhl; x(3) := zrll;
    y(0) := p0;   y(1) := p1;   y(2) := s0;   y(3) := t0;
    Add(x,y,z);
    zrhh := z(0); zrlh := z(1); zrhl := z(2); zrll := z(3);
   -- (3) compute xim*yim
    Double_Double_Basics.two_prod(xihh,yihh,p0,q0);
    Double_Double_Basics.two_prod(xihh,yilh,p1,q1);
    Double_Double_Basics.two_prod(xilh,yihh,p2,q2);
    Double_Double_Basics.two_prod(xihh,yihl,p3,q3);
    Double_Double_Basics.two_prod(xilh,yilh,p4,q4);
    Double_Double_Basics.two_prod(xihl,yihh,p5,q5);
    Quad_Double_Renormalizations.three_sum(p1,p2,q0);
    Quad_Double_Renormalizations.three_sum(p2,q1,q2);
    Quad_Double_Renormalizations.three_sum(p3,p4,p5);
    Double_Double_Basics.two_sum(p2,p3,s0,t0);
    Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s2 := q2 + p5;
    Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s2 := s2 + (t0 + t1);
    Double_Double_Basics.two_prod(xihh,yill,p6,q6);
    Double_Double_Basics.two_prod(xilh,yihl,p7,q7);
    Double_Double_Basics.two_prod(xihl,yilh,p8,q8);
    Double_Double_Basics.two_prod(xill,yihh,p9,q9);
    Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    Double_Double_Basics.two_sum(q4,q5,q4,q5);
    Double_Double_Basics.two_sum(p6,p7,p6,p7);
    Double_Double_Basics.two_sum(p8,p9,p8,p9);
    Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t1 := t1 + (q3 + q5);
    Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r1 := r1 + (p7 + p9); 
    Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q4 := q4 + (t1 + r1); 
    Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t1 := t1 + q4;
    t1 := t1 + xilh * yill + xihl * yihl + xill * yilh
        + q6 + q7 + q8 + q9 + s2;
    Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
   -- zhihi := p0; zlohi := p1; zhilo := s0; zlolo := t0;
   -- (4) update relevant part of zr variables
    x(0) := zrhh; x(1) := zrlh; x(2) := zrhl; x(3) := zrll;
    y(0) := -p0;  y(1) := -p1;  y(2) := -s0;  y(3) := -t0;
    Add(x,y,z);
    zrhh := z(0); zrlh := z(1); zrhl := z(2); zrll := z(3);
   -- zim = xre*yim + xim * yre
   -- (5) compute xre*yim
    Double_Double_Basics.two_prod(xrhh,yihh,p0,q0);
    Double_Double_Basics.two_prod(xrhh,yilh,p1,q1);
    Double_Double_Basics.two_prod(xrlh,yihh,p2,q2);
    Double_Double_Basics.two_prod(xrhh,yihl,p3,q3);
    Double_Double_Basics.two_prod(xrlh,yilh,p4,q4);
    Double_Double_Basics.two_prod(xrhl,yihh,p5,q5);
    Quad_Double_Renormalizations.three_sum(p1,p2,q0);
    Quad_Double_Renormalizations.three_sum(p2,q1,q2);
    Quad_Double_Renormalizations.three_sum(p3,p4,p5);
    Double_Double_Basics.two_sum(p2,p3,s0,t0);
    Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s2 := q2 + p5;
    Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s2 := s2 + (t0 + t1);
    Double_Double_Basics.two_prod(xrhh,yill,p6,q6);
    Double_Double_Basics.two_prod(xrlh,yihl,p7,q7);
    Double_Double_Basics.two_prod(xrhl,yilh,p8,q8);
    Double_Double_Basics.two_prod(xrll,yihh,p9,q9);
    Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    Double_Double_Basics.two_sum(q4,q5,q4,q5);
    Double_Double_Basics.two_sum(p6,p7,p6,p7);
    Double_Double_Basics.two_sum(p8,p9,p8,p9);
    Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t1 := t1 + (q3 + q5);
    Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r1 := r1 + (p7 + p9); 
    Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q4 := q4 + (t1 + r1); 
    Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t1 := t1 + q4;
    t1 := t1 + xrlh * yill + xrhl * yihl + xrll * yilh
        + q6 + q7 + q8 + q9 + s2;
    Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
   -- zhihi := p0; zlohi := p1; zhilo := s0; zlolo := t0;
   -- (6) update relevant part of zi variables
    x(0) := zihh; x(1) := zilh; x(2) := zihl; x(3) := zill;
    y(0) := p0;   y(1) := p1;   y(2) := s0;   y(3) := t0;
    Add(x,y,z);
    zihh := z(0); zilh := z(1); zihl := z(2); zill := z(3);
   -- (7) compute xim*yre
    Double_Double_Basics.two_prod(xihh,yrhh,p0,q0);
    Double_Double_Basics.two_prod(xihh,yrlh,p1,q1);
    Double_Double_Basics.two_prod(xilh,yrhh,p2,q2);
    Double_Double_Basics.two_prod(xihh,yrhl,p3,q3);
    Double_Double_Basics.two_prod(xilh,yrlh,p4,q4);
    Double_Double_Basics.two_prod(xihl,yrhh,p5,q5);
    Quad_Double_Renormalizations.three_sum(p1,p2,q0);
    Quad_Double_Renormalizations.three_sum(p2,q1,q2);
    Quad_Double_Renormalizations.three_sum(p3,p4,p5);
    Double_Double_Basics.two_sum(p2,p3,s0,t0);
    Double_Double_Basics.two_sum(q1,p4,s1,t1);
    s2 := q2 + p5;
    Double_Double_Basics.two_sum(s1,t0,s1,t0);
    s2 := s2 + (t0 + t1);
    Double_Double_Basics.two_prod(xihh,yrll,p6,q6);
    Double_Double_Basics.two_prod(xilh,yrhl,p7,q7);
    Double_Double_Basics.two_prod(xihl,yrlh,p8,q8);
    Double_Double_Basics.two_prod(xill,yrhh,p9,q9);
    Double_Double_Basics.two_sum(q0,q3,q0,q3); 
    Double_Double_Basics.two_sum(q4,q5,q4,q5);
    Double_Double_Basics.two_sum(p6,p7,p6,p7);
    Double_Double_Basics.two_sum(p8,p9,p8,p9);
    Double_Double_Basics.two_sum(q0,q4,t0,t1);
    t1 := t1 + (q3 + q5);
    Double_Double_Basics.two_sum(p6,p8,r0,r1); 
    r1 := r1 + (p7 + p9); 
    Double_Double_Basics.two_sum(t0,r0,q3,q4);
    q4 := q4 + (t1 + r1); 
    Double_Double_Basics.two_sum(q3,s1,t0,t1);
    t1 := t1 + q4;
    t1 := t1 + xilh * yrll + xihl * yrhl + xill * yrlh
        + q6 + q7 + q8 + q9 + s2;
    Quad_Double_Renormalizations.renorm5(p0,p1,s0,t0,t1);
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
                x,y,z : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the inner product of two quad double complex vectors,
  --   stores as sequences of doubles.

  -- REQUIRED : xr'range = xi'range = yr'range = yi'range.

  -- ON ENTRY :
  --   xr       real parts of the first vector in the inner product;
  --   xi       imaginary parts of the first vector in the inner product;
  --   yr       real parts of the second vector in the inner product;
  --   yi       imaginary parts of the second vector in the inner product;
  --   x,y,z    work space of range 0..3.

  -- ON RETURN :
  --   zrhh     highest double in the real part of the inner product;
  --   zihh     highest double in the imaginary part of the product;
  --   zrlh     second highest double in the real part of the inner product;
  --   zilh     second highest double in the imaginary part of the product;
  --   zrlh     second lowest double in the real part of the inner product;
  --   zilh     second lowest double in the imaginary part of the product;
  --   zrll     lowest double in the real part of the inner product;
  --   zill     lowest double in the imaginary part of the product.

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

  procedure Test_Add ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and tests their sum.

    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xrhh,xihh,xrlh,xilh : Standard_Floating_Vectors.Vector(1..dim);
    xrhl,xihl,xrll,xill : Standard_Floating_Vectors.Vector(1..dim);
    yrhh,yihh,yrlh,yilh : Standard_Floating_Vectors.Vector(1..dim);
    yrhl,yihl,yrll,yill : Standard_Floating_Vectors.Vector(1..dim);
    zrhh,zihh,zrlh,zilh : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    zrhl,zihl,zrll,zill : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    urhh,uihh,urlh,uilh : Standard_Floating_Vectors.Link_to_Vector;
    urhl,uihl,urll,uill : Standard_Floating_Vectors.Link_to_Vector;
    vrhh,vihh,vrlh,vilh : Standard_Floating_Vectors.Link_to_Vector;
    vrhl,vihl,vrll,vill : Standard_Floating_Vectors.Link_to_Vector;
    wrhh,wihh,wrlh,wilh : Standard_Floating_Vectors.Link_to_Vector;
    wrhl,wihl,wrll,will : Standard_Floating_Vectors.Link_to_Vector;
    xwrk,ywrk,zwrk : Standard_Floating_Vectors.Vector(0..3);

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    z1 := x + y;
    put_line("The sum of two random vectors :"); put_line(z1);
    Split(x,xrhh,xihh,xrlh,xilh,xrhl,xihl,xrll,xill);
    Split(y,yrhh,yihh,yrlh,yilh,yrhl,yihl,yrll,yill);
    urhh := new Standard_Floating_Vectors.Vector'(xrhh);
    uihh := new Standard_Floating_Vectors.Vector'(xihh);
    urlh := new Standard_Floating_Vectors.Vector'(xrlh);
    uilh := new Standard_Floating_Vectors.Vector'(xilh);
    urhl := new Standard_Floating_Vectors.Vector'(xrhl);
    uihl := new Standard_Floating_Vectors.Vector'(xihl);
    urll := new Standard_Floating_Vectors.Vector'(xrll);
    uill := new Standard_Floating_Vectors.Vector'(xill);
    vrhh := new Standard_Floating_Vectors.Vector'(yrhh);
    vihh := new Standard_Floating_Vectors.Vector'(yihh);
    vrlh := new Standard_Floating_Vectors.Vector'(yrlh);
    vilh := new Standard_Floating_Vectors.Vector'(yilh);
    vrhl := new Standard_Floating_Vectors.Vector'(yrhl);
    vihl := new Standard_Floating_Vectors.Vector'(yihl);
    vrll := new Standard_Floating_Vectors.Vector'(yrll);
    vill := new Standard_Floating_Vectors.Vector'(yill);
    wrhh := new Standard_Floating_Vectors.Vector'(zrhh);
    wihh := new Standard_Floating_Vectors.Vector'(zihh);
    wrlh := new Standard_Floating_Vectors.Vector'(zrlh);
    wilh := new Standard_Floating_Vectors.Vector'(zilh);
    wrhl := new Standard_Floating_Vectors.Vector'(zrhl);
    wihl := new Standard_Floating_Vectors.Vector'(zihl);
    wrll := new Standard_Floating_Vectors.Vector'(zrll);
    will := new Standard_Floating_Vectors.Vector'(zill);
    Add(wrhh,wihh,wrlh,wilh,wrhl,wihl,wrll,will,
        urhh,uihh,urlh,uilh,urhl,uihl,urll,uill,
        vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,xwrk,ywrk,zwrk);
    Merge(z2,wrhh.all,wihh.all,wrlh.all,wilh.all,
             wrhl.all,wihl.all,wrll.all,will.all);
    put_line("The recomputed sum :"); put_line(z2);
  end Test_Add;

  procedure Test_Two_Add ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and tests their sum.

    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zr,zi : constant Standard_Floating_Vectors.Vector(1..4*dim)
          := (1..4*dim => 0.0);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    z1 := x + y;
    put_line("The sum of two random vectors :"); put_line(z1);
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    Two_Add(wr,wi,ur,ui,vr,vi);
    Two_Merge(z2,wr.all,wi.all);
    put_line("The recomputed sum :"); put_line(z2);
  end Test_Two_Add;

  procedure Test_Inner_Product ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and tests their inner product.

    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : Complex_Number;
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill : double_float; 
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    xw,yw,zw : Standard_Floating_Vectors.Vector(0..3);

  begin
    z1 := Inner_Product(x,y);
    put_line("The inner product of two random vectors :");
    put(z1); new_line;
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    Inner_Product(zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill,
                  ur,ui,vr,vi,xw,yw,zw);
    Merge(z2,zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill);
    put_line("The recomputed inner product :");
    put(z2); new_line;
  end Test_Inner_Product;

  procedure Timing_Add ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and times the sum,
  --   for a frequency equal to the value of frq.

    timer : Timing_Widget;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xrhh,xihh,xrlh,xilh : Standard_Floating_Vectors.Vector(1..dim);
    xrhl,xihl,xrll,xill : Standard_Floating_Vectors.Vector(1..dim);
    yrhh,yihh,yrlh,yilh : Standard_Floating_Vectors.Vector(1..dim);
    yrhl,yihl,yrll,yill : Standard_Floating_Vectors.Vector(1..dim);
    zrhh,zihh,zrlh,zilh : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    zrhl,zihl,zrll,zill : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    urhh,uihh,urlh,uilh : Standard_Floating_Vectors.Link_to_Vector;
    urhl,uihl,urll,uill : Standard_Floating_Vectors.Link_to_Vector;
    vrhh,vihh,vrlh,vilh : Standard_Floating_Vectors.Link_to_Vector;
    vrhl,vihl,vrll,vill : Standard_Floating_Vectors.Link_to_Vector;
    wrhh,wihh,wrlh,wilh : Standard_Floating_Vectors.Link_to_Vector;
    wrhl,wihl,wrll,will : Standard_Floating_Vectors.Link_to_Vector;
    xwrk,ywrk,zwrk : Standard_Floating_Vectors.Vector(0..3);

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := x + y;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex addition");
    Split(x,xrhh,xihh,xrlh,xilh,xrhl,xihl,xrll,xill);
    Split(y,yrhh,yihh,yrlh,yilh,yrhl,yihl,yrll,yill);
    urhh := new Standard_Floating_Vectors.Vector'(xrhh);
    uihh := new Standard_Floating_Vectors.Vector'(xihh);
    urlh := new Standard_Floating_Vectors.Vector'(xrlh);
    uilh := new Standard_Floating_Vectors.Vector'(xilh);
    urhl := new Standard_Floating_Vectors.Vector'(xrhl);
    uihl := new Standard_Floating_Vectors.Vector'(xihl);
    urll := new Standard_Floating_Vectors.Vector'(xrll);
    uill := new Standard_Floating_Vectors.Vector'(xill);
    vrhh := new Standard_Floating_Vectors.Vector'(yrhh);
    vihh := new Standard_Floating_Vectors.Vector'(yihh);
    vrlh := new Standard_Floating_Vectors.Vector'(yrlh);
    vilh := new Standard_Floating_Vectors.Vector'(yilh);
    vrhl := new Standard_Floating_Vectors.Vector'(yrhl);
    vihl := new Standard_Floating_Vectors.Vector'(yihl);
    vrll := new Standard_Floating_Vectors.Vector'(yrll);
    vill := new Standard_Floating_Vectors.Vector'(yill);
    wrhh := new Standard_Floating_Vectors.Vector'(zrhh);
    wihh := new Standard_Floating_Vectors.Vector'(zihh);
    wrlh := new Standard_Floating_Vectors.Vector'(zrlh);
    wilh := new Standard_Floating_Vectors.Vector'(zilh);
    wrhl := new Standard_Floating_Vectors.Vector'(zrhl);
    wihl := new Standard_Floating_Vectors.Vector'(zihl);
    wrll := new Standard_Floating_Vectors.Vector'(zrll);
    will := new Standard_Floating_Vectors.Vector'(zill);
    tstart(timer);
    for k in 1..frq loop
      Add(wrhh,wihh,wrlh,wilh,wrhl,wihl,wrll,will,
          urhh,uihh,urlh,uilh,urhl,uihl,urll,uill,
          vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,xwrk,ywrk,zwrk);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"refitted addition");
    tstart(timer);
    for k in 1..frq loop
      Split(x,xrhh,xihh,xrlh,xilh,xrhl,xihl,xrll,xill);
      Split(y,yrhh,yihh,yrlh,yilh,yrhl,yihl,yrll,yill);
      Add(wrhh,wihh,wrlh,wilh,wrhl,wihl,wrll,will,
          urhh,uihh,urlh,uilh,urhl,uihl,urll,uill,
          vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,xwrk,ywrk,zwrk);
      Merge(z2,wrhh.all,wihh.all,wrlh.all,wilh.all,
               wrhl.all,wihl.all,wrll.all,will.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, add, merge");
  end Timing_Add;

  procedure Timing_Two_Add ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and times the sum,
  --   in the 2-vector representation.

    timer : Timing_Widget;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zr,zi : constant Standard_Floating_Vectors.Vector(1..4*dim)
          := (1..4*dim => 0.0);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := x + y;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex addition");
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    tstart(timer);
    for k in 1..frq loop
      Two_Add(wr,wi,ur,ui,vr,vi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"2-vector addition");
    tstart(timer);
    for k in 1..frq loop
      Two_Split(x,xr,xi); Two_Split(y,yr,yi);
      Two_Add(wr,wi,ur,ui,vr,vi);
      Two_Merge(z2,wr.all,wi.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, add, merge");
  end Timing_Two_Add;

  procedure Timing_Inner_Product ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and times the inner product,
  --   in the 2-vector representation.

    timer : Timing_Widget;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : Complex_Number;
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill : double_float; 
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    xw,yw,zw : Standard_Floating_Vectors.Vector(0..3);

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := Inner_Product(x,y);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex inner product");
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    tstart(timer);
    for k in 1..frq loop
      Inner_Product(zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill,
                    ur,ui,vr,vi,xw,yw,zw);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"2-vector inner product");
    tstart(timer);
    for k in 1..frq loop
      Two_Split(x,xr,xi); Two_Split(y,yr,yi);
      Inner_Product(zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill,
                    ur,ui,vr,vi,xw,yw,zw);
      Merge(z2,zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, inner product, merge");
  end Timing_Inner_Product;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of the vectors.

    dim,frq : integer32 := 0;
    ans,tst : character;

  begin
    new_line;
    put_line("MENU to test vector operations :");
    put_line("  1. test addition on 8-vector representation");
    put_line("  2. test addition on 2-vector representation");
    put_line("  3. test inner product");
    put("Type 1, 2, or 3 to select a test : ");
    Ask_Alternative(tst,"123");
    new_line;
    put("Give the dimension : "); get(dim);
    new_line;
    put("Interactive test ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then
      new_line;
      put("Give the frequency : "); get(frq);
    end if;
    if frq = 0 then
      case tst is
        when '1' => Test_Add(dim);
        when '2' => Test_Two_Add(dim);
        when '3' => Test_Inner_Product(dim);
        when others => null;
      end case;
    else
      case tst is
        when '1' => Timing_Add(dim,frq);
        when '2' => Timing_Two_Add(dim,frq);
        when '3' => Timing_Inner_Product(dim,frq);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_perfqdvc;
