with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;         use Standard_Complex_VecVecs_io;
with Standard_Random_Vectors;
with Standard_Vector_Splitters;           use Standard_Vector_Splitters;

procedure ts_perfade is

-- DESCRIPTION :
--   Tests better performing algorithmic differentiation and evaluation.

  procedure Forward ( x : in Standard_Complex_Vectors.Link_to_Vector;
                      f : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.

  -- REQUIRED : f'first = x'first = 1 and f'last >= x'last-1.

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
  end Forward;

  procedure Forward ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                      xi : in Standard_Floating_Vectors.Link_to_Vector;
                      fr : in Standard_Floating_Vectors.Link_to_Vector;
                      fi : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.

  -- REQUIRED :
  --   xr'range = xi'range, fr'first = xr'first = 1,
  --   and fi'last >= xi'last-1.

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
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.

  -- REQUIRED :
  --    f'first = x'first = 1 and f'last >= x'last-1,
  --    b'first = b'first = 1 and b'last >= x'last-2.

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

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.
  --   On single loop is applied.

  -- REQUIRED :
  --    f'first = x'first = 1 and f'last >= x'last-1.
  --    b'first = b'first = 1 and b'last >= x'last-2.

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

  procedure Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.

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

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.
  --   Applies loop fusion.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.

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
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.
  --   Computes all cross products of the values in x
  --   and stores the products in b.

  -- REQUIRED : x'last > 2, 
  --   f'first = x'first = 1 and f'last >= x'last-1,
  --   b'first = b'first = 1 and b'last >= x'last-2,
  --   c'first = c'first = 1 and c'last >= x'last-2.

  -- ON ENTRY :
  --   x        values for the variables;
  --   f        space allocated for forward products;
  --   b        space allocated for backward products;
  --   c        space allocated for cross products.

  -- ON RETURN : let n be x'last-1,
  --   f(n)     holds the product of all variables in x;
  --   f(n-1)   is the partial derivative of the product
  --            with respect to the last variable;
  --   b(n-2)   is the partial derivative of the product
  --            with respect to the first variable;
  --   c(k)     is the partial derivative of the product
  --            with respect to the (k+1)-th variable.

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
    if x'last > 2 then
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

  procedure Fused_Forward_Backward_Cross
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.
  --   Computes all cross products of the values in x
  --   and stores the products in b.
  --   Applies loop fusion.

  -- REQUIRED : x'last > 2, 
  --   f'first = x'first = 1 and f'last >= x'last-1,
  --   b'first = b'first = 1 and b'last >= x'last-2,
  --   c'first = c'first = 1 and c'last >= x'last-2.

  -- ON ENTRY :
  --   x        values for the variables;
  --   f        space allocated for forward products;
  --   b        space allocated for backward products;
  --   c        space allocated for cross products.

  -- ON RETURN : let n be x'last-1,
  --   f(n)     holds the product of all variables in x;
  --   f(n-1)   is the partial derivative of the product
  --            with respect to the last variable;
  --   b(n-2)   is the partial derivative of the product
  --            with respect to the first variable;
  --   c(k)     is the partial derivative of the product
  --            with respect to the (k+1)-th variable.

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
        lastend := x'last-3;
        firstend := lastend/2;
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

  procedure Forward_Backward_Cross
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.
  --   Computes all cross products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in cr and the imaginary parts in ci.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.
  --   cr'first = cr'first = 1, ci'last >= ci'last-2.

  -- ON ENTRY :
  --   xr       real parts of the values for the variables;
  --   xi       imaginary parts of the values for the variables;
  --   fr       space allocated for real parts of forward products;
  --   fi       space allocated for imaginary parts of forward products;
  --   br       space allocated for real parts of backward products;
  --   br       space allocated for imaginary parts of backward products;
  --   cr       space allocated for real parts of cross products;
  --   cr       space allocated for imaginary parts of cross products.

  -- ON RETURN : let n be x'last-1,
  --   fr(n)    the real part of the product of all variables in x;
  --   fi(n)    the imaginary part of the product of all variables in x;
  --   fr(n-1)  is the real part of the partial derivative of the product
  --            with respect to the last variable;
  --   fi(n-1)  is the imaginary part of the partial derivative
  --            of the product with respect to the last variable;
  --   br(n-2)  is the real part of the partial derivative
  --            of the product with respect to the first variable;
  --   bi(n-2)  is the imaginary part of the partial derivative
  --            of the product with respect to the first variable;
  --   cr(k)    is the real part of the partial derivative 
  --            of the product with respect to the (k+1)-th variable;
  --   ci(k)    is the imaginary part of the partial derivative 
  --            of the product with respect to the (k+1)-th variable.

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

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.
  --   Computes all cross products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in cr and the imaginary parts in ci.
  --   Applies loop fusion.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.
  --   cr'first = cr'first = 1, ci'last >= ci'last-2.

  -- ON ENTRY :
  --   xr       real parts of the values for the variables;
  --   xi       imaginary parts of the values for the variables;
  --   fr       space allocated for real parts of forward products;
  --   fi       space allocated for imaginary parts of forward products;
  --   br       space allocated for real parts of backward products;
  --   br       space allocated for imaginary parts of backward products;
  --   cr       space allocated for real parts of cross products;
  --   cr       space allocated for imaginary parts of cross products.

  -- ON RETURN : let n be x'last-1,
  --   fr(n)    the real part of the product of all variables in x;
  --   fi(n)    the imaginary part of the product of all variables in x;
  --   fr(n-1)  is the real part of the partial derivative of the product
  --            with respect to the last variable;
  --   fi(n-1)  is the imaginary part of the partial derivative
  --            of the product with respect to the last variable;
  --   br(n-2)  is the real part of the partial derivative
  --            of the product with respect to the first variable;
  --   bi(n-2)  is the imaginary part of the partial derivative
  --            of the product with respect to the first variable;
  --   cr(k)    is the real part of the partial derivative 
  --            of the product with respect to the (k+1)-th variable;
  --   ci(k)    is the imaginary part of the partial derivative 
  --            of the product with respect to the (k+1)-th variable.

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
        qr1 := xr(idx); qi1 := xr(idx);
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
        qr1 := xr(idx); qi1 := xr(idx);
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
        lastend := xr'last-3;
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
        plusidx := firstend+1;
       -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
        idx := plusidx-1;
        pr3 := fr(idx); pi3 := fi(idx);
        idx := xr'last-2-plusidx;
        qr3 := br(idx); qi3 := bi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
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
        qr1 := xr(idx); qi1 := xr(idx);
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
        qr1 := xr(idx); qi1 := xr(idx);
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

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of range mxe'range with space for
  --   complex vectors of range 1..mxe(k)-1, for k in mxe'range.

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

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Floating_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of range mxe'range with space for
  --   floating-point vectors of range 1..mxe(k)-1, for k in mxe'range.

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
                x : in Standard_Complex_Vectors.Link_to_Vector;
                pwt : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Computes the power table for the values of the variables in x.

  -- REQUIRED :
  --   mxe'range = x'range, pwt is allocated according to mxe,
  --   pwt'range = x'range and pwt(k)'range = 1..mxe(k)-1.

  -- ON ENTRY :
  --   mxe      highest exponents of the variables,
  --            mxe(k) is the highest exponent of the k-th variable
  --   x        values for all variables;
  --   pwt      allocated memory for all powers of the values in x.

  -- ON RETURN :
  --   pwt      power table, pwt(k)(i) equals x(k)**(i+1),
  --            for i in range 1..mxe(k)-1.

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

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Computes the power table for the values of the variables,
  --   with real parts in xr and imaginary parts in xi.

  -- REQUIRED :
  --   mxe'range = xr'range = xi'range,
  --   rpwt and ipwt are allocated according to mxe,
  --   rpwt'range = xr'range, rpwt(k)'range = 1..mxe(k)-1,
  --   ipwt'range = xr'range, ipwt(k)'range = 1..mxe(k)-1.

  -- ON ENTRY :
  --   mxe      highest exponents of the variables,
  --            mxe(k) is the highest exponent of the k-th variable
  --   xr       real parts of the values for all variables;
  --   xi       imaginary parts of the values for all variables;
  --   rpwt     allocated memory for the real parts of all powers
  --            of the values of the variables;
  --   ipwt     allocated memory for the imaginary parts of all powers
  --            of the values of the variables.

  -- ON RETURN :
  --   rpwt     real part of the power table,
  --            rpwt(k)(i) equals the real part of x(k)**(i+1),
  --            where x(k) is the complex value of the k-th variable,
  --            for i in range 1..mxe(k)-1;
  --   rpwt     imaginary part of the power table,
  --            rpwt(k)(i) equals the imaginary part of x(k)**(i+1),
  --            where x(k) is the complex value of the k-th variable,
  --            for i in range 1..mxe(k)-1.

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

  procedure Test_Forward ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward products.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    v : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Forward(x,f);
    put_line("the result : "); put_line(f);
    Forward(xr,xi,fr,fi);
    v := Make_Complex(fr,fi);
    put_line("recomputed : "); put_line(v);
  end Test_Forward;

  procedure Test_Forward_Backward ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward/backward products.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    v,w,v2,w2 : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Forward_Backward(x,f,b);
    Fused_Forward_Backward(x,f2,b2);
    Forward_Backward(xr,xi,fr,fi,br,bi);
    Fused_Forward_Backward(xr,xi,fr2,fi2,br2,bi2);
    v := Make_Complex(fr,fi); v2 := Make_Complex(fr2,fi2);
    w := Make_Complex(br,bi); w2 := Make_Complex(br2,bi2);
    put_line("the forward products : "); put_line(f);
    put_line("the forward products with loop fusion : "); put_line(f2);
    put_line("forward products recomputed : "); put_line(v);
    put_line("forward products recomputed with loop fusion : "); put_line(v2);
    put_line("the backward products : "); put_line(b);
    put_line("the backward products with loop fusion : "); put_line(b2);
    put_line("backward products recomputed : "); put_line(w);
    put_line("backward products recomputed with loop fusion : "); put_line(w2);
  end Test_Forward_Backward;

  procedure Test_Forward_Backward_Cross ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward/backward/cross products.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    c : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cc);
    c2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cc2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    cr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c);
    ci : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(c);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    cr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    ci2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(c2);
    u,v,w,u2,v2,w2 : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Forward_Backward_Cross(x,f,b,c);
    Fused_Forward_Backward_Cross(x,f2,b2,c2);
    Forward_Backward_Cross(xr,xi,fr,fi,br,bi,cr,ci);
    Fused_Forward_Backward_Cross(xr,xi,fr2,fi2,br2,bi2,cr2,ci2);
    u := Make_Complex(cr,ci); u2 := Make_Complex(cr2,ci2);
    v := Make_Complex(fr,fi); v2 := Make_Complex(fr2,fi2);
    w := Make_Complex(br,bi); w2 := Make_Complex(br2,bi2);
    put_line("the forward products : "); put_line(f);
    put_line("the forward products with loop fusion : "); put_line(f2);
    put_line("forward products recomputed : "); put_line(v);
    put_line("forward products recomputed with loop fusion : "); put_line(v2);
    put_line("the backward products : "); put_line(b);
    put_line("the backward products with loop fusion : "); put_line(b2);
    put_line("backward products recomputed : "); put_line(w);
    put_line("backward products recomputed with loop fusion : "); put_line(w2);
    put_line("the cross products : "); put_line(c);
    put_line("the cross products with loop fusion : "); put_line(c2);
    put_line("cross products recomputed : "); put_line(u);
    put_line("cross products recomputed wth loop fusion: "); put_line(u2);
  end Test_Forward_Backward_Cross;

  procedure Test_Power_Table ( dim,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Given the dimension dim and the highest power pwr,
  --   tests the computation of the power table for random values.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range) := Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range) := Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range) := Allocate(mxe);
    v : Standard_Complex_VecVecs.VecVec(x'range);

  begin
    Power_Table(mxe,x,pwt);
    put_line("The power table : "); put_line(pwt);
    Power_Table(mxe,xr,xi,rpwt,ipwt);
    v := Standard_Vector_Splitters.Make_Complex(rpwt,ipwt);
    put_line("The recomputed power table : "); put_line(v);
  end Test_Power_Table;

  procedure Timing_Test_Forward ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward product computations as freq
  --   on random vectors of dimension dim.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Forward(x,f);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward products");
    tstart(timer);
    for k in 1..frq loop
      Forward(xr,xi,fr,fi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward products");
  end Timing_Test_Forward;

  procedure Timing_Test_Forward_Backward ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward/backward product computations as freq
  --   on random vectors of dimension dim.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Forward_Backward(x,f,b);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward & backward products");
    tstart(timer);
    for k in 1..frq loop
      Fused_Forward_Backward(x,f2,b2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex with loop fusion");
    tstart(timer);
    for k in 1..frq loop
      Forward_Backward(xr,xi,fr,fi,br,bi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward & backward products");
    tstart(timer);
    for k in 1..frq loop
      Fused_Forward_Backward(xr,xi,fr2,fi2,br2,bi2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real with loop fusion");
  end Timing_Test_Forward_Backward;

  procedure Timing_Test_Forward_Backward_Cross ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward/backward/cross product computations as freq
  --   on random vectors of dimension dim.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    c : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cc);
    c2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cc2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    cr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c);
    ci : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(c);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    cr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    ci2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Forward_Backward_Cross(x,f,b,c);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward, backward, cross ");
    tstart(timer);
    for k in 1..frq loop
      Fused_Forward_Backward_Cross(x,f2,b2,c2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"fused complex forward, backward, cross");
    tstart(timer);
    for k in 1..frq loop
      Forward_Backward_Cross(xr,xi,fr,fi,br,bi,cr,ci);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward, backward, cross");
    tstart(timer);
    for k in 1..frq loop
      Fused_Forward_Backward_Cross(xr,xi,fr2,fi2,br2,bi2,cr2,ci2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"fused real forward, backward, cross");
  end Timing_Test_Forward_Backward_Cross;

  procedure Timing_Test_Power_Table ( dim,pwr,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Given the dimension dim, the highest power pwr,
  --   and the frequency frq, times the computation of the power table.

    timer : Timing_Widget;
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range) := Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range) := Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range) := Allocate(mxe);

  begin
    tstart(timer);
    for k in 1..frq loop
      Power_Table(mxe,x,pwt);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex power table");
    tstart(timer);
    for k in 1..frq loop
      Power_Table(mxe,xr,xi,rpwt,ipwt);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex power table");
  end Timing_Test_Power_Table;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension,
  --   the type of test, and then launches the test.

    dim,frq,pwr : integer32 := 0;
    ans,tst : character;

  begin
    new_line;
    put("Give the dimension of the vectors : "); get(dim);
    new_line;
    put_line("MENU for testing ADE :");
    put_line("  1. forward products");
    put_line("  2. forward and backward products");
    put_line("  3. forward, backward, and cross products");
    put_line("  4. power table");
    put("Type 1, 2, 3, or 4 to select the test : ");
    Ask_Alternative(tst,"1234");
    if tst = '4' then
      new_line;
      put("Give the highest power : "); get(pwr);
    end if;
    new_line;
    put("Interactive tests ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      case tst is
        when '1' => Test_Forward(dim);
        when '2' => Test_Forward_Backward(dim);
        when '3' => Test_Forward_Backward_Cross(dim);
        when '4' => Test_Power_Table(dim,pwr);
        when others => null;
      end case;
    else
      new_line;
      put("Give the frequency of the tests : "); get(frq);
      case tst is
        when '1' => Timing_Test_Forward(dim,frq);
        when '2' => Timing_Test_Forward_Backward(dim,frq);
        when '3' => Timing_Test_Forward_Backward_Cross(dim,frq);
        when '4' => Timing_Test_Power_Table(dim,pwr,frq);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_perfade;
