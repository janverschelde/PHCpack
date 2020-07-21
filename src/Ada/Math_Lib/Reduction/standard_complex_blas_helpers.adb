with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;

package body Standard_Complex_BLAS_Helpers is

  function Max0 ( a,b : integer32 ) return integer32 is
  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Max0;

  function dmax1 ( x1,x2 : double_float ) return double_float is
  begin
    if x1 > x2
     then return x1;
     else return x2;
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3 : double_float ) return double_float is
  begin
    if x1 > x2
     then return dmax1(x1,x3);
     else return dmax1(x2,x3);
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3,x4 : double_float ) return double_float is
  begin
    if x1 > x2
     then return dmax1(x1,x3,x4);
     else return dmax1(x2,x3,x4);
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3,x4,x5 : double_float ) return double_float is
  begin
    if x1 > x2
     then return dmax1(x1,x3,x4,x5);
     else return dmax1(x2,x3,x4,x5);
    end if;
  end dmax1;

  function cabs1 ( z : Complex_Number ) return double_float is
  begin
    return (abs(REAL_PART(z)) + abs(IMAG_PART(z)));
  end cabs1;

  function cdabs ( z : Complex_Number ) return double_float is
  begin
    return SQRT(REAL_PART(z)**2 + IMAG_PART(z)**2);
  end cdabs;

  function csign ( z1,z2 : Complex_Number ) return Complex_Number is
  begin
    return (Create(cdabs(z1)/cdabs(z2))*z2);
--  exception
--    when others => put_line("exception caught by csign");
--                   put("z1 = "); put(z1); new_line;
--                   put("z2 = "); put(z2); new_line;
--                   put("cdabs(z2) = "); put(cdabs(z2)); new_line;
--                   raise;
  end csign;

  function dsign ( a,b : double_float ) return double_float is
  begin
    if b >= 0.0
     then return abs(a);
     else return -abs(a);
    end if;
  end dsign;

  function dznrm2 ( n : integer32; x : Vector; ind,incx : integer32 )
                  return double_float is

    ix : integer32;
    norm,scale,ssq,temp : double_float;

  begin
    if n < 1 or incx < 1 then
      norm := 0.0;
    else
      scale := 0.0;
      ssq := 1.0;
      ix := ind;
      while ix <= ind + (n-1)*incx loop
        if REAL_PART(x(ix)) /= 0.0 then
          temp := abs(REAL_PART(x(ix)));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2;
            scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        if IMAG_PART(x(ix)) /= 0.0 then
          temp := abs(IMAG_PART(x(ix)));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2;
            scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        ix := ix + incx;
      end loop;
      norm := scale*SQRT(ssq);
    end if;
    return norm;
--  exception
--    when others => put_line("exception caught by 1st dznrm2"); raise;
  end dznrm2;

  function dznrm2 ( n : integer32; x : Matrix; row,col,incx : integer32 )
                  return double_float is

    ix : integer32;
    norm,scale,ssq,temp : double_float;

  begin
    if n < 1 or incx < 1 then
      norm := 0.0;
    else
      scale := 0.0; ssq := 1.0; ix := row;
      while ix <= row + (n-1)*incx loop
        if REAL_PART(x(ix,col)) /= 0.0 then
          temp := abs(REAL_PART(x(ix,col)));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2; scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        if IMAG_PART(x(ix,col)) /= 0.0 then
          temp := abs(IMAG_PART(x(ix,col)));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2; scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        ix := ix + incx;
      end loop;
      norm := scale*SQRT(ssq);
    end if;
    return norm;
--  exception
--    when others => put_line("exception caught by 2nd dznrm2"); raise;
  end dznrm2;

  procedure zscal ( n : in integer32; za : in Complex_Number;
                    zx : in out Vector; ind,incx : in integer32 ) is

    ix : integer32;

  begin
    if n <= 0 or incx <= 0 then
      null;
    elsif incx = 1 then
      for i in 0..n-1 loop
        zx(ind+i) := za*zx(ind+i);
      end loop;
    else
      ix := ind; 
      for i in 1..n loop
        zx(ix) := za*zx(ix);
        ix := ix + incx;
      end loop;
    end if;
--  exception
--    when others => put_line("exception caught by 1st zscal"); raise;
  end zscal;

  procedure zscal ( n : in integer32; za : in Complex_Number;
                    zx : in out Matrix; row,col,incx : in integer32 ) is

    ix : integer32;

  begin
    if n <= 0 or incx <= 0 then
      null;
    elsif incx = 1 then
      for i in 0..n-1 loop
        zx(row+i,col) := za*zx(row+i,col);
      end loop;
    else
      ix := row;
      for i in 1..n loop
        zx(ix,col) := za*zx(ix,col);
        ix := ix + incx;
      end loop;
    end if;
--  exception
--    when others => put_line("exception caught by 2nd zscal"); raise;
  end zscal;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Vector; ind,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

    ix,iy : integer32;

  begin
    if n > 0 and then AbsVal(z) /= 0.0 then
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
          y(rwy+i,cly) := y(rwy+i,cly) + z*x(ind+i);
        end loop;
      else
        if incx < 0
         then ix := (-n+1)*incx + ind;
         else ix := ind;
        end if;
        if incy < 0
         then iy := (-n+1)*incy + rwy;
         else iy := rwy;
        end if;
        for i in 0..n-1 loop
          y(iy,cly) := y(iy,cly) + z*x(ix);
          iy := iy + incy;
          ix := ix + incx;
        end loop;
      end if;
    end if;
--  exception
--    when others => put_line("exception caught by 1st zaxpy"); raise;
  end zaxpy;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Vector; ind,incy : in integer32 ) is

    ix,iy : integer32;

  begin
    if n > 0 and then AbsVal(z) /= 0.0 then
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
          y(ind+i) := y(ind+i) + z*x(rwx+i,clx);
        end loop;
      else
        if incx < 0
         then ix := (-n+1)*incx + rwx;
         else ix := rwx;
        end if;
        if incy < 0
         then iy := (-n+1)*incy + ind;
         else iy := ind;
        end if;
        for i in 0..n-1 loop
          y(iy) := y(iy) + z*x(ix,clx);
          iy := iy + incy;
          ix := ix + incx;
        end loop;
      end if;
    end if;
--  exception
--    when others => put_line("exception caught by 2nd zaxpy"); raise;
  end zaxpy;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

    ix,iy : integer32;

  begin
    if n > 0 and then AbsVal(z) /= 0.0 then
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
          y(rwy+i,cly) := y(rwy+i,cly) + z*x(rwx+i,clx);
        end loop;
      else
        if incx < 0 
         then ix := (-n+1)*incx + rwx;
         else ix := rwx;
        end if;
        if incy < 0
         then iy := (-n+1)*incy + rwy;
         else iy := rwy;
        end if;
        for i in 0..n-1 loop
          y(iy,cly) := y(iy,cly) + z*x(ix,clx);
          iy := iy + incy;
          ix := ix + incx;
        end loop;
      end if;
    end if;
--  exception
--    when others => put_line("exception caught by 3rd zaxpy"); raise;
  end zaxpy;

  function zdotc ( n : in integer32; x : in Matrix;
                   rwx,clx,incx : in integer32;
                   y : in Matrix; rwy,cly,incy : in integer32 )
                 return Complex_Number is

    ztemp : Complex_Number := Create(0.0);
    ix,iy : integer32;

  begin
    if incx = 1 and incy = 1 then
      for i in 0..n-1 loop
        ztemp := ztemp + Conjugate(x(rwx+i,clx))*y(rwy+i,cly);
      end loop;
    else 
      if incx < 0
       then ix := (-n+1)*incx + rwx;
       else ix := rwx;
      end if;
      if incy < 0
       then iy := (-n+1)*incy + rwy;
       else iy := rwy;
      end if;
      for i in 0..n-1 loop
        ztemp := ztemp + Conjugate(x(ix,clx))*y(iy,cly);
        ix := ix + incx;
        iy := iy + incy;
      end loop;
    end if;
    return ztemp;
--  exception
--    when others => put_line("exception caught by zdotc"); raise;
  end zdotc;

  procedure drotg ( da,db,c,s : in out double_float ) is

    roe,scale,r,z : double_float;

  begin
   -- put("in drotg da = "); put(da); new_line;
   -- put("         db = "); put(db); new_line;
    roe := db;
    if abs(da) > abs(db) then roe := da; end if;
    scale := abs(da) + abs(db);
    if 1.0 + scale = 1.0 then
      c := 1.0; s := 0.0;
      r := 0.0; z := 0.0;
    else
      r := scale*SQRT((da/scale)**2 + (db/scale)**2);
      r := dsign(1.0,roe)*r;
      c := da/r;
      s := db/r;
      z := 1.0;
      if abs(da) > abs(db) then z := s; end if;
      if abs(db) >= abs(da) and c /= 0.0
       then z := 1.0/c;
      end if;
    end if;
    da := r;
    db := z;
   -- put("out drotg c = "); put(c); new_line;
   -- put("          s = "); put(s); new_line;
   -- put("         da = "); put(da); new_line;
   -- put("         db = "); put(db); new_line;
--  exception
--    when others => put_line("exception caught by drotg"); raise;
  end drotg;

  procedure zdrot ( n : in integer32;
                    x : in out Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32;
                    c,s : in double_float ) is

    ix,iy : integer32;
    ztemp : Complex_Number;

  begin
    if n > 0 then
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
          ztemp := c*x(rwx+i,clx) + s*y(rwy+i,cly);
          y(rwy+i,cly) := c*y(rwy+i,cly) - s*x(rwx+i,clx);
          x(rwx+i,clx) := ztemp;
        end loop;
      else
        if incx < 0
         then ix := (-n+1)*incx + rwx;
         else ix := rwx;
        end if;
        if incy < 0
         then iy := (-n+1)*incy + rwy;
         else iy := rwy;
        end if;
        for i in 0..n-1 loop
          ztemp := c*x(ix,clx) + s*y(iy,cly);
          y(iy,cly) := c*y(iy,cly) - s*x(ix,clx);
          x(ix,clx) := ztemp;
          ix := ix + incx; iy := iy + incy;
        end loop;
      end if;
    end if;
--  exception
--    when others => put_line("exception caught by zdrot"); raise;
  end zdrot;

  procedure zswap ( n : in integer32;
                    x : in out Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

    ix,iy : integer32;
    ztemp : Complex_Number;
 
  begin
    if n > 0 then
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
          ztemp := x(rwx+i,clx);
          x(rwx+i,clx) := y(rwy+i,cly);
          y(rwy+i,cly) := ztemp;
        end loop;
      else 
        if incx < 0
         then ix := (-n+1)*incx + rwx;
         else ix := rwx;
        end if;
        if incy < 0
         then iy := (-n+1)*incy + rwy;
         else iy := rwy;
        end if;
        for i in 0..n-1 loop
          ztemp := x(ix,clx);
          x(ix,clx) := y(iy,cly);
          y(iy,cly) := ztemp;
          ix := ix + incx; iy := iy + incy;
        end loop;
      end if;
    end if;
--  exception
--    when others => put_line("exception caught by zswap"); raise;
  end zswap;

end Standard_Complex_BLAS_Helpers;
