with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Standard_Inlined_BLAS_Helpers is

  function cdabs ( zr,zi : double_float ) return double_float is
  begin
    return SQRT(zr*zr + zi*zi);
  end cdabs;

  procedure csign ( pr,pi,qr,qi : in double_float;
                    zr,zi : out double_float ) is

    f : constant double_float := cdabs(pr,pi)/cdabs(qr,qi);

  begin
    zr := f*qr;
    zi := f*qi;
  end csign;

  function dznrm2 ( n : integer32;
                    xre : Standard_Floating_Vectors.Link_to_Vector;
                    xim : Standard_Floating_Vectors.Link_to_Vector;
                    ind,incx : integer32 ) return double_float is

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
        if xre(ix) /= 0.0 then
          temp := abs(xre(ix));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2;
            scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        if xim(ix) /= 0.0 then
          temp := abs(xim(ix));
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
  end dznrm2;

  function dznrm2 ( n : integer32;
                    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
                    row,col,incx : integer32 ) return double_float is

    ix : integer32;
    norm,scale,ssq,temp : double_float;
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if n < 1 or incx < 1 then
      norm := 0.0;
    else
      scale := 0.0; ssq := 1.0; ix := row;
      xre := rvv(col); xim := ivv(col);
      while ix <= row + (n-1)*incx loop
        if xre(ix) /= 0.0 then
          temp := abs(xre(ix));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2; scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        if xim(ix) /= 0.0 then
          temp := abs(xim(ix));
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
  end dznrm2;

  procedure zscal ( n : in integer32; zr,zi : in double_float;
                    xre : in Standard_Floating_Vectors.Link_to_Vector;
                    xim : in Standard_Floating_Vectors.Link_to_Vector;
                    ind,incx : in integer32 ) is

    ix : integer32;
    pr,pi : double_float;

  begin
    if n <= 0 or incx <= 0 then
      null;
    elsif incx = 1 then
      for i in 0..n-1 loop
        pr := xre(ind+i);
        pi := xim(ind+i);
        xre(ind+i) := zr*pr - zi*pi;
        xim(ind+i) := zr*pi + zi*pr;
      end loop;
    else
      ix := ind; 
      for i in 1..n loop
        pr := xre(ix);
        pi := xim(ix);
        xre(ix) := zr*pr - zi*pi;
        xim(ix) := zr*pi + zi*pr;
        ix := ix + incx;
      end loop;
    end if;
  end zscal;

  procedure zscal ( n : in integer32; zr,zi : in double_float;
                    rvv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    ivv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    row,col,incx : in integer32 ) is

    ix : integer32;
    pr,pi : double_float;
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if n <= 0 or incx <= 0 then
      null;
    elsif incx = 1 then
      xre := rvv(col);
      xim := ivv(col);
      for i in 0..n-1 loop
        pr := xre(row+i);
        pi := xim(row+i);
        xre(row+i) := zr*pr - zi*pi;
        xim(row+i) := zr*pi + zi*pr;
      end loop;
    else
      xre := rvv(col);
      xim := ivv(col);
      ix := row;
      for i in 1..n loop
        pr := xre(ix);
        pi := xim(ix);
        xre(ix) := zr*pr - zi*pi;
        xim(ix) := zr*pi + zi*pr;
        ix := ix + incx;
      end loop;
    end if;
  end zscal;

  procedure zaxpy ( n : in integer32; zr,zi : in double_float;
                    xre : in Standard_Floating_Vectors.Link_to_Vector;
                    xim : in Standard_Floating_Vectors.Link_to_Vector;
                    ind,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32 ) is

    ix,iy : integer32;
    abz : constant double_float := abs(zr) + abs(zi);
    yre,yim : Standard_Floating_Vectors.Link_to_Vector;
    pr,pi : double_float;

  begin
    if n > 0 and then abz /= 0.0 then
      yre := yrv(cly); yim := yiv(cly);
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
         -- y(rwy+i,cly) := y(rwy+i,cly) + z*x(ind+i);
          pr := xre(ind+i); pi := xim(ind+i);
          yre(rwy+i) := yre(rwy+i) + zr*pr - zi*pi;
          yim(rwy+i) := yim(rwy+i) + zr*pi + zi*pr;
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
         -- y(iy,cly) := y(iy,cly) + z*x(ix);
          pr := xre(ix); pi := xim(ix);
          yre(iy) := yre(iy) + zr*pr - zi*pi;
          yim(iy) := yim(iy) + zr*pi + zi*pr;
          iy := iy + incy;
          ix := ix + incx;
        end loop;
      end if;
    end if;
  end zaxpy;

  procedure zaxpy ( n : in integer32; zr,zi : in double_float;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yre : in  Standard_Floating_Vectors.Link_to_Vector;
                    yim : in  Standard_Floating_Vectors.Link_to_Vector;
                    ind,incy : in integer32 ) is

    ix,iy : integer32;
    abz : constant double_float := abs(zr) + abs(zi);
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    pr,pi : double_float;

  begin
    if n > 0 and then abz /= 0.0 then
      xre := xrv(clx); xim := xiv(clx);
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
         -- y(ind+i) := y(ind+i) + z*x(rwx+i,clx);
          pr := xre(rwx+i); pi := xim(rwx+i);
          yre(ind+i) := yre(ind+i) + zr*pr - zi*pi;
          yim(ind+i) := yim(ind+i) + zr*pi + zi*pr;
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
         -- y(iy) := y(iy) + z*x(ix,clx);
          pr := xre(ix); pi := xim(ix);
          yre(iy) := yre(iy) + zr*pr - zi*pi;
          yim(iy) := yim(iy) + zr*pi + zi*pr;
          iy := iy + incy;
          ix := ix + incx;
        end loop;
      end if;
    end if;
  end zaxpy;

  procedure zaxpy ( n : in integer32; zr,zi : in double_float;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32 ) is

    ix,iy : integer32;
    abz : constant double_float := abs(zr) + abs(zi);
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    yre,yim : Standard_Floating_Vectors.Link_to_Vector;
    pr,pi : double_float;

  begin
    if n > 0 and then abz /= 0.0 then
      xre := xrv(clx); xim := xiv(clx);
      yre := yrv(cly); yim := yiv(cly);
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
         -- y(rwy+i,cly) := y(rwy+i,cly) + z*x(rwx+i,clx);
          pr := xre(rwx+i); pi := xim(rwx+i);
          yre(rwy+i) := yre(rwy+i) + zr*pr - zi*pi;
          yim(rwy+i) := yim(rwy+i) + zr*pi + zi*pr;
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
         -- y(iy,cly) := y(iy,cly) + z*x(ix,clx);
          pr := xre(ix); pi := xim(ix);
          yre(iy) := yre(iy) + zr*pr - zi*pi;
          yre(iy) := yre(iy) + zr*pi + zi*pr;
          iy := iy + incy;
          ix := ix + incx;
        end loop;
      end if;
    end if;
  end zaxpy;

  procedure zdotc ( n : in integer32;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32;
                    zr,zi : out double_float ) is

    ix,iy : integer32;
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    yre,yim : Standard_Floating_Vectors.Link_to_Vector;
    pr,pi,qr,qi : double_float;

  begin
    zr := 0.0; zi := 0.0;
    xre := xrv(clx); xim := xiv(clx);
    yre := yrv(cly); yim := yiv(cly);
    if incx = 1 and incy = 1 then
      for i in 0..n-1 loop
       -- ztemp := ztemp + Conjugate(x(rwx+i,clx))*y(rwy+i,cly);
        pr := xre(rwx+i); pi := -xim(rwx+i); -- conjugate
        qr := yre(rwx+i); qi := yim(rwx+i);
        zr := zr + pr*qr - pi*qi;
        zi := zi + pr*qi + pi*qr;
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
       -- ztemp := ztemp + Conjugate(x(ix,clx))*y(iy,cly);
        pr := xre(ix); pi := -xim(ix); -- conjugate
        qr := yre(iy); qi := yim(iy);
        zr := zr + pr*qr - pi*qi;
        zi := zi + pr*qi + pi*qr;
        ix := ix + incx;
        iy := iy + incy;
      end loop;
    end if;
  end zdotc;

  procedure zdrot ( n : in integer32;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32;
                    c,s : in double_float ) is

    ix,iy : integer32;
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    yre,yim : Standard_Floating_Vectors.Link_to_Vector;
    pr,pi,qr,qi,zr,zi : double_float;

  begin
    if n > 0 then
      xre := xrv(clx); xim := xiv(clx);
      yre := yrv(cly); yim := yiv(cly);
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
         -- ztemp := c*x(rwx+i,clx) + s*y(rwy+i,cly);
          pr := xre(rwx+i); pi := xim(rwx+i);
          qr := yre(rwy+i); qi := yim(rwy+i);
          zr := c*pr + s*qr;
          zi := c*pi + s*qi;
         -- y(rwy+i,cly) := c*y(rwy+i,cly) - s*x(rwx+i,clx);
          yre(rwy+i) := c*qr - s*pr;
          yim(rwy+i) := c*qi - s*pi;
         -- x(rwx+i,clx) := ztemp;
          xre(rwx+i) := zr;
          xim(rwx+i) := zi;
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
         -- ztemp := c*x(ix,clx) + s*y(iy,cly);
          pr := xre(ix); pi := xim(ix);
          qr := yre(iy); qi := yim(iy);
          zr := c*pr + s*qr;
          zi := c*pi + s*qi;
         -- y(iy,cly) := c*y(iy,cly) - s*x(ix,clx);
          yre(iy) := c*qr - s*pr;
          yim(iy) := c*qi - s*pi;
         -- x(ix,clx) := ztemp;
          xre(ix) := zr;
          xim(ix) := zi;
          ix := ix + incx; iy := iy + incy;
        end loop;
      end if;
    end if;
  end zdrot;

  procedure zswap ( n : in integer32;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32 ) is

    ix,iy : integer32;
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    yre,yim : Standard_Floating_Vectors.Link_to_Vector;
    zr,zi : double_float;
 
  begin
    if n > 0 then
      xre := xrv(clx); xim := xiv(clx);
      yre := yrv(cly); yim := yiv(cly);
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
          zr := xre(rwx+i); zi := xim(rwx+i);   -- ztemp := x(rwx+i,clx);
          xre(rwx+i) := yre(rwy+i);
          xim(rwx+i) := yim(rwy+i);      -- x(rwx+i,clx) := y(rwy+i,cly);
          yre(rwy+i) := zr; yim(rwy+i) := zi;   -- y(rwy+i,cly) := ztemp;
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
          zr := xre(ix); zi := xim(ix);      -- ztemp := x(ix,clx);
          xre(ix) := yre(iy);
          xim(ix) := yim(iy);                -- x(ix,clx) := y(iy,cly);
          yre(iy) := zr; yim(iy) := zi;      -- y(iy,cly) := ztemp;
          ix := ix + incx; iy := iy + incy;
        end loop;
      end if;
    end if;
  end zswap;

end Standard_Inlined_BLAS_Helpers;
