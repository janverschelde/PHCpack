-- i/o for debugging purposes
-- with text_io,integer_io;                use text_io,integer_io;
-- with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
-- with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
-- with Multprec_Complex_Matrices_io;      use Multprec_Complex_Matrices_io;

with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Mathematical_Functions;   use Multprec_Mathematical_Functions;

package body Multprec_Complex_Singular_Values is

-- AUXILIARY BLAS ROUTINES :

  function Min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.
  --   Note that this should correspond to the fortran min0 function.

  begin
    if a < b 
     then return a;
     else return b;
    end if;
  end Min0;

  function Max0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the maximum of a and b.
  --   Note that this should correspond to the fortran max0 function. 

  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Max0;

  function dmax1 ( x1,x2 : Floating_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   Returns the maximum of x1 and x2.

    res : Floating_Number;

  begin
    if x1 > x2
     then Copy(x1,res);
     else Copy(x2,res);
    end if;
    return res;
  end dmax1;

  function dmax1 ( x1,x2,x3 : Floating_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   Returns the maximum of the three floating-point numbers.

  begin
    if x1 > x2
     then return dmax1(x1,x3);
     else return dmax1(x2,x3);
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3,x4 : Floating_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   Returns the maximum of the four floating-point numbers.

  begin
    if x1 > x2
     then return dmax1(x1,x3,x4);
     else return dmax1(x2,x3,x4);
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3,x4,x5 : Floating_Number )
                 return Floating_Number is

  -- DESCRIPTION :
  --   Returns the maximum of the five floating-point numbers.

  begin
    if x1 > x2
     then return dmax1(x1,x3,x4,x5);
     else return dmax1(x2,x3,x4,x5);
    end if;
  end dmax1;

  function cabs1 ( z : Complex_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of real and imaginary
  --   part of the complex number z.  Translation of
  --     complex*16 zdum
  --     double precision cabs1
  --     cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))

    res : Floating_Number;
    rpz : Floating_Number := REAL_PART(z);
    ipz : Floating_Number := IMAG_PART(z);
    acc : Floating_Number := AbsVal(rpz);

  begin
    Copy(acc,res); Clear(acc);
    acc := AbsVal(ipz);
    Add(res,acc); Clear(acc);
    Clear(rpz); Clear(ipz);
    return res;
 -- exception
 --   when others => put_line("exception raised in cabs1"); raise;
  end cabs1;

  function cdabs ( z : Complex_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   Computes the modulus of the complex number, let us hope this
  --   corresponds to the `cdabs' fortran function.

    res : Floating_Number;
    rpz : Floating_Number := REAL_PART(z);
    ipz : Floating_Number := IMAG_PART(z);
    sum : Floating_Number := rpz*rpz;
    acc : Floating_Number := ipz*ipz;

  begin
    Clear(rpz); Clear(ipz);
    Add(sum,acc); Clear(acc);
    res := SQRT(sum); Clear(sum);
    return res;
  end cdabs;  

  function csign ( z1,z2 : Complex_Number ) return Complex_Number is

  -- DESCRIPTION : translated from
  --       csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2)) 

    res : Complex_Number;
    az1 : Floating_Number := cdabs(z1);
    az2 : Floating_Number := cdabs(z2);
    acc : Floating_Number := az1/az2;

  begin
   -- put("in csign with z1 = "); put(z1); new_line;
   -- put("in csign with z2 = "); put(z2); new_line;
    res := Create(acc);
    Mul(res,z2);
    Clear(az1); Clear(az2); Clear(acc);
    return res;
 -- exception 
 --   when others => put_line("exception raised in csign"); raise;
  end csign;

  function dsign ( a,b : Floating_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   The implementation of this routine is written from web page
  --   documentation of sign...

    res : Floating_Number := AbsVal(a);

  begin
    if b < 0.0
     then Min(res);
    end if;
    return res;
  end dsign;

  function dznrm2 ( n : integer32; x : Vector; ind,incx : integer32 )
                  return Floating_Number is

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector x, starting at x(ind)
  --   and continueing n steps with increment incx.

    ix : integer32;
    norm,scale,ssq,temp,wrk,acc: Floating_Number;

  begin
    if n < 1 or incx < 1 then
      norm := Create(0.0);
    else
      scale := Create(0.0);
      ssq := Create(1.0);
      ix := ind;
      while ix <= ind + (n-1)*incx loop
        wrk := REAL_PART(x(ix));
        if not Equal(wrk,0.0) then
          temp := AbsVal(wrk);
          if scale < temp then -- ssq := 1.0 + ssq*(scale/temp)**2;
            acc := scale/temp;
            Mul(acc,acc);
            Mul(acc,ssq);
            Clear(ssq);
            ssq := 1.0 + acc;
            Clear(acc);
            Copy(temp,scale);
          else -- ssq := ssq + (temp/scale)**2;
            acc := temp/scale;
            Mul(acc,acc);
            Add(ssq,acc);
            Clear(acc);
          end if;
          Clear(temp);
        end if;
        Clear(wrk);
        wrk := IMAG_PART(x(ix));
        if not Equal(wrk,0.0) then
          temp := AbsVal(wrk);
          if scale < temp then -- ssq := 1.0 + ssq*(scale/temp)**2;
            acc := scale/temp;
            Mul(acc,acc);
            Mul(acc,ssq);
            Clear(ssq);
            ssq := 1.0 + acc;
            Clear(acc);
            Copy(temp,scale);
          else -- ssq := ssq + (temp/scale)**2;
            acc := temp/scale;
            Mul(acc,acc);
            Add(ssq,acc);
            Clear(acc);
          end if;
          Clear(temp);
        end if;
        Clear(wrk);
        ix := ix + incx;
      end loop;
      -- norm := scale*SQRT(ssq);
      norm := SQRT(ssq);
      Mul(norm,scale);
      Clear(scale); Clear(ssq);
    end if;
    return norm;
 -- exception 
 --   when others => put_line("Exception raised in dznrm2"); raise;
  end dznrm2;

  function dznrm2 ( n : integer32; x : Matrix; row,col,incx : integer32 )
                  return Floating_Number is

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector x, starting at x(row,col)
  --   and continueing n steps with increment incx in the same column.

    ix : integer32;
    norm,scale,ssq,temp,wrk,acc : Floating_Number;

  begin
    if n < 1 or incx < 1 then
      norm := Create(0.0);
    else
      scale := Create(0.0);
      ssq := Create(1.0);
      ix := row;
      while ix <= row + (n-1)*incx loop
        wrk := REAL_PART(x(ix,col));
        if not Equal(wrk,0.0) then
          temp := AbsVal(wrk);
          if scale < temp then -- ssq := 1.0 + ssq*(scale/temp)**2;
            acc := scale/temp;
            Mul(acc,acc);
            Mul(acc,ssq);
            Clear(ssq);
            ssq := 1.0 + acc;
            Clear(acc);
            Copy(temp,scale);
          else -- ssq := ssq + (temp/scale)**2;
            acc := temp/scale;
            Mul(acc,acc);
            Add(ssq,acc);
            Clear(acc);
          end if;
          Clear(temp);
        end if;
        Clear(wrk);
        wrk := IMAG_PART(x(ix,col));
        if not Equal(wrk,0.0) then
          temp := AbsVal(wrk);
          if scale < temp then -- ssq := 1.0 + ssq*(scale/temp)**2;
            acc := scale/temp;
            Mul(acc,acc);
            Mul(acc,ssq);
            Clear(ssq);
            ssq := 1.0 + acc;
            Clear(acc);
            Copy(temp,scale);
          else -- ssq := ssq + (temp/scale)**2;
            acc := temp/scale;
            Mul(acc,acc);
            Add(ssq,acc);
            Clear(acc);
          end if;
          Clear(temp);
        end if;
        Clear(wrk);
        ix := ix + incx;
      end loop;
      norm := scale*SQRT(ssq);
    end if;
    return norm;
  end dznrm2;

  procedure zscal ( n : in integer32; za : in Complex_Number;
                    zx : in out Vector; ind,incx : in integer32 ) is

  -- DESCRIPTION :
  --   Scales the vector starting at zx(ind) with the constant za.

    ix : integer32;

  begin
    if n <= 0 or incx <= 0 then
      null;
    elsif incx = 1 then
      for i in 0..n-1 loop
        Mul(zx(ind+i),za);
      end loop;
    else
      ix := ind; 
      for i in 1..n loop
        Mul(zx(ix),za);
        ix := ix + incx;
      end loop;
    end if;
  end zscal;

  procedure zscal ( n : in integer32; za : in Complex_Number;
                    zx : in out Matrix; row,col,incx : in integer32 ) is

  -- DESCRIPTION :
  --   Scales the vector starting at zx(row,col) with the constant za.

    ix : integer32;

  begin
    if n <= 0 or incx <= 0 then
      null;
    elsif incx = 1 then
      for i in 0..n-1 loop
        Mul(zx(row+i,col),za);
      end loop;
    else 
      ix := row; 
      for i in 1..n loop
        Mul(zx(ix,col),za);
        ix := ix + incx;
      end loop;
    end if;
  end zscal;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Vector; ind,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at ind and in y
  --   at (rwy,cly), using increments incx and incy to advance.

    ix,iy : integer32;
    avz : Floating_Number;
    acc : Complex_Number;

  begin
    if n > 0 then
      avz := AbsVal(z);
      if avz > 0.0 then  -- if not Equal(avz,0.0)
        if incx = 1 and incy = 1 then
          for i in 0..n-1 loop
             -- y(rwy+i,cly) := y(rwy+i,cly) + z*x(ind+i);
            acc := z*x(ind+i);
            Add(y(rwy+i,cly),acc);
            Clear(acc);
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
            acc := z*x(ix);
            Add(y(iy,cly),acc);
            Clear(acc);
            iy := iy + incy;
            ix := ix + incx;
          end loop;
        end if;
      end if;
      Clear(avz);
    end if;
  end zaxpy;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Vector; ind,incy : in integer32 ) is

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at (rwx,clx)
  --   and in y at ind, using increments incx and incy to advance.

    ix,iy : integer32;
    avz : Floating_Number;
    acc : Complex_Number;

  begin
    if n > 0 then
      avz := AbsVal(z);
      if avz > 0.0 then -- if not Equal(avz,0.0) then
        if incx = 1 and incy = 1 then
          for i in 0..n-1 loop
            -- y(ind+i) := y(ind+i) + z*x(rwx+i,clx);
            acc := z*x(rwx+i,clx);
            Add(y(ind+i),acc);
            Clear(acc);
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
            acc := z*x(ix,clx);
            Add(y(iy),acc);
            Clear(acc);
            iy := iy + incy;
            ix := ix + incx;
          end loop;
        end if;
      end if;
      Clear(avz);
    end if;
  end zaxpy;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x and y at the
  --   respective (rows,columns): (rwx,clx) and (rwy,cly), using
  --   increments incx and incy to advance in the rows.

    ix,iy : integer32;
    avz : Floating_Number;
    acc : Complex_Number;

  begin
    if n > 0 then
       avz := AbsVal(z);
       if avz > 0.0 then -- if not Equal(avz,0.0)
         if incx = 1 and incy = 1 then
           for i in 0..n-1 loop
             -- y(rwy+i,cly) := y(rwy+i,cly) + z*x(rwx+i,clx);
             acc := z*x(rwx+i,clx);
             Add(y(rwy+i,cly),acc);
             Clear(acc);
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
             acc := z*x(ix,clx);
             Add(y(iy,cly),acc);
             Clear(acc);
             iy := iy + incy;
             ix := ix + incx;
           end loop;
         end if;
       end if;
       Clear(avz);
    end if;
  end zaxpy;

  function zdotc ( n : in integer32; x : in Matrix;
                   rwx,clx,incx : in integer32;
                   y : in Matrix; rwy,cly,incy : in integer32 )
                 return Complex_Number is

  -- DESCRIPTION :
  --   Returns the dot product of two vectors in two columns of matrices
  --   x and y, starting at rows rwx and rwy respectively, with respective
  --   increments in incx and incy.

    ztemp : Complex_Number := Multprec_Complex_Numbers.Create(integer(0));
    wrk : Complex_Number;
    ix,iy : integer32;

  begin
    if incx = 1 and incy = 1 then
      for i in 0..n-1 loop
        -- ztemp := ztemp + Conjugate(x(rwx+i,clx))*y(rwy+i,cly);
        wrk := Conjugate(x(rwx+i,clx));
        Mul(wrk,y(rwy+i,cly));
        Add(ztemp,wrk);
        Clear(wrk);
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
         wrk := Conjugate(x(ix,clx));
         Mul(wrk,y(iy,cly));
         Add(ztemp,wrk);
         Clear(wrk);
         ix := ix + incx;
         iy := iy + incy;
       end loop;
    end if;
    return ztemp;
  end zdotc;

  procedure drotg ( da,db,c,s : in out Floating_Number ) is

  -- DESCRIPTION :
  --   Constructs Givens plane rotation.

    roe,scale,r,z,wrk,acc,one : Floating_Number;
    absda : Floating_Number := AbsVal(da);
    absdb : Floating_Number := AbsVal(db);

  begin
   -- put("in drotg da = "); put(da); new_line;
   -- put("         db = "); put(db); new_line;
    Copy(db,roe);
    if absda > absdb
     then Copy(da,roe);
    end if;
    scale := absda + absdb;
    if Equal(scale,0.0) then
      c := Create(1.0); s := Create(0.0);
      r := Create(0.0); z := Create(0.0);
    else -- r := scale*SQRT((da/scale)**2 + (db/scale)**2);
      wrk := da/scale;
      Mul(wrk,wrk);
      acc := db/scale;
      Mul(acc,acc);
      Add(wrk,acc);
      Clear(acc);
      r := SQRT(wrk);
      Clear(wrk);
      Mul(r,scale);
      -- r := dsign(Create(1.0),roe)*r;
      one := Create(1.0);
      wrk := dsign(one,roe);
      Mul(r,wrk);
      Clear(wrk);
      c := da/r;
      s := db/r;
      if absda > absdb then
        Copy(s,z);
      else
        Copy(one,z);
        if not Equal(c,0.0)
         then Div(z,c);
        end if;
      end if;
      Clear(roe);
      Clear(one);
    end if;
    Clear(absda); Clear(absdb);
    Clear(scale);
    Clear(da); da := r;
    Clear(db); db := z;
   -- put("out drotg c = "); put(c); new_line;
   -- put("          s = "); put(s); new_line;
   -- put("         da = "); put(da); new_line;
   -- put("         db = "); put(db); new_line;
  end drotg;

  procedure zdrot ( n : in integer32;
                    x : in out Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32;
                    c,s : in Floating_Number ) is

  -- DESCRIPTION :
  --   Applies a plane rotation where the cos and sin are c and s
  --   and the vectors are in the columns of the matrices x and y,
  --   starting at (rwx,clx) and (rwy,cly) advancing in the rows with
  --   increments incx and incy respectively.

    ix,iy : integer32;
    ztemp,acc : Complex_Number;

  begin
    if n > 0 then
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop
          ztemp := c*x(rwx+i,clx);
          acc := s*y(rwy+i,cly);
          Add(ztemp,acc); Clear(acc);
          Mul(y(rwy+i,cly),c);
          acc := s*x(rwx+i,clx);
          Sub(y(rwy+i,cly),acc); Clear(acc);
          Copy(ztemp,x(rwx+i,clx)); Clear(ztemp);
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
          ztemp := c*x(ix,clx);
          acc := s*y(iy,cly);
          Add(ztemp,acc); Clear(acc);
          Mul(y(iy,cly),c);
          acc := s*x(ix,clx);
          Sub(y(iy,cly),acc); Clear(acc);
          Copy(ztemp,x(ix,clx)); Clear(ztemp);
          ix := ix + incx; iy := iy + incy;
        end loop;
      end if;
    end if;
  end zdrot;

  procedure zswap ( n : in integer32;
                    x : in out Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

  -- DESCRIPTION :
  --   Interchanges two vectors in the columns of the matrices x and y,
  --   respectively starting at (rwx,clx) and (rwy,cly), and advancing
  --   in the rows with the respective increments incx and incy.

    ix,iy : integer32;
    ztemp : Complex_Number;
 
  begin
    if n > 0 then
      if incx = 1 and incy = 1 then
        for i in 0..n-1 loop 
          Copy(x(rwx+i,clx),ztemp);
          Copy(y(rwy+i,cly),x(rwx+i,clx));
          Copy(ztemp,y(rwy+i,cly));
          Clear(ztemp);
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
          Copy(x(ix,clx),ztemp);
          Copy(y(iy,cly),x(ix,clx));
          Copy(ztemp,y(iy,cly));
          ix := ix + incx; iy := iy + incy;
        end loop;
      end if;
    end if;
  end zswap;

-- TARGET PROCEDURE :

  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix; 
                  job : in integer32; info : out integer32 ) is

    work : Vector(1..n);
    iter,jobu,kase,kk,ll,lm1,lp1,ls,lu,m,maxit : integer32;
    mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1 : integer32;
    t,r : Complex_Number;
    b,c,cs,el,emm1,f,g,scale,shift,sl,sm,sn : Floating_Number;
    smm1,t1,test,ztest,wrk,tmp : Floating_Number;
    acc,acc2 : Complex_Number;
    wantu,wantv : boolean;
    one : Floating_Number := Create(1.0);
    c_one : Complex_Number := Create(integer(1));
    f_min_one : Floating_Number := Create(-1.0);
    c_min_one : Complex_Number := Create(f_min_one);
    leave_loop : boolean;

  begin
   -- put_line("The matrix x on entry : "); put(x);
    maxit := 30;     -- maximal number of iterations
    wantu := false;  -- determine what is to be computed
    wantv := false;
    jobu := (job mod 100)/10;
    ncu := n;
    if jobu > 1 then ncu := min0(n,p); end if;
    if jobu /= 0 then wantu := true; end if;
    if job mod 10 /= 0 then wantv := true; end if;
   -- put_line("reduce x to bidiagonal form, storing the diagonal" 
   --          & " elements in s");
   -- put_line("and the super diagonal elements in e");
    info := 0;
    nct := min0(n-1,p);
    nrt := max0(0,min0(p-2,n));
    lu := max0(nct,nrt);
   -- put("lu = "); put(lu,1); new_line;
    if lu >= 1 then
      for l in 1..lu loop
       -- put("l = "); put(l,1); put("  nct = "); put(nct,1); new_line;
        lp1 := l+1;
        if l <= nct then
          --put_line("compute the transformation for the l-th column");
          --put_line("and place the l-th diagonal in s(l)");
          wrk := dznrm2(n-l+1,x,l,l,1);
          Clear(s(l));
          s(l) := Create(wrk);
          Clear(wrk);
          -- put("s("); put(l,1); put(") = "); put(s(l)); new_line;
          wrk := cabs1(s(l));
          -- put("wrk = "); put(wrk); new_line;
          if wrk > 0.0 then -- if not Equal(wrk,0.0)
            tmp := cdabs(x(l,l));
            -- put("tmp = "); put(tmp); new_line;
            if tmp > 0.0 then -- if not Equal(tmp,0.0)
              --put_line("Calling csign");
              acc := csign(s(l),x(l,l));
              Copy(acc,s(l));
              Clear(acc);
            end if;
            -- put_line("clearing tmp");
            Clear(tmp);
            -- put_line("Moving into Div 1...");
            -- declare
            -- begin
            --    put("s("); put(l,1); put(") = "); put(s(l)); new_line;
            acc := c_one/s(l);
            --    put("inverse : "); put(acc); new_line;
            -- exception
            --    when others => put_line("exception Div 1"); raise;
            -- end;
            -- put_line("Ending Div 1.");
            zscal(n-l+1,acc,x,l,l,1);
            Clear(acc);
            Add(x(l,l),one);
            -- put("x("); put(l,1); put(","); put(l,1);
            -- put(") = "); put(x(l,l)); new_line;
          end if;
          Clear(wrk);
          Min(s(l));
        end if;
        -- put_line("The matrix x : "); put(x);
        -- put("p = "); put(p,1); put("  lp1 = "); put(lp1,1); new_line;
        if p >= lp1 then
          for j in lp1..p loop
            if l <= nct then
              wrk := cabs1(s(l));
              if wrk > 0.0 then -- if (not Equal(wrk,0.0))
                -- apply the transformation
                -- t := -zdotc(n-l+1,x,l,l,1,x,l,j,1)/x(l,l);
                t := zdotc(n-l+1,x,l,l,1,x,l,j,1);
                Min(t);
                -- declare
                -- begin
                --    put("x("); put(l,1); put(","); put(l,1);
                --    put(") : "); put(x(l,l)); new_line;
                Div(t,x(l,l));
                --    put(" t = "); put(t); new_line;
                -- exception
                --    when others => put_line("Div 2 exception");
                --                   raise;
                -- end;
                zaxpy(n-l+1,t,x,l,l,1,x,l,j,1);
                Clear(t);
              end if;
              Clear(wrk);
            end if;
            -- place the l-th row of x into e for the subsequent
            -- calculation of the row transformation
            acc := Conjugate(x(l,j));
            Copy(acc,e(j));
            Clear(acc);
          end loop;
        end if;
        -- put_line("The matrix x : "); put(x);
        if wantu and l <= nct then
          -- place the transformation in u for subsequent
          -- back multiplication
          for i in l..n loop
            Copy(x(i,l),u(i,l));
          end loop;
        end if;
        if l <= nrt then
          -- compute the l-th row transformation 
          -- and place the l-th super diagonal in e(l)
          wrk := dznrm2(p-l,e,lp1,1);
          Clear(e(l));
          e(l) := Create(wrk);
          Clear(wrk);
          wrk := cabs1(e(l));
          if wrk > 0.0 then -- if not Equal(wrk,0.0)
            tmp := cdabs(e(lp1));
            if tmp > 0.0 then -- if not Equal(tmp,0.0)
              acc := csign(e(l),e(lp1));
              Copy(acc,e(l));
              Clear(acc);
            end if;
            Clear(tmp);
            -- put("e("); put(l,1); put(") = "); put(e(l)); new_line;
            acc := c_one/e(l);
            -- put("inverse : "); put(acc); new_line;
            zscal(p-l,acc,e,lp1,1);
            Clear(acc);
            Add(e(lp1),one);
          end if;
          Clear(wrk);
          acc := Conjugate(e(l));
          Min(acc);
          Copy(acc,e(l));
          Clear(acc);
          wrk := cabs1(e(l));
          if lp1 <= n and (wrk > 0.0) then -- not Equal(wrk,0.0)
            -- apply the transformation
            for i in lp1..n loop
              work(i) := Create(integer(0));
            end loop;
            for j in lp1..p loop
              zaxpy(n-l,e(j),x,lp1,j,1,work,lp1,1);
            end loop;
            for j in lp1..p loop
              acc := -e(j);
              -- put("e("); put(lp1,1); put(") = ");
              -- put(e(lp1)); new_line;
              Div(acc,e(lp1));
              -- put("acc : "); put(acc); new_line;
              acc2 := Conjugate(acc);
              zaxpy(n-l,acc2,work,lp1,1,x,lp1,j,1);
              Clear(acc); Clear(acc2);
            end loop;
            for i in lp1..n loop
              Clear(work(i));
            end loop;
          end if;
          Clear(wrk);
          if wantv then
            -- place the transformation in v 
            -- for subsequent back multiplication
            for i in lp1..p loop
              Copy(e(i),v(i,l));
            end loop;
          end if;
        end if;
      end loop;
    end if;
   -- put_line("set up the final bidiagonal matrix of order m");
    m := min0(p,n+1);
    nctp1 := nct+1;
    nrtp1 := nrt+1;
    if nct < p
     then Copy(x(nctp1,nctp1),s(nctp1));
    end if;
    if n < m then
      Clear(s(m));
      s(m) := Create(integer(0));
    end if;
    if nrtp1 < m
     then Copy(x(nrtp1,m),e(nrtp1));
    end if;
    Clear(e(m));
    e(m) := Create(integer(0));
   -- put_line("if required, generate u");
    if wantu then
      if ncu >= nctp1 then
         for j in nctp1..ncu loop
           for i in 1..n loop
             Clear(u(i,j));
             u(i,j) := Create(integer(0));
           end loop;
           Clear(u(j,j));
           u(j,j) := Create(integer(1));
         end loop;
      end if;
      if nct >= 1 then
        for l in 1..nct loop
          ll := nct - l + 1;
          wrk := cabs1(s(ll));
          if Equal(wrk,0.0) then
            for i in 1..n loop
              Clear(u(i,ll));
              u(i,ll) := Create(integer(0));
            end loop;
            Clear(u(ll,ll));
            u(ll,ll) := Create(integer(1));
          else
            lp1 := ll + 1;
            if ncu >= lp1 then
              for j in lp1..ncu loop
                t := zdotc(n-ll+1,u,ll,ll,1,u,ll,j,1);
                Min(t);
                -- put("u("); put(ll,1); put(","); put(ll,1);
                -- put(") = "); put(u(ll,ll)); new_line;
                Div(t,u(ll,ll));
                -- put("t = "); put(t); new_line;
                zaxpy(n-ll+1,t,u,ll,ll,1,u,ll,j,1);
                Clear(t);
              end loop;
            end if;
            zscal(n-ll+1,c_min_one,u,ll,ll,1);
            Add(u(ll,ll),one);
            lm1 := ll - 1;
            if lm1 >= 1 then
              for i in 1..lm1 loop
                Clear(u(i,ll));
                u(i,ll) := Create(integer(0));
              end loop;
            end if;
          end if;
          Clear(wrk);
        end loop;
      end if;
    end if;
   -- put_line("The matrix u : "); put(u);
   -- put_line("if required, generate v");
    if wantv then
      for l in 1..p loop
        ll := p - l + 1;
        lp1 := ll + 1;
        if ll <= nrt then
          wrk := cabs1(e(ll));
          if wrk > 0.0 then -- if not Equal(wrk,0.0)
            for j in lp1..p loop
              t := zdotc(p-ll,v,lp1,ll,1,v,lp1,j,1);
              Min(t);
              Div(t,v(lp1,ll));
              zaxpy(p-ll,t,v,lp1,ll,1,v,lp1,j,1);
              Clear(t);
            end loop;
          end if;
          Clear(wrk);
        end if;
        for i in 1..p loop
          Clear(v(i,ll));
          v(i,ll) := Create(integer(0));
        end loop;
        Clear(v(ll,ll));
        v(ll,ll) := Create(integer(1));
      end loop;
    end if;
   -- put_line("The matrix v : "); put(v);
   -- put_line("transform s and e so that they are double precision");
    for i in 1..m loop
      wrk := cabs1(s(i));
      if wrk > 0.0 then -- if not Equal(wrk,0.0)
        tmp := cdabs(s(i));
        Clear(t); t := Create(tmp); Clear(tmp);
        Clear(r); r := s(i)/t;
        Copy(t,s(i));
        if i < m
         then Div(e(i),r);
        end if;
        if wantu then zscal(n,r,u,1,i,1); end if;
        Clear(r); Clear(t);
      end if;
      Clear(wrk);
      exit when (i = m);
      wrk := cabs1(e(i));
      if wrk > 0.0 then -- if not Equal(wrk,0.0)
        tmp := cdabs(e(i));
        Clear(t); t := Create(tmp); Clear(tmp);
        Clear(r); r := t/e(i);
        Copy(t,e(i));
        Mul(s(i+1),r);
        if wantv then zscal(p,r,v,1,i+1,1); end if;
        Clear(r); Clear(t);
      end if;
      Clear(wrk);
    end loop;
   -- put_line(" main iteration loop for the singular values");
    mm := m;
    iter := 0;
    loop
      exit when m = 0; -- quit if all the singular values have been found
      if iter > maxit  -- too many iterations have been performed
       then info := m; -- set the flag
            exit;      -- return
      end if;
     -- This section of the program inspects for negligible elements in 
     -- the s and e arrays.  On completion the variables kase and l are
     -- set as follows:
     --   kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
     --   kase = 2     if s(l) is negligible and l.lt.m
     --   kase = 3     if e(l-1) is negligible, l.lt.m, and
     --                  s(l), ..., s(m) are not negligible (qr step).
     --   kase = 4     if e(m-1) is negligible (convergence).
      for l in 1..m loop
        ll := m - l;
        exit when (ll = 0);
       -- test := cdabs(s(ll)) + cdabs(s(ll+1));
        Clear(test);
        test := cdabs(s(ll));
        wrk := cdabs(s(ll+1));
        Add(test,wrk);
        Clear(wrk);
       -- ztest := test + cdabs(e(ll));
        wrk := cdabs(e(ll));
        Clear(ztest);
        ztest := test + wrk;
        Clear(wrk);
        if Equal(ztest,test) then
          Clear(e(ll));
	  e(ll) := Create(integer(0));
          exit;
        end if;
      end loop;
      if ll = m - 1 then
        kase := 4;
      else 
        lp1 := ll + 1;
        mp1 := m + 1;
        for lls in lp1..mp1 loop
          ls := m - lls + lp1;
          exit when (ls = ll);
          Clear(test); test := Create(0.0);
          if ls /= n then
            wrk := cdabs(e(ls));
            Add(test,wrk);
            Clear(wrk);
          end if;
          if ls /= ll+1 then
            wrk := cdabs(e(ls-1));
            Add(test,wrk);
            Clear(wrk);
          end if;
          wrk := cdabs(s(ls));
          Clear(ztest);
          ztest := test + wrk;
          Clear(wrk);
          if Equal(ztest,test) then
            Clear(s(ls));
            s(ls) := Create(integer(0));
            exit;
          end if;
        end loop;
        if ls = ll then
          kase := 3;
        elsif ls = m then
          kase := 1;
        else
          kase := 2;
          ll := ls;
        end if;
      end if; 
      ll := ll + 1;
     -- perform the task indicated by kase
      case kase is
        when 1 => -- deflate negligible s(m)
                  mm1 := m-1;
                  Clear(f); f := REAL_PART(e(m-1));
                  Clear(e(m-1)); e(m-1) := Create(integer(0));
                  for k in ll..mm1 loop
                    kk := mm1 - k + ll;
                    Clear(t1); t1 := REAL_PART(s(kk));
                    drotg(t1,f,cs,sn);
                    Clear(s(kk)); s(kk) := Create(t1);
                    if kk /= ll
                     then -- f := -sn*REAL_PART(e(kk-1));
                          f := REAL_PART(e(kk-1));
                          Mul(f,sn);
                          Min(f);
                          Mul(e(kk-1),cs);
                    end if;
                    if wantv
                     then zdrot(p,v,1,kk,1,v,1,m,1,cs,sn);
                    end if;
                  end loop;
        when 2 => -- split at negligible s(ll)
                  Clear(f); f := REAL_PART(e(ll-1));
                  Clear(e(ll-1)); e(ll-1) := Create(integer(0));
                  for k in ll..m loop
                    Clear(t1); t1 := REAL_PART(s(k));
                    drotg(t1,f,cs,sn);
                    Clear(s(k)); s(k) := Create(t1);
                   -- f := -sn*REAL_PART(e(k));
                    Clear(f); f := REAL_PART(e(k));
                    Mul(f,sn);
                    Min(f);
                    Mul(e(k),cs);
                    if wantu
                     then zdrot(n,u,1,k,1,u,1,ll-1,1,cs,sn);
                    end if;
                  end loop;
        when 3 => -- perform one qr step
                  -- 1) calculate the shift
                  scale := dmax1(cdabs(s(m)),cdabs(s(m-1)),cdabs(e(m-1)),
                                 cdabs(s(ll)),cdabs(e(ll)));
                 -- sm := REAL_PART(s(m))/scale;
                  Clear(sm);
                  sm := REAL_PART(s(m));
                  Div(sm,scale);
                 -- smm1 := REAL_PART(s(m-1))/scale;
                  Clear(smm1);
                  smm1 := REAL_PART(s(m-1));
                  Div(smm1,scale);
                 -- emm1 := REAL_PART(e(m-1))/scale;
                  Clear(emm1);
                  emm1 := REAL_PART(e(m-1));
                  Div(emm1,scale);
                 -- sl := REAL_PART(s(ll))/scale;
                  Clear(sl);
                  sl := REAL_PART(s(ll));
                  Div(sl,scale);
                 -- el := REAL_PART(e(ll))/scale;
                  Clear(el);
                  el := REAL_PART(e(ll));
                  Div(el,scale);
                 -- b := ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0;
                  Copy(smm1,b);
                  Add(b,sm);
                  wrk := smm1 - sm;
                  Mul(b,wrk);
                  Copy(emm1,wrk);
                  Mul(wrk,wrk);
                  Add(b,wrk);
                  Div(b,2.0);
                 -- c := (sm*emm1)**2;
                  Copy(sm,c);
                  Mul(c,emm1);
                  Mul(c,c);
                  Clear(shift);
                  shift := Create(0.0);
                  if Equal(b,0.0) or Equal(c,0.0) then
                    -- shift := SQRT(b**2+c);
                    Copy(b,wrk);
                    Mul(wrk,b);
                    Add(wrk,c);
                    Clear(shift);
                    shift := SQRT(wrk);
                    Clear(wrk);
                    if b < 0.0 
                     then Min(shift);
                    end if;
                    -- shift := c/(b + shift);
                    Copy(b,wrk);
                    Add(wrk,shift);
                    Copy(c,shift);
                    Div(shift,wrk);
                    Clear(wrk);
                  end if;
                 -- f := (sl + sm)*(sl - sm) + shift;
                  Copy(sl,f);
                  Add(f,sm);
                  wrk := sl - sm;
                  Mul(f,wrk); Clear(wrk);
                  Add(f,shift);
                 -- g := sl*el;
                  Copy(sl,g);
                  Mul(g,el);
                  -- 2) chase zeros
                  mm1 := m - 1;
                  for k in ll..mm1 loop
                    drotg(f,g,cs,sn);
                    if k /= ll
                     then Clear(e(k-1)); e(k-1) := Create(f);
                    end if;
                   -- f := cs*REAL_PART(s(k)) + sn*REAL_PART(e(k));
                    Clear(f);
                    f := REAL_PART(s(k));
                    Mul(f,cs);
                    wrk := REAL_PART(e(k));
                    Mul(wrk,sn);
                    Add(f,wrk); Clear(wrk);
                   -- e(k) := cs*e(k) - sn*s(k);
                    Mul(e(k),cs);
                    acc := sn*s(k);
                    Sub(e(k),acc); Clear(acc);
                   -- g := sn*REAL_PART(s(k+1));
                    Clear(g);
                    g := REAL_PART(s(k+1));
                    Mul(g,sn);
                    Mul(s(k+1),cs);
                    if wantv
                     then zdrot(p,v,1,k,1,v,1,k+1,1,cs,sn);
                    end if;
                    drotg(f,g,cs,sn);
                    Clear(s(k)); s(k) := Create(f);
                   -- f := cs*REAL_PART(e(k)) + sn*REAL_PART(s(k+1));
                    Clear(f);
                    f := REAL_PART(e(k));
                    Mul(f,cs);
                    wrk := REAL_PART(s(k+1));
                    Mul(wrk,sn);
                    Add(f,wrk); Clear(wrk);
                   -- s(k+1) := -sn*e(k) + cs*s(k+1);
                    Mul(s(k+1),cs);
                    acc := sn*e(k);
                    Sub(s(k+1),acc); Clear(acc);
                   -- g := sn*REAL_PART(e(k+1));
                    Clear(g);
                    g := REAL_PART(e(k+1));   
                    Mul(g,sn);
                    Mul(e(k+1),cs);
                    if wantu and k < n
                     then zdrot(n,u,1,k,1,u,1,k+1,1,cs,sn);
                    end if;
                  end loop;
                  Clear(e(m-1)); e(m-1) := Create(f);
                  iter := iter + 1;
        when 4 => -- convergence
                  -- 1) make the singular value positive
                  wrk := REAL_PART(s(ll));
                  if wrk < 0.0
                   then Min(s(ll));
                        if wantv 
                         then zscal(p,c_min_one,v,1,ll,1);
                        end if;
                  end if;
                  -- 2) order the singular values
                  while ll /= mm loop
                    wrk := REAL_PART(s(ll));
                    tmp := REAL_PART(s(ll+1));
                    if wrk < tmp then
                      Copy(s(ll),t);
                      Copy(s(ll+1),s(ll));
                      Copy(t,s(ll+1)); Clear(t);
                      if wantv and ll < p
                       then zswap(p,v,1,ll,1,v,1,ll+1,1);
                      end if;
                      if wantu and ll < n
                       then zswap(n,u,1,ll,1,u,1,ll+1,1);
                      end if;
                      ll := ll+ 1;
                      leave_loop := false;
                    else
                      leave_loop := true;
                    end if;
                    Clear(wrk); Clear(tmp);
                    exit when leave_loop;
                  end loop;
                  iter := 0;
                  m := m-1;
        when others => null;
      end case;
    end loop;
    Clear(one); Clear(c_one);
    Clear(f_min_one); Clear(c_min_one);
 -- exception 
 --   when others => put_line("Exception raised in SVD.");
 --                 -- put_line("The current x matrix is "); put(x);
 --                  raise;
  end SVD;

  function Rank ( s : Vector ) return integer32 is

    val : Floating_Number;

  begin
    for i in s'range loop 
      val := AbsVal(s(i));
      if Equal(val,0.0) then
        Clear(val);
        return i-1;
      end if;
      Clear(val);
    end loop;
    return s'length;
  end Rank;

  function Rank ( s : Vector; tol : double_float ) return integer32 is

    val : Floating_Number;

  begin
    for i in s'range loop 
      val := AbsVal(s(i));
      if val < tol then
        Clear(val);
        return i-1;
      end if;
      Clear(val);
    end loop;
    return s'length;
  end Rank;

  function Inverse_Condition_Number
             ( s : Vector ) return double_float is

    res : double_float;
    mres : Floating_Number;
    smax : Floating_Number := AbsVal(s(s'first));
    smin,val : Floating_Number;

  begin
    if Equal(smax,0.0) then
      Clear(smax);
      res := 0.0;
    else
      Copy(smax,smin);
      for i in s'first+1..s'last loop
        val := AbsVal(s(i));
        exit when Equal(val,0.0);
        Copy(val,smin); Clear(val);
      end loop;
      mres := smin/smax;
      res := Round(mres);
      Clear(smin); Clear(smax); Clear(mres);       
    end if;
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( s : Vector; tol : double_float ) return double_float is

    res : double_float;
    smax : Floating_Number := AbsVal(s(s'first));
    smin,val : Floating_Number;

  begin
    if smax < tol then
      res := 0.0;
    else
      Copy(smax,smin);
      for i in s'first+1..s'last loop
        val := AbsVal(s(i));
        exit when (val < tol);
        Copy(val,smin);
      end loop;
      Div(smin,smax);
      res := Round(smin);
    end if;
    Clear(smin); Clear(smax); Clear(val);
    return res;
  end Inverse_Condition_Number;

  function Conjugate_Transpose ( z : Matrix ) return Matrix is

    res : Matrix(z'range(2),z'range(1));

  begin
    for i in z'range(1) loop
      for j in z'range(2) loop
        res(j,i) := Conjugate(z(i,j));
      end loop;
    end loop;
    return res;
  end Conjugate_Transpose;

  function Inverse ( u,v : Matrix; s : Vector ) return Matrix is

    ut : Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    su : Matrix(u'range(2),u'range(1));
    res : Matrix(v'range(1),su'range(2));
    val : Floating_Number;

  begin
    for i in s'range loop
      val := AbsVal(s(i));
      if Equal(val,0.0)
       then Clear(val); exit;
      end if;
      for j in ut'range(2) loop
        su(i,j) := ut(i,j)/s(i);
      end loop;
      Clear(val);
    end loop;
    res := v*su;
    Clear(ut); Clear(su);
    return res;
  end Inverse;

  function Inverse ( u,v : Matrix; s : Vector; tol : double_float )
                   return Matrix is

    ut : Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    su : Matrix(u'range(2),u'range(1));
    res : Matrix(v'range(1),su'range(2));
    val : Floating_Number;

  begin
    for i in s'range loop
      val := AbsVal(s(i));
      if val < tol
       then Clear(val); exit;
      end if;
      for j in ut'range(2) loop
        su(i,j) := ut(i,j)/s(i);
      end loop;
      Clear(val);
    end loop;
    res := v*su;
    Clear(ut); Clear(su);
    return res;
  end Inverse;

  function Solve ( u,v : Matrix; s,b : Vector ) return Vector is

    res : Vector(v'range(1));
    ut : Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    utb : Vector(u'range(2)) := ut*b;
    sub : Vector(v'range(1)) := (v'range(1) => Create(integer(0)));
    val : Floating_Number;

  begin
    for i in s'range loop
      val := AbsVal(s(i));
      if Equal(val,0.0)
       then Clear(val); exit;
      end if;
      exit when ((i > sub'last) or (i > utb'last));
      sub(i) := utb(i)/s(i);
      Clear(val);
    end loop;
    res := v*sub;
    Clear(ut); Clear(utb); Clear(sub);
    return res;
  end Solve;

  function Solve ( u,v : Matrix; s,b : Vector; tol : double_float )
                 return Vector is

    res : Vector(v'range(1));
    ut : Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    utb : Vector(u'range(2)) := ut*b;
    sub : Vector(v'range(1)) := (v'range(1) => Create(integer(0)));
    val : Floating_Number;

  begin
    for i in s'range loop
      val := AbsVal(s(i));
      if val < tol
       then Clear(val); exit;
      end if;
      sub(i) := utb(i)/s(i);
      Clear(val);
    end loop;
    Clear(ut); Clear(utb);
    res := v*sub;
    return res;
  end Solve;

end Multprec_Complex_Singular_Values;
