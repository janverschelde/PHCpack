-- with text_io,integer_io;                use text_io,integer_io;
-- with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
-- with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
-- with Penta_Double_Numbers_io;           use Penta_Double_Numbers_io;
-- with PentDobl_Complex_Numbers_io;       use PentDobl_Complex_Numbers_io;
-- with PentDobl_Complex_Vectors_io;       use PentDobl_Complex_Vectors_io;
-- with PentDobl_Complex_Matrices_io;      use PentDobl_Complex_Matrices_io;

with PentDobl_Complex_Numbers;          use PentDobl_Complex_Numbers;
with PentDobl_Mathematical_Functions;   use PentDobl_Mathematical_Functions;

package body PentDobl_Complex_Singular_Values is

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

  function dmax1 ( x1,x2 : penta_double ) return penta_double is

  -- DESCRIPTION :
  --   Returns the maximum of x1 and x2.

  begin
    if x1 > x2
     then return x1;
     else return x2;
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3 : penta_double ) return penta_double is

  -- DESCRIPTION :
  --   Returns the maximum of the three floating-point numbers.

  begin
    if x1 > x2
     then return dmax1(x1,x3);
     else return dmax1(x2,x3);
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3,x4 : penta_double ) return penta_double is

  -- DESCRIPTION :
  --   Returns the maximum of the four floating-point numbers.

  begin
    if x1 > x2
     then return dmax1(x1,x3,x4);
     else return dmax1(x2,x3,x4);
    end if;
  end dmax1;

  function dmax1 ( x1,x2,x3,x4,x5 : penta_double ) return penta_double is

  -- DESCRIPTION :
  --   Returns the maximum of the five floating-point numbers.

  begin
    if x1 > x2
     then return dmax1(x1,x3,x4,x5);
     else return dmax1(x2,x3,x4,x5);
    end if;
  end dmax1;

  function cabs1 ( z : Complex_Number ) return penta_double is

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of real and imaginary
  --   part of the complex number z.  Translation of
  --     complex*16 zdum
  --     double precision cabs1
  --     cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))

  begin
    return (abs(REAL_PART(z)) + abs(IMAG_PART(z)));
  end cabs1;

  function cdabs ( z : Complex_Number ) return penta_double is

  -- DESCRIPTION :
  --   Computes the modulus of the complex number, let us hope this
  --   corresponds to the `cdabs' fortran function.

  -- NOTE : 
  --   The SQRT crashes for very tiny numbers because of overflow
  --   in Newton's method, which is applied without scaling.
  --   The patch is to check for zero imaginary part.
  --   For efficiency, this should have been done anyway ...

    rpz : constant penta_double := REAL_PART(z);
    ipz : constant penta_double := IMAG_PART(z);

  begin
    if is_zero(ipz)
     then return AbsVal(rpz);
     else return SQRT(sqr(rpz) + sqr(ipz));
    end if;
 -- exception
 --   when others =>
 --     put_line("exception raised in cdabs ...");
 --     put("z = "); put(z); new_line;
 --     raise;
  end cdabs;

  function csign ( z1,z2 : Complex_Number ) return Complex_Number is

  -- DESCRIPTION : translated from
  --       csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2)) 

  begin
    return (Create(cdabs(z1)/cdabs(z2))*z2);
 -- exception
 --   when others => put_line("exception caught by csign");
 --                  put("z1 = "); put(z1); new_line;
 --                  put("z2 = "); put(z2); new_line;
 --                 -- put("cdabs(z2) = "); put(cdabs(z2)); new_line;
 --                  raise;
  end csign;

  function dsign ( a,b : penta_double ) return penta_double is

  -- DESCRIPTION :
  --   The implementation of this routine is written from web page
  --   documentation of sign...

  begin
    if b >= 0.0
     then return abs(a);
     else return -abs(a);
    end if;
  end dsign;

  function dznrm2 ( n : integer32; x : Vector; ind,incx : integer32 )
                  return penta_double is

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector x, starting at x(ind)
  --   and continueing n steps with increment incx.

    ix : integer32;
    zero : constant penta_double := create(0.0);
    one : constant penta_double := create(1.0);
    norm,scale,ssq,temp : penta_double;

  begin
    if n < 1 or incx < 1 then
      norm := zero;
    else
      scale := zero;
      ssq := one;
      ix := ind;
      while ix <= ind + (n-1)*incx loop
        if REAL_PART(x(ix)) /= zero then
          temp := abs(REAL_PART(x(ix)));
          if scale < temp then
            ssq := 1.0 + ssq*sqr(scale/temp);
            scale := temp;
          else
            ssq := ssq + sqr(temp/scale);
          end if;
        end if;
        if IMAG_PART(x(ix)) /= zero then
          temp := abs(IMAG_PART(x(ix)));
          if scale < temp then
            ssq := 1.0 + ssq*sqr(scale/temp);
            scale := temp;
          else
            ssq := ssq + sqr(temp/scale);
          end if;
        end if;
        ix := ix + incx;
      end loop;
      norm := scale*SQRT(ssq);
    end if;
    return norm;
 -- exception
 --   when others => put_line("exception caught by 1st dznrm2"); raise;
  end dznrm2;

  function dznrm2 ( n : integer32; x : Matrix; row,col,incx : integer32 )
                  return penta_double is

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector x, starting at x(row,col)
  --   and continueing n steps with increment incx in the same column.

    ix : integer32;
    zero : constant penta_double := create(0.0);
    one : constant penta_double := create(1.0);
    norm,scale,ssq,temp : penta_double;

  begin
    if n < 1 or incx < 1 then
      norm := zero;
    else
      scale := zero; ssq := one; ix := row;
      while ix <= row + (n-1)*incx loop
        if REAL_PART(x(ix,col)) /= zero then
          temp := abs(REAL_PART(x(ix,col)));
          if scale < temp then
            ssq := 1.0 + ssq*sqr(scale/temp); scale := temp;
          else
            ssq := ssq + sqr(temp/scale);
          end if;
        end if;
        if IMAG_PART(x(ix,col)) /= zero then
          temp := abs(IMAG_PART(x(ix,col)));
          if scale < temp then
            ssq := 1.0 + ssq*sqr(scale/temp); scale := temp;
          else
            ssq := ssq + sqr(temp/scale);
          end if;
        end if;
        ix := ix + incx;
      end loop;
      norm := scale*SQRT(ssq);
    end if;
    return norm;
 -- exception
 --   when others => put_line("exception caught by 2nd dznrm2"); raise;
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
        zx(ind+i) := za*zx(ind+i);
      end loop;
    else
      ix := ind; 
      for i in 1..n loop
        zx(ix) := za*zx(ix);
        ix := ix + incx;
      end loop;
    end if;
 -- exception
 --   when others => put_line("exception caught by 1st zscal"); raise;
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
        zx(row+i,col) := za*zx(row+i,col);
      end loop;
    else
      ix := row;
      for i in 1..n loop
        zx(ix,col) := za*zx(ix,col);
        ix := ix + incx;
      end loop;
    end if;
 -- exception
 --   when others => put_line("exception caught by 2nd zscal"); raise;
  end zscal;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Vector; ind,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at ind and in y
  --   at (rwy,cly), using increments incx and incy to advance.

    ix,iy : integer32;
    zero : constant penta_double := create(0.0);

  begin
    if n > 0 and then AbsVal(z) /= zero then
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
 -- exception
 --   when others => put_line("exception caught by 1st zaxpy"); raise;
  end zaxpy;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Vector; ind,incy : in integer32 ) is

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at (rwx,clx)
  --   and in y at ind, using increments incx and incy to advance.

    zero : constant penta_double := create(0.0);
    ix,iy : integer32;

  begin
    if n > 0 and then AbsVal(z) /= zero then
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
 -- exception
 --   when others => put_line("exception caught by 2nd zaxpy"); raise;
  end zaxpy;

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 ) is

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x and y at the
  --   respective (rows,columns): (rwx,clx) and (rwy,cly), using
  --   increments incx and incy to advance in the rows.

    zero : constant penta_double := create(0.0);
    ix,iy : integer32;

  begin
    if n > 0 and then AbsVal(z) /= zero then
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
 -- exception
 --   when others => put_line("exception caught by 3rd zaxpy"); raise;
  end zaxpy;

  function zdotc ( n : in integer32; x : in Matrix;
                   rwx,clx,incx : in integer32;
                   y : in Matrix; rwy,cly,incy : in integer32 )
                 return Complex_Number is

  -- DESCRIPTION :
  --   Returns the dot product of two vectors in two columns of matrices
  --   x and y, starting at rows rwx and rwy respectively, with respective
  --   increments in incx and incy.

    zero : constant penta_double := create(0.0);
    ztemp : Complex_Number := Create(zero);
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
 -- exception
 --   when others => put_line("exception caught by zdotc"); raise;
  end zdotc;

  procedure drotg ( da,db,c,s : in out penta_double ) is

  -- DESCRIPTION :
  --   Constructs Givens plane rotation.

    zero : constant penta_double := create(0.0);
    one : constant penta_double := create(1.0);
    roe,scale,r,z : penta_double;

  begin
   -- put("in drotg da = "); put(da); new_line;
   -- put("         db = "); put(db); new_line;
    roe := db;
    if abs(da) > abs(db) then roe := da; end if;
    scale := abs(da) + abs(db);
    if one + scale = one then
      c := one;  s := zero;
      r := zero; z := zero;
    else
      r := scale*SQRT(sqr(da/scale) + sqr(db/scale));
      r := dsign(one,roe)*r;
      c := da/r;
      s := db/r;
      z := one;
      if abs(da) > abs(db) then z := s; end if;
      if abs(db) >= abs(da) and c /= zero
       then z := 1.0/c;
      end if;
    end if;
    da := r;
    db := z;
   -- put("out drotg c = "); put(c); new_line;
   -- put("          s = "); put(s); new_line;
   -- put("         da = "); put(da); new_line;
   -- put("         db = "); put(db); new_line;
 -- exception
 --   when others => put_line("exception caught by drotg"); raise;
  end drotg;

  procedure zdrot ( n : in integer32;
                    x : in out Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32;
                    c,s : in penta_double ) is

  -- DESCRIPTION :
  --   Applies a plane rotation where the cos and sin are c and s
  --   and the vectors are in the columns of the matrices x and y,
  --   starting at (rwx,clx) and (rwy,cly) advancing in the rows with
  --   increments incx and incy respectively.

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
 -- exception
 --   when others => put_line("exception caught by zdrot"); raise;
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
 -- exception
 --   when others => put_line("exception caught by zswap"); raise;
  end zswap;

-- TARGET PROCEDURES :

  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix; 
                  job : in integer32; info : out integer32 ) is

    work : Vector(1..n);

  begin
    SVD(x,n,p,s,e,u,v,job,info,work);
  end SVD;

  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix; 
                  job : in integer32; info : out integer32;
                  work : in out Vector ) is

    iter,jobu,kase,kk,ll,lm1,lp1,ls,lu,m,maxit : integer32;
    mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1 : integer32;
    t,r : Complex_Number;
    b,c,cs,el,emm1,f,g,scale,shift,sl,sm,sn : penta_double;
    smm1,t1,test,ztest : penta_double;
    wantu,wantv : boolean;
    zero : constant penta_double := create(0.0);
    one : constant penta_double := create(1.0);
    minusone : constant penta_double := create(-1.0);

    cntkase2 : integer32 := 0; -- patch to fix with NaN inputs...
 
  begin
   -- put("entering SVD ...");
   -- put_line("the input matrix x : "); put(x);
    maxit := 30;     -- maximal number of iterations
    wantu := false;  -- determine what is to be computed
    wantv := false;
    jobu := (job mod 100)/10;
    ncu := n;
    if jobu > 1 then ncu := min0(n,p); end if;
    if jobu /= 0 then wantu := true; end if;
    if job mod 10 /= 0 then wantv := true; end if;
   -- reduce x to bidiagonal form, storing the diagonal elements in s
   -- and the super diagonal elements in e
    info := 0;
    nct := min0(n-1,p);
    nrt := max0(0,min0(p-2,n));
    lu := max0(nct,nrt);
    if lu >= 1 then
      for l in 1..lu loop
        lp1 := l+1;
        if l <= nct then
          -- compute the transformation for the l-th column
          -- and place the l-th diagonal in s(l)
          s(l) := Create(dznrm2(n-l+1,x,l,l,1),zero);
          -- put("s("); put(l,1); put(") = "); put(s(l)); new_line;
          if cabs1(s(l)) /= zero then
            if cdabs(x(l,l)) /= zero
             then s(l) := csign(s(l),x(l,l));
            end if;
            -- put("s("); put(l,1); put(") = "); put(s(l)); new_line;
            -- put("inverse : "); put(Create(one)/s(l)); new_line;
            zscal(n-l+1,Create(one)/s(l),x,l,l,1);
            x(l,l) := Create(one) + x(l,l);
            -- put("x("); put(l,1); put(","); put(l,1);
            -- put(") = "); put(x(l,l)); new_line;
          end if;
          s(l) := -s(l);
        end if;
        -- put_line("The matrix x : "); put(x);
        if p >= lp1 then
          for j in lp1..p loop
            if l <= nct then
              if (cabs1(s(l)) /= zero)
               then -- apply the transformation
                    t := -zdotc(n-l+1,x,l,l,1,x,l,j,1)/x(l,l);
                    zaxpy(n-l+1,t,x,l,l,1,x,l,j,1);
              end if;
            end if;
            -- place the l-th row of x into e for the subsequent
            -- calculation of the row transformation
            e(j) := Conjugate(x(l,j));
          end loop;
        end if;
        -- put_line("The matrix x : "); put(x);
        if wantu and l <= nct then
          -- place the transformation in u for subsequent
          -- back multiplication
          for i in l..n loop
            u(i,l) := x(i,l);
          end loop;
        end if;
        if l <= nrt then
           -- compute the l-th row transformation
           -- and place the l-th super diagonal in e(l)
          e(l) := Create(dznrm2(p-l,e,lp1,1),zero);
          if cabs1(e(l)) /= zero then
            if cdabs(e(lp1)) /= zero
             then e(l) := csign(e(l),e(lp1));
            end if;
            zscal(p-l,Create(one)/e(l),e,lp1,1);
            e(lp1) := Create(one) + e(lp1);
          end if;
          e(l) := -Conjugate(e(l));
          if lp1 <= n and cabs1(e(l)) /= zero then
             -- apply the transformation
            for i in lp1..n loop
              work(i) := Create(zero);
            end loop;
            for j in lp1..p loop
              zaxpy(n-l,e(j),x,lp1,j,1,work,lp1,1);
            end loop;
            for j in lp1..p loop
              zaxpy(n-l,Conjugate(-e(j)/e(lp1)),work,lp1,1,x,lp1,j,1);
            end loop;
          end if;
          if wantv then
            -- place the transformation in v
            -- for subsequent back multiplication
            for i in lp1..p loop
              v(i,l) := e(i);
            end loop;
          end if;
        end if;
      end loop;
    end if;
   -- set up the final bidiagonal matrix of order m
    m := min0(p,n+1);
    nctp1 := nct+1;
    nrtp1 := nrt+1;
    if nct < p then s(nctp1) := x(nctp1,nctp1); end if;
    if n < m then s(m) := Create(zero); end if;
    if nrtp1 < m then e(nrtp1) := x(nrtp1,m); end if;
    e(m) := Create(zero);
   -- if required, generate u
    if wantu then
      if ncu >= nctp1 then
        for j in nctp1..ncu loop
          for i in 1..n loop
            u(i,j) := Create(zero);
          end loop;
          u(j,j) := Create(one);
        end loop;
      end if;
      if nct >= 1 then
        for l in 1..nct loop
          ll := nct - l + 1;
          if cabs1(s(ll)) = zero then
            for i in 1..n loop
              u(i,ll) := Create(zero);
            end loop;
            u(ll,ll) := Create(one);
          else
            lp1 := ll + 1;
            if ncu >= lp1 then
              for j in lp1..ncu loop
                t := -zdotc(n-ll+1,u,ll,ll,1,u,ll,j,1)/u(ll,ll);
                zaxpy(n-ll+1,t,u,ll,ll,1,u,ll,j,1);
              end loop;
            end if;
            zscal(n-ll+1,Create(minusone),u,ll,ll,1);
            u(ll,ll) := Create(one) + u(ll,ll);
            lm1 := ll - 1;
            if lm1 >= 1 then
              for i in 1..lm1 loop
                u(i,ll) := Create(zero);
              end loop;
            end if;
          end if;
        end loop;
      end if;
    end if;
   -- put_line("The matrix u : "); put(u);
   -- if required, generate v
    if wantv then
      for l in 1..p loop
        ll := p - l + 1;
        lp1 := ll + 1;
        if ll <= nrt then
          if cabs1(e(ll)) /= zero then
            for j in lp1..p loop
              t := -zdotc(p-ll,v,lp1,ll,1,v,lp1,j,1)/v(lp1,ll);
              zaxpy(p-ll,t,v,lp1,ll,1,v,lp1,j,1);
            end loop;
          end if;
        end if;
        for i in 1..p loop
          v(i,ll) := Create(zero);
        end loop;
        v(ll,ll) := Create(one);
      end loop;
    end if;
   -- put_line("The matrix v = "); put(v);
   -- transform s and e so that they are double precision
    for i in 1..m loop
      if cabs1(s(i)) /= zero then
        t := Create(cdabs(s(i)),zero);
        r := s(i)/t;
        s(i) := t;
        if i < m then e(i) := e(i)/r; end if;
        if wantu then zscal(n,r,u,1,i,1); end if;
      end if;
      exit when (i = m);
      if cabs1(e(i)) /= zero then
        t := Create(cdabs(e(i)),zero);
        r := t/e(i);
        e(i) := t;
        s(i+1) := s(i+1)*r;
        if wantv then zscal(p,r,v,1,i+1,1); end if;
      end if;
    end loop;
   -- main iteration loop for the singular values
    mm := m;
    iter := 0;
    loop
     -- put("inside SVD loop with m = "); put(m,1);
     -- put(" and iter = "); put(iter,1); new_line;
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
        test := cdabs(s(ll)) + cdabs(s(ll+1));
        ztest := test + cdabs(e(ll));
        if ztest = test
         then e(ll) := Create(zero); exit;
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
          test := zero;
          if ls /= n then test := test + cdabs(e(ls)); end if;
          if ls /= ll+1 then test := test + cdabs(e(ls-1)); end if;
          ztest := test + cdabs(s(ls));
          if ztest = test 
           then s(ls) := Create(zero); exit;
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
     -- put(" kase = "); put(kase,1); new_line;
      case kase is
        when 1 => -- deflate negligible s(m)
          mm1 := m-1;
          f := REAL_PART(e(m-1));
          e(m-1) := Create(zero);
          for k in ll..mm1 loop
            kk := mm1 - k + ll;
            t1 := REAL_PART(s(kk));
            cs := zero; sn := zero;
            drotg(t1,f,cs,sn);
            s(kk) := Create(t1);
            if kk /= ll then
              f := -sn*REAL_PART(e(kk-1));
              e(kk-1) := cs*e(kk-1);
            end if;
            if wantv
             then zdrot(p,v,1,kk,1,v,1,m,1,cs,sn);
            end if;
          end loop;
        when 2 => -- split at negligible s(ll)
          f := REAL_PART(e(ll-1));
          e(ll-1) := Create(zero);
          for k in ll..m loop
            t1 := REAL_PART(s(k));
            drotg(t1,f,cs,sn);
            s(k) := Create(t1);
            f := -sn*REAL_PART(e(k));
            e(k) := cs*e(k);
            if wantu
             then zdrot(n,u,1,k,1,u,1,ll-1,1,cs,sn);
            end if;
          end loop;
          cntkase2 := cntkase2 + 1;  -- patch for fixing NaN inputs ...
          if cntkase2 > 100 then return; end if;
        when 3 => -- perform one qr step
                  -- 1) calculate the shift
           scale := dmax1(cdabs(s(m)),cdabs(s(m-1)),cdabs(e(m-1)),
                          cdabs(s(ll)),cdabs(e(ll)));
           sm := REAL_PART(s(m))/scale;
           smm1 := REAL_PART(s(m-1))/scale;
           emm1 := REAL_PART(e(m-1))/scale;
           sl := REAL_PART(s(ll))/scale;
           el := REAL_PART(e(ll))/scale;
           b := ((smm1 + sm)*(smm1 - sm) + sqr(emm1))/2.0;
           c := sqr(sm*emm1);
           shift := zero;
           if b = zero or c = zero then
             shift := SQRT(sqr(b)+c);
             if b < 0.0 then shift := -shift; end if;
             shift := c/(b + shift);
           end if;
           f := (sl + sm)*(sl - sm) + shift;
           g := sl*el;
          -- 2) chase zeros
           mm1 := m - 1;
          -- put("chasing zeros, ll = "); put(ll,1);
          -- put(", mm1 = "); put(mm1,1); new_line;
           for k in ll..mm1 loop
             drotg(f,g,cs,sn);
             if k /= ll then e(k-1) := Create(f); end if;
             f := cs*REAL_PART(s(k)) + sn*REAL_PART(e(k));
             e(k) := cs*e(k) - sn*s(k);
             g := sn*REAL_PART(s(k+1));
             s(k+1) := cs*s(k+1);
             if wantv
              then zdrot(p,v,1,k,1,v,1,k+1,1,cs,sn);
             end if;
             drotg(f,g,cs,sn);
             s(k) := Create(f);
             f := cs*REAL_PART(e(k)) + sn*REAL_PART(s(k+1));
             s(k+1) := -sn*e(k) + cs*s(k+1);
             g := sn*REAL_PART(e(k+1));
             e(k+1) := cs*e(k+1);
             if wantu and k < n
              then zdrot(n,u,1,k,1,u,1,k+1,1,cs,sn);
             end if;
           end loop;
           e(m-1) := Create(f);
           iter := iter + 1;
        when 4 => -- convergence
                  -- 1) make the singular value positive
          if REAL_PART(s(ll)) < 0.0 then
            s(ll) := -s(ll);
            if wantv 
             then zscal(p,Create(minusone),v,1,ll,1);
            end if;
          end if;
         -- 2) order the singular values
          while ll /= mm loop
            exit when REAL_PART(s(ll)) >= REAL_PART(s(ll+1));
            t := s(ll);
            s(ll) := s(ll+1);
            s(ll+1) := t;
            if wantv and ll < p
             then zswap(p,v,1,ll,1,v,1,ll+1,1);
            end if;
            if wantu and ll < n
             then zswap(n,u,1,ll,1,u,1,ll+1,1);
            end if;
            ll := ll+ 1;
          end loop;
          iter := 0;
          m := m-1;
        when others => null;
      end case;
    end loop;
   -- put_line("... leaving SVD");
 -- exception
 --   when others 
 --     => put_line("exception caught by SVD"); raise;
  end SVD;

  function Rank ( s : Vector ) return integer32 is

    one : constant penta_double := create(1.0);

  begin
    for i in s'range loop 
      if AbsVal(s(i)) + one = one
       then return i-1;
      end if;
    end loop;
    return s'length;
  end Rank;

  function Rank ( s : Vector; tol : double_float ) return integer32 is
  begin
    for i in s'range loop 
      if AbsVal(s(i)) < tol
       then return i-1;
      end if;
    end loop;
    return s'length;
  end Rank;

  function Inverse_Condition_Number
             ( s : Vector ) return penta_double is

    smax : constant penta_double := AbsVal(s(s'first));
    one : constant penta_double := create(1.0);
    smin,val : penta_double;

  begin
    if smax + one = one then
      return create(0.0);
    else
      smin := smax;
      for i in s'first+1..s'last loop
        val := AbsVal(s(i));
        exit when (val + one = one);
        smin := val;
      end loop;
      return smin/smax;
    end if;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( s : Vector; tol : double_float ) return penta_double is

    smax : constant penta_double := AbsVal(s(s'first));
    smin,val : penta_double;

  begin
    if smax < tol then
      return create(0.0);
    else
      smin := smax;
      for i in s'first+1..s'last loop
        val := AbsVal(s(i));
        exit when (val < tol);
        smin := val;
      end loop;
      return smin/smax;
    end if;
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

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    su : Matrix(u'range(2),u'range(1));
    one : constant penta_double := create(1.0);

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) + one = one);
      for j in ut'range(2) loop
        su(i,j) := ut(i,j)/s(i);
      end loop;
    end loop;
    return v*su;
  end Inverse;

  function Inverse ( u,v : Matrix; s : Vector; tol : double_float )
                   return Matrix is

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    su : Matrix(u'range(2),u'range(1));

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) < tol);
      for j in ut'range(2) loop
        su(i,j) := ut(i,j)/s(i);
      end loop;
    end loop;
    return v*su;
  end Inverse;

  function Solve ( u,v : Matrix; s,b : Vector ) return Vector is

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    utb : constant Vector(u'range(2)) := ut*b;
    zero : constant penta_double := create(0.0);
    one : constant penta_double := create(1.0);
    sub : Vector(v'range(1)) := (v'range(1) => Create(zero));

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) + one = one);
      exit when ((i > sub'last) or (i > utb'last));
      sub(i) := utb(i)/s(i);
    end loop;
    return v*sub;
  end Solve;

  procedure Solve ( ut,v : in Matrix; s,b : in Vector;
                    utb,sub : in out Vector; sol : out Vector ) is

    one : constant penta_double := create(integer(1));

  begin
    utb := ut*b;
    sub := (v'range(1) => Create(integer(0)));
    for i in s'range loop
      exit when (AbsVal(s(i)) + one = one);
      exit when ((i > sub'last) or (i > utb'last));
      sub(i) := utb(i)/s(i);
    end loop;
    sol := v*sub;
  end Solve;

  function Solve ( u,v : Matrix; s,b : Vector; tol : double_float )
                 return Vector is

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    utb : constant Vector(u'range(2)) := ut*b;
    zero : constant penta_double := create(0.0);
    sub : Vector(v'range(1)) := (v'range(1) => Create(zero));

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) < tol);
      sub(i) := utb(i)/s(i);
    end loop;
    return v*sub;
  end Solve;

end PentDobl_Complex_Singular_Values;
