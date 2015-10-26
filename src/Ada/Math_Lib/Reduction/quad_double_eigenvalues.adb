with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with QuadDobl_Mathematical_Functions;     use QuadDobl_Mathematical_Functions;

package body Quad_Double_Eigenvalues is

-- ACKNOWLEDGMENT :
--   The fortran code contained the following:
--    (1) "this subroutine is a translation of the algol procedure ..."
--    (2) "questions and comments should be directed to burton s. garbow,
--         mathematics and computer science div, argonne national laboratory"

-- TOP LEVEL WRAPPERS :

  procedure Eigenvalues ( A : in Matrix; ierr : out integer32;
                          L : out QuadDobl_Complex_Vectors.Vector ) is

    wr,wi : Quad_Double_Vectors.Vector(A'range(1));
    z : Matrix(A'range(1),A'range(2));

  begin
    RG(A'last(2),A'last(1),A,0,wr,wi,z,ierr);
    L := Create(wr,wi);
  end Eigenvalues;

  procedure Eigenvectors ( A : in Matrix; ierr : out integer32;
                           L : out QuadDobl_Complex_Vectors.Vector;
                           v : out QuadDobl_Complex_VecVecs.VecVec ) is

    wr,wi : Quad_Double_Vectors.Vector(A'range(1));
    z : Matrix(A'range(1),A'range(2));

  begin
    RG(A'last(2),A'last(1),A,1,wr,wi,z,ierr);
    L := Create(wr,wi);
    v := Create(z,wi);
  end Eigenvectors;

-- AUXILIARIES FOR WRAPPERS :

  function Create ( wr,wi : Quad_Double_Vectors.Vector )
                  return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(wr'range);

  begin
    for k in res'range loop
      res(k) := Create(wr(k),wi(k));
    end loop;
    return res;
  end Create;

  function Create ( z : in Quad_Double_Matrices.Matrix;
                    wi : in Quad_Double_Vectors.Vector ) 
                  return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(wi'range);
    ind : integer32 := res'first - 1;
    one : constant quad_double := create(1.0);

  begin
    while ind < res'last loop
      ind := ind + 1;
      res(ind) := new QuadDobl_Complex_Vectors.Vector(z'range(1));
      if wi(ind) + one = one then
        for i in z'range(1) loop
          res(ind)(i) := Create(z(i,ind));
        end loop;
      else
        res(ind+1) := new QuadDobl_Complex_Vectors.Vector(z'range(1));
        for i in z'range(1) loop
          res(ind)(i) := Create(z(i,ind),z(i,ind+1));
          res(ind+1)(i) := Create(z(i,ind),-z(i,ind+1));
        end loop;
        ind := ind + 1;
      end if;
    end loop;
    return res;
  end Create;

-- AUXILIARIES FOR TARGET ROUTINES :

  function min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end min0;

  function dsign ( a,b : quad_double ) return quad_double is

  -- DESCRIPTION :
  --   The implementation of this routine is written from web page
  --   documentation of sign...

  begin
    if b >= 0.0
     then return abs(a);
     else return -abs(a);
    end if;
  end dsign;

  function dmax1 ( a,b : quad_double ) return quad_double is

  -- DESCRIPTION : returns max(a,b), what else could it be?

  begin
    if a >= b 
     then return a;
     else return b;
    end if;
  end dmax1;

-- TARGET ROUTINES :

  procedure RG ( nm,n : in integer32; A : in Matrix; matz : in integer32;
                 wr,wi : out Quad_Double_Vectors.Vector;
                 z : out Matrix; ierr : out integer32 ) is

    wrk : Quad_Double_Matrices.Matrix(A'range(1),A'range(2));
    fv1 : Quad_Double_Vectors.Vector(1..n);
    iv1 : Standard_Integer_Vectors.Vector(1..n);
    is1,is2 : integer32;

  begin
    if n > nm then
      ierr := 10*n;
    else
      ierr := 0;
      wrk := A;
     -- put("calling balanc ...");
      balanc(nm,n,wrk,is1,is2,fv1);
     -- put_line("the matrix after balanc : "); put(wrk);
     -- put_line(" done");
     -- put("calling elmhes ...");
      elmhes(nm,n,is1,is2,wrk,iv1);
     -- put_line("the matrix after elmhes : "); put(wrk);
     -- put_line(" done");
      if matz = 0 then  -- find eigenvalues only
       -- put("calling hqr ...");
        hqr(nm,n,is1,is2,wrk,wr,wi,ierr);
       -- put_line(" done");
      else              -- find both eigenvalues and eigenvectors
        eltran(nm,n,is1,is2,wrk,iv1,z);
       -- put_line("the matrix after eltran : "); put(wrk);
       -- put_line("the matrix z after eltran : "); put(z);
        hqr2(nm,n,is1,is2,wrk,wr,wi,z,ierr);
       -- put_line("the matrix after hqr2 : "); put(wrk);
       -- put_line("the matrix z after hqr2 : "); put(z);
        if (ierr = 0) then
          balbak(nm,n,is1,is2,fv1,n,z);
         -- put_line("the matrix z after balbak : "); put(z);
        end if;
      end if;
    end if;
  end RG;

-- SUBROUTINES :

  procedure balanc ( nm,n : in integer32; A : in out Matrix;
                     low,igh : out integer32;
                     scale : out Quad_Double_Vectors.Vector ) is

  -- NOTE : the fortran code contained the following reference:
  --   This subroutine is a translation of the algol procedure balance,
  --   num. math. 13, 293-304(1969) by parlett and reinsch.
  --   handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

    noconv,nonzero,nextloop : boolean;
    j,k,L,m,iexc : integer32;
    c,f,g,r,s,b2,radix : quad_double;
    zero : constant quad_double := create(0.0);

    procedure Exchange is -- performs row and column exchange
    begin
      scale(m) := Quad_Double_Numbers.create(j);
      if (j /= m) then
        for ii in 1..L loop
          f := A(ii,j);
          A(ii,j) := A(ii,m);
          A(ii,m) := f;
        end loop;
        for ii in k..n loop
          f := A(j,ii);
          A(j,ii) := A(m,ii);
          A(m,ii) := f;
        end loop;
      end if;
      if iexc = 1 then
        if L = 1
         then low := k; igh := L; return;
        end if;  -- search for rows isolating an eigenvalue and push them down
        L := L - 1;
        nextloop := false;
      else
        k := k + 1;
        nextloop := true;
      end if;
     -- put_line(" done");
    end Exchange;

  begin
    radix := create(16.0);
    b2 := radix*radix;
    k := 1; L := n;
    nextloop := false;
    if not nextloop then
      for jj in 1..L loop    -- for j = L step -1 until 1 do
        j := L + 1 - jj;
        nonzero := false;
        for ii in 1..L loop
          if (ii /= j) then
            if (A(j,ii) /= zero) then
              nonzero := true; exit;
            end if;
          end if;
        end loop;
        if not nonzero then
          m := L;
          iexc := 1;
          Exchange;
        end if;
        exit when nextloop;
      end loop;
    end if;
    for j in k..L loop
      nonzero := false;
      for ii in k..L loop
        if (ii /= j) then
          if (A(ii,j) /= zero) then
            nonzero := true; exit;
          end if;
        end if;
      end loop;
      if not nonzero then
        m := k;
        iexc := 2;
        Exchange;
      end if;
    end loop;
  -- now balance the submatrix in rows k to L
    for ii in k..L loop
      scale(ii) := create(1.0);
    end loop;
    loop                         -- iterative loop for norm reduction 
      noconv := false;
      for ii in k..L loop
        c := zero; r := zero;
        for jj in k..L loop
          if jj /= ii then
            c := c + abs(A(jj,ii));
            r := r + abs(A(ii,jj));
          end if;
        end loop; -- guard against zero c or r due to underflow
        if ((c /= zero) and (r /= zero)) then
          g := r / radix;
          f := create(1.0);
          s := c + r;
          loop
            exit when (c >= g);
            f := f * radix;
            c := c * b2;
          end loop;
          g := r * radix;
          loop
            exit when (c < g);
            f := f / radix;
            c := c / b2;
          end loop;
          if ((c + r) / f < 0.95 * s) then -- now balance 
            g := 1.0/f;
            scale(ii) := scale(ii)*f;
            noconv := true;
            for jj in k..n loop
              A(ii,jj) := A(ii,jj)*g;
            end loop;
            for jj in 1..L loop
              A(jj,ii) := A(jj,ii)*f;
            end loop;
          end if;
        end if;
      end loop;
      exit when (not noconv);
    end loop;
    low := k; igh := L;
  end balanc;

  procedure elmhes ( nm,n,low,igh : in integer32; A : in out Matrix;
                     int : out Standard_Integer_Vectors.Vector) is

  -- NOTE : the fortran code contained the following reference:
  --   this subroutine is a translation of the algol procedure elmhes,
  --   num. math. 12, 349-368(1968) by martin and wilkinson.
  --   handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

    la,kp1,i,mm1,mp1 : integer32;
    x,y : quad_double;
    zero : constant quad_double := create(0.0);

  begin
    la := igh - 1;
    kp1 := low + 1;
    if la >= kp1 then
      for m in kp1..la loop
        mm1 := m - 1;
        x := zero;
        i := m;
        for j in m..igh loop
          if (abs(A(j,mm1)) > abs(x))
           then x := A(j,mm1); i := j;
          end if;
        end loop;
        int(m) := i;
        if (i /= m) then
          for j in mm1..n loop   -- interchange rows and columns of A
            y := A(i,j);
            A(i,j) := A(m,j);
            A(m,j) := y;
          end loop;
          for j in 1..igh loop
            y := A(j,i);
            A(j,i) := A(j,m);
            A(j,m) := y;
          end loop;              -- end interchange
        end if;
        if (x /= zero) then
          mp1 := m + 1;
          for ii in mp1..igh loop
            y := a(ii,mm1);
            if y /= zero then
              y := y / x;
              A(ii,mm1) := y;
              for j in m..n loop
                A(ii,j) := A(ii,j) - y*A(m,j);
              end loop;
              for j in 1..igh loop
                A(j,m) := A(j,m) + y*A(j,ii);
              end loop;
            end if;
          end loop;
        end if;
      end loop;
    end if;
  end elmhes;

  procedure hqr ( nm,n,low,igh : in integer32; H : in out Matrix;
                  wr,wi : out Quad_Double_Vectors.Vector;
                  ierr : out integer32 ) is

  -- NOTE : the fortran code contained the following reference:
  --   this subroutine is a translation of the algol procedure hqr,
  --   num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
  --   handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).

    notlas : boolean;
    jj,k,L,m,en,na,itn,its,mp2,enm2 : integer32;
    p,q,r,s,t,w,x,y,zz,norm,tst1,tst2 : quad_double;
    zero : constant quad_double := create(0.0);

  begin
    ierr := 0; norm := zero;
    k := 1; -- store roots isolated by balanc and compute matrix norm 
    for i in 1..n loop
      for j in k..n loop
        norm := norm + abs(H(i,j));
      end loop;
      k := i;
      if ((i < low) or (i > igh)) then
        wr(i) := H(i,i); wi(i) := zero;
      end if;
    end loop;
    en := igh; t := zero; itn := 30*n;
    loop -- main outer loop
      if (en < low) then return; end if; -- search for next eigenvalues
      its := 0; na := en - 1; enm2 := na - 1;
     -- look for single small sub-diagonal element for l=en step -1 until low
      loop -- main inner loop
        for LL in low..en loop
          L := en + low - LL;
          exit when (L = low);
          s := abs(H(L-1,L-1)) + abs(H(L,L));
          if (s = zero) then s := norm; end if;
          tst1 := s;
          tst2 := tst1 + abs(H(L,L-1));
          exit when (tst2 = tst1);
        end loop;
        x := H(en,en); -- form shift
        if (L = en) then -- one root found
          wr(en) := x + t; wi(en) := zero; en := na;
          exit; -- leave main inner loop
        else
          y := H(na,na);
          w := H(en,na) * H(na,en);
          if (L = na) then -- two roots found
            p := (y - x) / 2.0;
            q := p * p + w;
            zz := sqrt(abs(q));
            x := x + t;
            if (q < 0.0) then
              wr(na) := x + p; wi(na) := zz;     -- complex pair
              wr(en) := x + p; wi(en) := -zz;
            else
              zz := p + dsign(zz,p);  -- real pair
              wr(na) := x + zz; wi(na) := zero;
              wr(en) := wr(na); wi(en) := zero;
              if (zz /= zero) then wr(en) := x - w / zz; end if;
            end if;
            en := enm2;
            exit; -- leave main inner loop
          else
            if (itn = 0) then  -- set error: all eigenvalues have not
	      ierr := en;      --   converged after 30*n iterations
	      return;
            end if;
            if ((its = 10) or (its = 20)) then
              t := t + x; -- form exceptional shift
              for i in low..en loop
                H(i,i) := H(i,i) - x;
              end loop;
              s := abs(H(en,na)) + abs(H(na,enm2));
              x := 0.75 * s; y := x;
              w := -0.4375 * s * s;
            end if;
            its := its + 1; itn := itn - 1;
           -- look for two consecutive small sub-diagonal elements
           --  for m=en-2 step -1 until l do --
            for mm in L..enm2 loop
              m := enm2 + L - mm;
              zz := H(m,m); r := x - zz; s := y - zz;
              p := (r * s - w) / H(m+1,m) + H(m,m+1);
              q := H(m+1,m+1) - zz - r - s;
              r := H(m+2,m+1);
              s := abs(p) + abs(q) + abs(r);
              p := p / s; q := q / s; r := r / s;
              exit when (m = L);
              tst1 := abs(p)*(abs(H(m-1,m-1)) + abs(zz) + abs(H(m+1,m+1)));
              tst2 := tst1 + abs(H(m,m-1))*(abs(q) + abs(r));
              exit when (tst2 = tst1);
            end loop;
            mp2 := m + 2;
            for i in mp2..en loop
              H(i,i-2) := zero;
              if (i /= mp2) 
               then H(i,i-3) := zero;
              end if;
            end loop;
           -- double qr step involving rows l to en and columns m to en
            for k in m..na loop
              notlas := (k /= na);
              if k /= m then
                p := H(k,k-1); q := H(k+1,k-1); r := zero;
                if notlas then r := H(k+2,k-1); end if;
                x := abs(p) + abs(q) + abs(r);
              end if;
              if ((k = m) or else (x /= zero)) then
                if k /= m then
                  p := p / x; q := q / x; r := r / x;
                end if;
                s := dsign(sqrt(p*p+q*q+r*r),p);
                if (k = m) then
                  if (L /= m) then H(k,k-1) := -H(k,k-1); end if;
                else
                  H(k,k-1) := -s*x;
                end if;
                p := p + s; x := p / s; y := q / s;
                zz := r / s; q := q / p; r := r / p;
                if (not notlas) then 
                  for j in k..en loop -- row modification
                    p := H(k,j) + q * H(k+1,j);
                    H(k,j) := H(k,j) - p * x;
                    H(k+1,j) := H(k+1,j) - p * y;
                  end loop;
                  jj := min0(en,k+3);
                  for i in L..jj loop -- column modification
                    p := x * H(i,k) + y * H(i,k+1);
                    H(i,k) := H(i,k) - p;
                    H(i,k+1) := H(i,k+1) - p * q;
                  end loop;
                else
                  for j in k..en loop -- row modification
                    p := H(k,j) + q * H(k+1,j) + r * H(k+2,j);
                    H(k,j) := H(k,j) - p * x;
                    H(k+1,j) := H(k+1,j) - p * y;
                    H(k+2,j) := H(k+2,j) - p * zz;
                  end loop;
                  jj := min0(en,k+3);
                  for i in L..jj loop  -- column modification
                    p := x * H(i,k) + y * H(i,k+1) + zz * H(i,k+2);
                    H(i,k) := H(i,k) - p;
                    H(i,k+1) := H(i,k+1) - p * q;
                    H(i,k+2) := H(i,k+2) - p * r;
                  end loop;
                end if;
              end if;
            end loop;
          end if;
        end if;
      end loop;
    end loop; 
  end hqr;

  procedure eltran ( nm,n,low,igh : in integer32; A : in Matrix;
                     int : in Standard_Integer_Vectors.Vector;
                     z : out Matrix ) is

  -- NOTE : the fortran code contained the following reference
  --   this subroutine is a translation of the algol procedure elmtrans,
  --   num. math. 16, 181-204(1970) by peters and wilkinson.
  --   handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).

    kl,mp,mp1,ii : integer32;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for j in 1..n loop       -- initialize z to the identity matrix
      for i in 1..n loop
        z(i,j) := zero;
      end loop;
      z(j,j) := one;
    end loop;
    kl := igh - low - 1;
    if kl >= 1 then
      for mm in 1..kl loop   -- for mp=igh-1 step -1 until low+1 do 
        mp := igh - mm;
        mp1 := mp + 1;
        for i in mp1..igh loop
          z(i,mp) := a(i,mp-1);
        end loop;
        ii := int(mp);
        if ii /= mp then
          for j in mp..igh loop
            z(mp,j) := z(ii,j);
            z(ii,j) := zero;
          end loop;
          z(ii,mp) := one;
        end if;
      end loop;
    end if;
  end eltran;

  procedure hqr2 ( nm,n,low,igh : in integer32; H : in out Matrix;
                   wr,wi : out Quad_Double_Vectors.Vector;
                   z : in out Matrix; ierr : out integer32 ) is

  -- NOTE : the fortran code contained the following reference:
  --   this subroutine is a translation of the algol procedure hqr2,
  --   num. math. 16, 181-204(1970) by peters and wilkinson.
  --   handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).

    ind,jj,k,L,m,en,na,itn,its,mp2,enm2 : integer32;
    p,q,r,s,t,w,x,y,ra,sa,vr,vi,zz,norm,tst1,tst2 : quad_double;
    notlas : boolean;
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    ierr := 0;
    norm := zero;
    k := 1;
    for i in 1..n loop     -- store roots isolated by balanc
      for j in k..n loop   -- and compute matrix norm
        norm := norm + abs(H(i,j));
      end loop;
      k := i;
      if ((i < low) or (i > igh)) then
        wr(i) := H(i,i);
        wi(i) := zero;
      end if;
    end loop;
    en := igh; t := zero; itn := 30*n;
    loop -- main outer loop
      exit when (en <low); -- search for next eigenvalues
      its := 0; na := en - 1; enm2 := na - 1;
      loop -- main inner loop
        -- look for single small sub-diagonal element
        for ll in low..en loop --  for l=en step -1 until low do
          L := en + low - ll;
          exit when (L = low);
          s := abs(H(L-1,L-1)) + abs(H(L,L));
          if (s = zero) then s := norm; end if;
          tst1 := s;
          tst2 := tst1 + abs(H(L,L-1));
          exit when (tst2 = tst1);
        end loop;
        x := H(en,en);  -- form shift
        if (L = en) then
          H(en,en) := x + t;    -- one root found
          wr(en) := H(en,en);
          wi(en) := zero; en := na;
          exit;                 -- leave main inner loop
        else
          y := H(na,na);
          w := H(en,na) * H(na,en);
          if (L = na) then      -- two roots found
            p := (y - x) / 2.0;
            q := p * p + w;
            zz := sqrt(abs(q));
            H(en,en) := x + t;
            x := H(en,en);
            H(na,na) := y + t;
            if (q >= zero) then
              zz := p + dsign(zz,p);  -- real pair
              wr(na) := x + zz;
              wr(en) := wr(na);
              if (zz /= zero) then wr(en) := x - w / zz; end if;
              wi(na) := zero;
              wi(en) := zero;
              x := H(en,na);
              s := abs(x) + abs(zz);
              p := x / s;
              q := zz / s;
              r := sqrt(p*p+q*q);
              p := p / r;
              q := q / r;
              for j in na..n loop -- row modification
                zz := H(na,j);
                H(na,j) := q * zz + p * H(en,j);
                H(en,j) := q * H(en,j) - p * zz;
              end loop;
              for i in 1..en loop -- column modification
                zz := H(i,na);
                H(i,na) := q * zz + p * H(i,en);
                H(i,en) := q * H(i,en) - p * zz;
              end loop;
              for i in low..igh loop  -- accumulate transformations
                zz := z(i,na);
                z(i,na) := q * zz + p * z(i,en);
                z(i,en) := q * z(i,en) - p * zz;
              end loop;
            else -- complex pair
              wr(na) := x + p; wr(en) := x + p;
              wi(na) := zz; wi(en) := -zz;
            end if;
            en := enm2;
            exit;  -- leave main inner loop
          else
            if (itn = 0) then  -- set error: all eigenvalues have not
              ierr := en;      --   converged after 30*n iterations
              return;
            end if;
            if ((its = 10) or (its = 20)) then
              t := t + x;  -- form exceptional shift
              for i in low..en loop
                H(i,i) := H(i,i) - x;
              end loop;
              s := abs(H(en,na)) + abs(H(na,enm2));
              x := 0.75 * s;
              y := x;
              w := -0.4375 * s * s;
            end if;
            its := its + 1;
            itn := itn - 1;
          -- look for two consecutive small sub-diagonal elements.
            for mm in L..enm2 loop  -- for m=en-2 step -1 until l do
              m := enm2 + L - mm;
              zz := H(m,m);
              r := x - zz;
              s := y - zz;
              p := (r * s - w) / H(m+1,m) + H(m,m+1);
              q := H(m+1,m+1) - zz - r - s;
              r := H(m+2,m+1);
              s := abs(p) + abs(q) + abs(r);
              p := p / s;
              q := q / s;
              r := r / s;
              exit when (m = L);
              tst1 := abs(p)*(abs(H(m-1,m-1)) + abs(zz) + abs(H(m+1,m+1)));
              tst2 := tst1 + abs(H(m,m-1))*(abs(q) + abs(r));
              exit when (tst2 = tst1);
            end loop;
            mp2 := m + 2;
            for i in mp2..en loop
              H(i,i-2) := zero;
              if (i /= mp2) then H(i,i-3) := zero; end if;
            end loop;
           -- double qr step involving rows l to en and columns m to en
            for k in m..na loop
              notlas := (k /= na);
              if k /= m then
                p := H(k,k-1);
                q := H(k+1,k-1);
                r := zero;
                if (notlas) then r := H(k+2,k-1); end if;
                x := abs(p) + abs(q) + abs(r);
              end if;
              if ((k = m) or else (x /= zero)) then
                if k /= m then
                  p := p / x; q := q / x; r := r / x;
                end if;
                s := dsign(sqrt(p*p+q*q+r*r),p);
                if (k = m) then
                  if (L /= m) then H(k,k-1) := -H(k,k-1); end if;
                else
                  H(k,k-1) := -s * x;
                end if;
                p := p + s; x := p / s; y := q / s;
                zz := r / s; q := q / p; r := r / p;
                if (not notlas) then
                  for j in k..n loop  -- row modification
                    p := H(k,j) + q * H(k+1,j);
                    H(k,j) := H(k,j) - p * x;
                    H(k+1,j) := H(k+1,j) - p * y;
                  end loop;
                  jj := min0(en,k+3);
                  for i in 1..jj loop  -- column modification
                    p := x * H(i,k) + y * H(i,k+1);
                    H(i,k) := H(i,k) - p;
                    H(i,k+1) := H(i,k+1) - p * q;
                  end loop;
                  for i in low..igh loop  -- accumulate transformations
                    p := x * z(i,k) + y * z(i,k+1);
                    z(i,k) := z(i,k) - p;
                    z(i,k+1) := z(i,k+1) - p * q;
                  end loop;
                else
                  for j in k..n loop  -- row modification
                    p := H(k,j) + q * H(k+1,j) + r * H(k+2,j);
                    H(k,j) := H(k,j) - p * x;
                    H(k+1,j) := H(k+1,j) - p * y;
                    H(k+2,j) := H(k+2,j) - p * zz;
                  end loop;
                  jj := min0(en,k+3);
                  for i in 1..jj loop -- column modification
                    p := x * H(i,k) + y * H(i,k+1) + zz * H(i,k+2);
                    H(i,k) := H(i,k) - p;
                    H(i,k+1) := H(i,k+1) - p * q;
                    H(i,k+2) := H(i,k+2) - p * r;
                  end loop;
                  for i in low..igh loop -- accumulate transformations
                    p := x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2);
                    z(i,k) := z(i,k) - p;
                    z(i,k+1) := z(i,k+1) - p * q;
                    z(i,k+2) := z(i,k+2) - p * r;
                  end loop;
                end if;
              end if;
            end loop;
          end if;
        end if;
      end loop; -- end inner loop
    end loop; -- end main loop
   -- all roots found, backsubstitute to find vectors of upper triangular form 
   -- put_line("after main loop, before backsubstitute, H is "); put(H);
   -- put_line("after main loop, before backsubstitute, z is "); put(z);
    if norm /= zero then
      for nn in 1..n loop   --  for en=n step -1 until 1 do
        en := n + 1 - nn;
        p := wr(en); q := wi(en);
        na := en - 1;
        if q < zero then -- complex vector
          m := na; -- last vector component chosen imaginary so that
                   --   eigenvector matrix is triangular 
          if (abs(H(en,na)) > abs(H(na,en))) then
            H(na,na) := q / H(en,na);
            H(na,en) := -(H(en,en) - p) / H(en,na);
          else
            cdiv(zero,-H(na,en),H(na,na)-p,q,H(na,na),H(na,en));
          end if;
          H(en,na) := zero; H(en,en) := one;
          enm2 := na - 1;
          if (enm2 /= 0) then
            for ii in 1..enm2 loop -- for i=en-2 step -1 until 1 do
              ind := na - ii;
              w := H(ind,ind) - p;
              ra := zero; sa := zero;
              for j in m..en loop
                ra := ra + H(ind,j) * H(j,na);
                sa := sa + H(ind,j) * H(j,en);
              end loop;
              if (wi(ind) < zero) then
                zz := w; r := ra; s := sa;
              else
                m := ind;
                if (wi(ind) = zero) then
                  cdiv(-ra,-sa,w,q,H(ind,na),H(ind,en));
                else
                  x := H(ind,ind+1);  -- solve complex equations
                  y := H(ind+1,ind);
                  vr := (wr(ind) - p)*(wr(ind) - p) + wi(ind)*wi(ind) - q*q;
                  vi := (wr(ind) - p) * 2.0 * q;
                  if ((vr = zero) and (vi = zero)) then
                    tst1 := norm*(abs(w)+abs(q)+abs(x)+abs(y)+abs(zz));
                    vr := tst1;
                    loop
                      vr := create(0.01)*vr;
                      tst2 := tst1 + vr;
                      exit when (tst2 <= tst1);
                    end loop;
                  end if;
                  cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,
                       H(ind,na),H(ind,en));
                  if (abs(x) > abs(zz) + abs(q)) then
                    H(ind+1,na) := (-ra - w * H(ind,na) + q * H(ind,en)) / x;
                    H(ind+1,en) := (-sa - w * H(ind,en) - q * H(ind,na)) / x;
                  else
                    cdiv(-r-y*H(ind,na),-s-y*H(ind,en),zz,q,
                         H(ind+1,na),H(ind+1,en));
                  end if;
                end if; -- overflow control
                t := dmax1(abs(H(ind,na)),abs(H(ind,en)));
                if (t /= zero) then
                  tst1 := t;
                  tst2 := tst1 + 1.0/tst1;
                  if (tst2 <= tst1) then
                    for j in ind..en loop
                      H(j,na) := H(j,na)/t;
                      H(j,en) := H(j,en)/t;
                    end loop;
                  end if;
                end if;
              end if;
            end loop; -- end complex vector 
          end if;
        elsif q = zero then -- real vector
         -- put_line("in real vector, H is"); put(H);
          m := en;
          H(en,en) := one;
         -- put("setting H("); put(en,1); put(","); put(en,1);
         -- put_line(") to one");
          if (na /= 0) then
            for ii in 1..na loop -- for i=en-1 step -1 until 1 do
              ind := en - ii;
              w := H(ind,ind) - p;
              r := zero;
              for j in m..en loop
                r := r + H(ind,j) * H(j,en);
              end loop;
              if (wi(ind) < zero) then
                zz := w;
                s := r;
              else
                m := ind;
               -- put("wi("); put(ind,1); put(") = "); put(wi(ind)); new_line;
                if (wi(ind) /= zero) then -- solve real equations
                 -- put_line("in solve real equations");
                  x := H(ind,ind+1);
                  y := H(ind+1,ind);
                  q := (wr(ind) - p) * (wr(ind) - p) + wi(ind) * wi(ind);
                  t := (x * s - zz * r) / q;
                  H(ind,en) := t;
                  if (abs(x) > abs(zz)) then
                    H(ind+1,en) := (-r - w * t) / x;
                  else
                    H(ind+1,en) := (-s - y * t) / zz;
                  end if;
                else
                  t := w;
                  if (t = zero) then
                    tst1 := norm;
                    t := tst1;
                    loop
                      t := create(0.01) * t;
                      tst2 := norm + t;
                      exit when (tst2 <= tst1);
                    end loop;
                  end if;
                  H(ind,en) := -r / t;
                end if;
                t := abs(H(ind,en)); -- overflow control 
                if t /= zero then
                  tst1 := t;
                  tst2 := tst1 + 1.0/tst1;
                  if tst2 <= tst1 then
                    for j in ind..en loop
                      H(j,en) := H(j,en)/t;
                    end loop;
                  end if;
                end if;
              end if;
            end loop;
          end if;
         -- put_line("at the end of real vector, H is"); put(H);
        end if;
      end loop;
     -- end back substitution, vectors of isolated roots
      for i in 1..n loop
        if ((i < low) or (i > igh)) then
          for j in 1..n loop
            z(i,j) := H(i,j);
          end loop;
        end if;
      end loop;
  -- multiply by transformation matrix to give vectors of original full matrix
      for j in low..n loop --  for j=n step -1 until low do
        jj := n + low - j;
        m := min0(jj,igh);
        for i in low..igh loop
          zz := zero;
          for k in low..m loop
            zz := zz + z(i,k) * H(k,jj);
          end loop;
          z(i,jj) := zz;
        end loop;
      end loop;
    end if;
  end hqr2;

  procedure balbak ( nm,n,low,igh : in integer32;
                     scale : in Quad_Double_Vectors.Vector;
                     m : in integer32; z : in out Matrix ) is

  -- NOTE : the fortran code contained the following reference
  --   this subroutine is a translation of the algol procedure balbak,
  --   num. math. 13, 293-304(1969) by parlett and reinsch.
  --   handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

    s : quad_double;
    ind,k : integer32;

  begin
    if m > 0 then
      if igh > low then
        for i in low..igh loop -- left hand eigenvectors are back transformed
          s := scale(i);       -- if the foregoing statement is replaced by
          for j in 1..m loop   --     s = 1.0d0/scale(i).
            z(i,j) := z(i,j)*s;
          end loop;
        end loop;
      end if;
      for i in 1..n loop       -- for i=low-1 step -1 until 1, 
        ind := i;              -- igh+1 step 1 until n do
        if ind < low or ind > igh then
          if (ind < low) then ind := low - i; end if;
          k := integer32(hihi_part(scale(ind)));
          if (k /= ind) then
            for j in 1..m loop
              s := z(ind,j);
              z(ind,j) := z(k,j);
              z(k,j) := s;
            end loop;
          end if;
        end if;
      end loop;
    end if;
  end balbak;

-- AUXILIARY ROUTINE for hqr2 :

  procedure cdiv ( ar,ai,br,bi : in quad_double;
                   cr,ci : out quad_double ) is

    s,ars,ais,brs,bis : quad_double;

  begin
    s := abs(br) + abs(bi);
    ars := ar/s;
    ais := ai/s;
    brs := br/s;
    bis := bi/s;
    s := sqr(brs) + sqr(bis);    -- brs**2 + bis**2;
    cr := (ars*brs + ais*bis)/s;
    ci := (ais*brs - ars*bis)/s;
  end cdiv;

end Quad_Double_Eigenvalues;
