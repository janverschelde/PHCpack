with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Mathematical_Functions;
with Standard_Vector_Splitters;
with Standard_Matrix_Splitters;
with Standard_Complex_BLAS_Helpers;
with Standard_Inlined_BLAS_Helpers;

package body Standard_Inlined_Singular_Values is

  function Min0 ( a,b : integer32 ) return integer32 is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Min0;

  procedure SVD ( xrv,xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  n,p : in integer32;
                  sr,si : in Standard_Floating_Vectors.Link_to_Vector;
                  er,ei : in Standard_Floating_Vectors.Link_to_Vector;
                  urv,uiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  vrv,viv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  job : in integer32; info : out integer32;
                  rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector ) is

    maxit : constant integer32 := 30; -- maximum number of iterations
    jobu : constant integer32 := (job mod 100)/10; -- what to compute?
    wantu : constant boolean := (jobu /= 0);       -- want U matrix
    wantv : constant boolean := (job mod 10 /= 0); -- want V matrix
    ncu : integer32 := n;
    nct : constant integer32 := Min0(n-1,p);
    nrt : constant integer32
        := Standard_Complex_BLAS_Helpers.Max0(0,Min0(p-2,n));
    lu : constant integer32
       := Standard_Complex_BLAS_Helpers.Max0(nct,nrt);
    lp1,m,nctp1,nrtp1,ll,lm1,mm,iter,kase,ls,mp1,mm1,kk : integer32;
    b,f,g,tr,ti,zr,zi,test,ztest,t1,cs,sn : double_float;
    scale,sm,smm1,emm1,sl,el,c,shift : double_float;
    cntkase2 : integer32 := 0; -- patch to fix with NaN inputs...

    use Standard_Inlined_BlAS_Helpers;

  begin
    info := 0;
    if jobu > 1
     then ncu := min0(n,p);
    end if;
   -- reduce (xrv, xiv) to bidiagonal form, store diagonal elements in s
   -- and store the super diagonal elements in e
    if lu >= 1 then
      for ell in 1..lu loop
        lp1 := ell+1;
        if ell <= nct then
          -- compute the transformation for the ell-th column
          -- and place the ell-th diagonal in s(l)
          sr(ell) := dznrm2(n-ell+1,xrv,xiv,ell,ell,1);
          si(ell) := 0.0;
          if abs(sr(ell)) /= 0.0 then
            if cdabs(xrv(ell)(ell),xiv(ell)(ell)) /= 0.0 then
              csign(sr(ell),si(ell),xrv(ell)(ell),xiv(ell)(ell),
                    sr(ell),si(ell));
            end if;
            -- compute 1/s(ell)
            b := sr(ell)*sr(ell) + si(ell)*si(ell);
            b := 1.0/b; f := b*sr(ell); g := -b*si(ell);
            zscal(n-ell+1,f,g,xrv,xiv,ell,ell,1);
            xrv(ell)(ell) := 1.0 + xrv(ell)(ell);
          end if;
          sr(ell) := -sr(ell); si(ell) := -si(ell);
        end if;
        if p >= lp1 then
          for j in lp1..p loop
            if ell <= nct then
              if (abs(sr(ell)) + abs(si(ell))) /= 0.0 then
               -- apply the transformation
               -- t := -zdotc(n-l+1,x,l,l,1,x,l,j,1)/x(l,l);
                zdotc(n-ell+1,xrv,xiv,ell,ell,1,xrv,xiv,ell,j,1,tr,ti);
               -- divide t by x(ell)(ell) and flip sign
                f := xrv(ell)(ell); g := xiv(ell)(ell);
                b := f*f + g*g; zr := tr*f + ti*g; zi := ti*f - tr*g;
                tr := -zr/b; ti := -zi/b;
                zaxpy(n-ell+1,tr,ti,xrv,xiv,ell,ell,1,xrv,xiv,ell,j,1);
              end if;
              -- place the ell-th row of x into e for subsequent
              -- calculation of the row transformation
              er(j) := xrv(j)(ell); -- e(j) = Conjugate(x(l,j))
              ei(j) := -xiv(j)(ell);
            end if;
          end loop;
        end if;
        if wantu and ell <= nct then
         -- place the transformation in u for subsequent back multiplication
          for i in ell..n loop -- u(i,l) := x(i,l);
            urv(ell)(i) := xrv(ell)(i);
            uiv(ell)(i) := xiv(ell)(i);
          end loop;
        end if;
        if ell <= nrt then
           -- compute the ell-th row transformation
           -- and place the ell-th super diagonal in e(ell)
          er(ell) := dznrm2(p-ell,er,ei,lp1,1);
          ei(ell) := 0.0;
          if abs(er(ell)) /= 0.0 then
            if cdabs(er(lp1),ei(lp1)) /= 0.0
             then csign(er(ell),ei(ell),er(lp1),ei(lp1),er(ell),ei(ell));
            end if;
           -- compute 1/e(ell)
            b := er(ell)*er(ell) + ei(ell)*ei(ell);
            b := 1.0/b; f := b*er(ell); g := -b*ei(ell);
            zscal(p-ell,f,g,er,ei,lp1,1);
            er(lp1) := 1.0 + er(lp1);
          end if;
          er(ell) := -er(ell);
          ei(ell) := ei(ell); -- e(l) := -Conjugate(e(l));
          if lp1 <= n and ((abs(er(ell)) + abs(ei(ell))) /= 0.0) then
             -- apply the transformation
            for i in lp1..n loop
              rwrk(i) := 0.0; iwrk(i) := 0.0;
            end loop;
            for j in lp1..p loop
              zaxpy(n-ell,er(j),ei(j),xrv,xiv,lp1,j,1,rwrk,iwrk,lp1,1);
            end loop;
            for j in lp1..p loop
             -- zaxpy(n-l,Conjugate(-e(j)/e(lp1)),work,lp1,1,x,lp1,j,1);
             -- compute Conjugate(-e(j)/e(lp1)) and store in zr and zi
              f := er(lp1); g := ei(lp1); b := f*f + g*g;
              tr := -er(j); ti := -ei(j); -- take -e(j)
              zr := tr*f + ti*g; zi := ti*f - tr*g;
              zr := zr/b; zi := -zi/b; -- take conjugate
              zaxpy(n-ell,zr,zi,rwrk,iwrk,lp1,1,xrv,xiv,lp1,j,1);
            end loop;
          end if;
          if wantv then
           -- place the transformation in v for subsequent back multiplication
            for i in lp1..p loop -- v(i,l) := e(i);
              vrv(ell)(i) := er(i);
              viv(ell)(i) := ei(i);
            end loop;
          end if;
        end if;
      end loop;
    end if;
   -- set up the final bidiagonal matrix of order m
    m := min0(p,n+1);
    nctp1 := nct+1;
    nrtp1 := nrt+1;
    if nct < p then
      sr(nctp1) := xrv(nctp1)(nctp1);
      si(nctp1) := xiv(nctp1)(nctp1);
    end if;
    if n < m
     then sr(m) := 0.0; si(m) := 0.0;
    end if;
    if nrtp1 < m then
      er(nrtp1) := xrv(m)(nrtp1);
      ei(nrtp1) := xiv(m)(nrtp1);
    end if;
    er(m) := 0.0; ei(m) := 0.0;
   -- if required, generate u
    if wantu then
      if ncu >= nctp1 then
        for j in nctp1..ncu loop
          for i in 1..n loop 
            urv(j)(i) := 0.0; uiv(j)(i) := 0.0; -- u(i,j) := Create(0.0);
          end loop;
          urv(j)(j) := 1.0; uiv(j)(j) := 0.0;   -- u(j,j) := Create(1.0);
        end loop;
      end if;
      if nct >= 1 then
        for lidx in 1..nct loop
          ll := nct - lidx + 1;
          if (abs(sr(ll)) + abs(si(ll))) = 0.0 then
            for i in 1..n loop
              urv(ll)(i) := 0.0; uiv(ll)(i) := 0.0; -- u(i,ll) := Create(0.0);
            end loop;
            urv(ll)(ll) := 1.0; uiv(ll)(ll) := 0.0; -- u(ll,ll) := Create(1.0);
          else
            lp1 := ll + 1;
            if ncu >= lp1 then
              for j in lp1..ncu loop
               -- t := -zdotc(n-ll+1,u,ll,ll,1,u,ll,j,1)/u(ll,ll);
                zdotc(n-ll+1,urv,uiv,ll,ll,1,urv,uiv,ll,j,1,tr,ti);
               -- divide t by urv(ll)(ll) and flip sign
                f := urv(ll)(ll); g := uiv(ll)(ll);
                b := f*f + g*g; zr := tr*f + ti*g; zi := ti*f - tr*g;
                tr := -zr/b; ti := -zi/b;
               -- zaxpy(n-ll+1,t,u,ll,ll,1,u,ll,j,1);
                zaxpy(n-ll+1,tr,ti,urv,uiv,ll,ll,1,urv,uiv,ll,j,1);
              end loop;
            end if;
            zscal(n-ll+1,-1.0,0.0,urv,uiv,ll,ll,1);
            urv(ll)(ll) := 1.0 + urv(ll)(ll);
            lm1 := ll - 1;
            if lm1 >= 1 then
              for i in 1..lm1 loop -- u(i,ll) := Create(0.0);
                urv(ll)(i) := 0.0; uiv(ll)(i) := 0.0;
              end loop;
            end if;
          end if;
        end loop;
      end if;
    end if;
   -- if required, generate v
    if wantv then
      for lidx in 1..p loop
        ll := p - lidx + 1;
        lp1 := ll + 1;
        if ll <= nrt then
          if (abs(er(ll)) + abs(ei(ll))) /= 0.0 then
            for j in lp1..p loop
             -- t := -zdotc(p-ll,v,lp1,ll,1,v,lp1,j,1)/v(lp1,ll);
              zdotc(p-ll,vrv,viv,lp1,ll,1,vrv,viv,lp1,j,1,tr,ti);
             -- divide t by v(ll)(lp1) and flip sign
              f := vrv(ll)(lp1); g := viv(ll)(lp1);
              b := f*f + g*g; zr := tr*f + ti*g; zi := ti*f - tr*g;
              tr := -zr/b; ti := -zi/b;
             -- zaxpy(p-ll,t,v,lp1,ll,1,v,lp1,j,1);
              zaxpy(p-ll,tr,ti,vrv,viv,lp1,ll,1,vrv,viv,lp1,j,1);
            end loop;
          end if;
        end if;
        for i in 1..p loop
          vrv(ll)(i) := 0.0; viv(ll)(i) := 0.0; -- v(i,ll) := Create(0.0);
        end loop;
        vrv(ll)(ll) := 1.0; viv(ll)(ll) := 0.0; -- v(ll,ll) := Create(1.0);
      end loop;
    end if;
   -- transform s and e so that they are double precision
    for i in 1..m loop
      if (abs(sr(i)) + abs(si(i))) /= 0.0 then
        tr := cdabs(sr(i),si(i));       -- t := Create(cdabs(s(i)),0.0);
        zr := sr(i)/tr; zi := si(i)/tr; -- r := s(i)/t;
        sr(i) := tr; si(i) := 0.0;      -- s(i) := t;
        if i < m then -- e(i) := e(i)/r;
          b := zr*zr + zi*zi; f := er(i); g := ei(i);
          er(i) := (f*zr + g*zi)/b;
          ei(i) := (g*zr - f*zi)/b;
        end if;
        if wantu
         then zscal(n,zr,zi,urv,uiv,1,i,1);
        end if;
      end if;
      exit when (i = m);
      if (abs(er(i)) + abs(ei(i))) /= 0.0 then
        tr := cdabs(er(i),ei(i));   -- Create(cdabs(e(i)),0.0);
        f := er(i); g := ei(i); b := f*f + g*g; b := tr/b;
        zr := b*f; zi := -b*g;      -- r := t/e(i);
        er(i) := tr;  ei(i) := 0.0; -- e(i) := t;
        f := sr(i+1); g := si(i+1);
        sr(i+1) := f*zr - g*zi;     -- s(i+1) := s(i+1)*r;
        si(i+1) := f*zi + g*zr;
        if wantv
         then zscal(p,zr,zi,vrv,viv,1,i+1,1);
        end if;
      end if;
    end loop;
   -- main iteration loop for the singular values
    mm := m; iter := 0;
    loop
      exit when m = 0; -- quit if all the singular values have been found
      if iter > maxit  -- too many iterations have been performed
       then info := m; exit; -- set the flag and exit the loop
      end if;
     -- This section of the program inspects for negligible elements in 
     -- the s and e arrays.  On completion the variables kase and l are
     -- set as follows:
     --   kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
     --   kase = 2     if s(l) is negligible and l.lt.m
     --   kase = 3     if e(l-1) is negligible, l.lt.m, and
     --                  s(l), ..., s(m) are not negligible (qr step).
     --   kase = 4     if e(m-1) is negligible (convergence).
      for ell in 1..m loop
        ll := m - ell;
        exit when (ll = 0);
        test := cdabs(sr(ll),si(ll)) + cdabs(sr(ll+1),si(ll+1));
        ztest := test + cdabs(er(ll),ei(ll));
        if ztest = test
         then er(ll) := 0.0; ei(ll) := 0.0; exit;
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
          test := 0.0;
          if ls /= n
           then test := test + cdabs(er(ls),ei(ls));
          end if;
          if ls /= ll+1
           then test := test + cdabs(er(ls-1),ei(ls-1));
          end if;
          ztest := test + cdabs(sr(ls),si(ls));
          if ztest = test 
           then sr(ls) := 0.0; si(ls) := 0.0; exit;
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
          f := er(m-1); er(m-1) := 0.0; ei(m-1) := 0.0;
          for k in ll..mm1 loop
            kk := mm1 - k + ll;
            t1 := sr(kk); cs := 0.0; sn := 0.0;
            Standard_Complex_BLAS_Helpers.drotg(t1,f,cs,sn);
            sr(kk) := t1; si(kk) := 0.0;
            if kk /= ll then
              f := -sn*er(kk-1);
              er(kk-1) := cs*er(kk-1);
              ei(kk-1) := cs*ei(kk-1);
            end if;
            if wantv
             then zdrot(p,vrv,viv,1,kk,1,vrv,viv,1,m,1,cs,sn);
            end if;
          end loop;
        when 2 => -- split at negligible s(ll)
          f := er(ll-1); er(ll-1) := 0.0; ei(ll-1) := 0.0;
          for k in ll..m loop
            t1 := sr(k);
            Standard_Complex_BLAS_Helpers.drotg(t1,f,cs,sn);
            sr(k) := t1; si(k) := 0.0;
            f := -sn*er(k);
            er(k) := cs*er(k);
            ei(k) := cs*ei(k);
            if wantu
             then zdrot(n,urv,uiv,1,k,1,urv,uiv,1,ll-1,1,cs,sn);
            end if;
          end loop;
          cntkase2 := cntkase2 + 1;  -- patch for fixing NaN inputs ...
          if cntkase2 > maxit then return; end if;
        when 3 => -- perform one qr step
          -- 1) calculate the shift
           scale := Standard_Complex_BLAS_Helpers.dmax1
                      (cdabs(sr(m),si(m)),cdabs(sr(m-1),si(m-1)),
                       cdabs(er(m-1),ei(m-1)),cdabs(sr(ll),si(ll)),
                       cdabs(er(ll),ei(ll)));
           sm := sr(m)/scale;
           smm1 := sr(m-1)/scale;
           emm1 := er(m-1)/scale;
           sl := sr(ll)/scale;
           el := er(ll)/scale;
           b := ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0;
           c := (sm*emm1)**2;
           shift := 0.0;
           if b = 0.0 or c = 0.0 then
             shift := Standard_Mathematical_Functions.SQRT(b**2+c);
             if b < 0.0
              then shift := -shift;
             end if;
             shift := c/(b + shift);
           end if;
           f := (sl + sm)*(sl - sm) + shift;
           g := sl*el;
          -- 2) chase zeros
           mm1 := m - 1;
           for k in ll..mm1 loop
             Standard_Complex_BLAS_Helpers.drotg(f,g,cs,sn);
             if k /= ll 
              then er(k-1) := f; ei(k-1) := 0.0;
             end if;
             f := cs*sr(k) + sn*er(k);
             er(k) := cs*er(k) - sn*sr(k);
             ei(k) := cs*ei(k) - sn*si(k);
             g := sn*sr(k+1);
             sr(k+1) := cs*sr(k+1);
             si(k+1) := cs*si(k+1);
             if wantv
              then zdrot(p,vrv,viv,1,k,1,vrv,viv,1,k+1,1,cs,sn);
             end if;
             Standard_Complex_BLAS_Helpers.drotg(f,g,cs,sn);
             sr(k) := f; si(k) := 0.0;
             f := cs*er(k) + sn*sr(k+1);
             sr(k+1) := -sn*er(k) + cs*sr(k+1);
             si(k+1) := -sn*ei(k) + cs*si(k+1);
             g := sn*er(k+1);
             er(k+1) := cs*er(k+1);
             ei(k+1) := cs*ei(k+1);
             if wantu and k < n
              then zdrot(n,urv,uiv,1,k,1,urv,uiv,1,k+1,1,cs,sn);
             end if;
           end loop;
           er(m-1) := f; ei(m-1) := 0.0;
           iter := iter + 1;
        when 4 => -- convergence
         -- 1) make the singular value positive
          if sr(ll) < 0.0 then
            sr(ll) := -sr(ll); si(ll) := -si(ll);
            if wantv 
             then zscal(p,-1.0,0.0,vrv,viv,1,ll,1);
            end if;
          end if;
         -- 2) order the singular values
          while ll /= mm loop
            exit when (sr(ll) >= sr(ll+1));
            tr := sr(ll);       ti := si(ll);
            sr(ll) := sr(ll+1); si(ll) := si(ll+1);
            sr(ll+1) := tr;     si(ll+1) := ti;
            if wantv and ll < p
             then zswap(p,vrv,viv,1,ll,1,vrv,viv,1,ll+1,1);
            end if;
            if wantu and ll < n
             then zswap(n,urv,uiv,1,ll,1,urv,uiv,1,ll+1,1);
            end if;
            ll := ll+ 1;
          end loop;
          iter := 0;
          m := m-1;
        when others => null;
      end case;
    end loop;
  end SVD;

  procedure SVD ( x : in out Standard_Complex_Matrices.Matrix;
                  n,p : in integer32;
                  s,e : out Standard_Complex_Vectors.Vector;
                  u : out Standard_Complex_Matrices.Matrix;
                  v : out Standard_Complex_Matrices.Matrix;
                  job : in integer32; info : out integer32 ) is

    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec;
    urv,uiv : Standard_Floating_VecVecs.Link_to_VecVec;
    vrv,viv : Standard_Floating_VecVecs.Link_to_VecVec;
    mm : constant integer32 := Min0(n+1,p);
    sr,si,er,ei,wr,wi : Standard_Floating_Vectors.Link_to_Vector;

  begin
    xrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    xiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    urv := Standard_Vector_Splitters.Allocate(n,n,1,1);
    uiv := Standard_Vector_Splitters.Allocate(n,n,1,1);
    vrv := Standard_Vector_Splitters.Allocate(p,p,1,1);
    viv := Standard_Vector_Splitters.Allocate(p,p,1,1);
    sr := new Standard_Floating_Vectors.Vector'(1..mm => 0.0);
    si := new Standard_Floating_Vectors.Vector'(1..mm => 0.0);
    er := new Standard_Floating_Vectors.Vector'(1..p => 0.0);
    ei := new Standard_Floating_Vectors.Vector'(1..p => 0.0);
    wr := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    wi := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    Standard_Matrix_Splitters.Complex_Parts(x,xrv,xiv);
    SVD(xrv,xiv,n,p,sr,si,er,ei,urv,uiv,vrv,viv,job,info,wr,wi);
    Standard_Vector_Splitters.Complex_Merge(sr,si,s);
    Standard_Vector_Splitters.Complex_Merge(er,ei,e);
    Standard_Matrix_Splitters.Complex_Merge(xrv,xiv,x);
    Standard_Matrix_Splitters.Complex_Merge(urv,uiv,u);
    Standard_Matrix_Splitters.Complex_Merge(vrv,viv,v);
    Standard_Floating_Vectors.Clear(sr);
    Standard_Floating_Vectors.Clear(si);
    Standard_Floating_Vectors.Clear(er);
    Standard_Floating_Vectors.Clear(ei);
    Standard_Floating_Vectors.Clear(wr);
    Standard_Floating_Vectors.Clear(wi);
    Standard_Floating_VecVecs.Deep_Clear(xrv);
    Standard_Floating_VecVecs.Deep_Clear(xiv);
    Standard_Floating_VecVecs.Deep_Clear(urv);
    Standard_Floating_VecVecs.Deep_Clear(uiv);
    Standard_Floating_VecVecs.Deep_Clear(vrv);
    Standard_Floating_VecVecs.Deep_Clear(viv);
  end SVD;

end Standard_Inlined_Singular_Values;
