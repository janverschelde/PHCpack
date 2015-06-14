with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Matrices;

package body Multprec_Complex_Linear_Solvers is

-- AUXLILIARIES :

  function dconjg ( x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the complex conjugate of x.

    res : Complex_Number;
    re : Floating_Number := REAL_PART(x);
    im : Floating_Number := IMAG_PART(x);

  begin
    Min(im);
    res := Create(re,im);
    Clear(re); Clear(im);
    return res;
  end dconjg;

  function csign ( x,y : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns |x|*y/|y|.

    res : Complex_Number;
    fltacc : Floating_Number;
    cmpacc : Complex_Number;

  begin
    fltacc := AbsVal(x);
    res := Create(fltacc);
    Clear(fltacc);
    Mul(res,y);
    fltacc := AbsVal(y);
    cmpacc := Create(fltacc);
    Div(res,cmpacc);
    Clear(fltacc);
    Clear(cmpacc);
    return res;
  end csign;

  function Inverse_Abs_Sum ( z : Multprec_Complex_Vectors.Vector )
                           return Floating_Number is

  -- DESCRIPTION :
  --   Returns the reciprocal of the sum of the absolute values in z.

    res,sum,acc : Floating_Number;

  begin
    sum := Create(integer(0));
    for i in z'range loop
      acc := AbsVal(z(i));
      Add(sum,acc);
      Clear(acc);
    end loop;
    res := Create(integer(1));
    Div(res,sum);
    Clear(sum);
    return res;
  end Inverse_Abs_Sum;

-- TARGET ROUTINES :

  procedure Scale ( a : in out Multprec_Complex_Matrices.Matrix;
                    b : in out Multprec_Complex_Vectors.Vector ) is

    fac : Complex_Number;

    function Maximum ( a : in Multprec_Complex_Matrices.Matrix;
                       i : in integer32 )
                     return Complex_Number is

      res : integer32 := a'first(2);
      max : Floating_Number := AbsVal(a(i,res));
      tmp : Floating_Number;

    begin
      for j in a'first(2)+1..a'last(2) loop
        tmp := AbsVal(a(i,j));
        if tmp > max
         then max := tmp; res := j;
        end if;
      end loop;
      return a(i,res);
    end Maximum;

    procedure Divide ( a : in out Multprec_Complex_Matrices.Matrix;
                       b : in out Multprec_Complex_Vectors.Vector;
                       i : in integer32; fac : in Complex_Number ) is
    begin
      for j in a'range(2) loop
        a(i,j) := a(i,j)/fac;
      end loop;
      b(i) := b(i)/fac;
    end Divide;
  
  begin
    for i in a'range(1) loop
      fac := Maximum(a,i);
      Divide(a,b,i,fac);
    end loop;
  end Scale;

  function Norm1 ( a : Multprec_Complex_Matrices.Matrix )
                 return Floating_Number is

    res : Floating_Number := Create(0.0);
    sum,acc : Floating_Number;

  begin
    for j in a'range(2) loop
      sum := Create(0.0);
      for i in a'range(1) loop
        acc := AbsVal(a(i,j));
        Add(sum,acc); Clear(acc);
      end loop;
      if sum > res
       then Copy(sum,res);
      end if; 
      Clear(sum);
    end loop;
    return res;
  end Norm1;

  function Norm1 ( a : Multprec_Complex_VecVecs.VecVec )
                 return Floating_Number is

    res : Floating_Number := Create(0.0);
    sum,acc : Floating_Number;
    aj : Multprec_Complex_Vectors.Link_to_Vector;

  begin
    for j in a'range loop
      sum := Create(0.0);
      aj := a(j);
      for i in aj'range loop
        acc := AbsVal(aj(i));
        Add(sum,acc); Clear(acc);
      end loop;
      if sum > res
       then Copy(sum,res);
      end if; 
      Clear(sum);
    end loop;
    return res;
  end Norm1;

  procedure lufac ( a : in out Multprec_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

    kp1,ell,nm1 : integer32;
    smax,fltacc : Floating_Number;
    temp,cmpacc : Complex_Number;

  begin
    info := 0;
    nm1 := n - 1;
    if nm1 >= 1 then
      for k in 1..nm1 loop
        kp1 := k + 1;
        ell := k;
        smax := AbsVal(a(k,k));                  -- find the pivot index ell
        for i in kp1..n loop
          fltacc := AbsVal(a(i,k));
          if fltacc > smax then
            ell := i;
            Copy(fltacc,smax);
          end if;
          Clear(fltacc);
        end loop;
        ipvt(k) := ell;
        if Equal(smax,0.0) then     -- this column is already triangularized
          info := k;
        else
          if ell /= k then                       -- interchange if necessary
            Copy(a(ell,k),temp);
            Copy(a(k,k),a(ell,k));
            Copy(temp,a(k,k)); Clear(temp);
          end if;
          cmpacc := Create(integer(1));               -- compute multipliers
          Div(cmpacc,a(k,k));
          Min(cmpacc);
          for i in kp1..n loop
            Mul(a(i,k),cmpacc);
          end loop;
          Clear(cmpacc);
          for j in kp1..n loop                           -- row elimination
            Copy(a(ell,j),temp);
            if ell /= k then
              Copy(a(k,j),a(ell,j));
              Copy(temp,a(k,j));
            end if;
            for i in kp1..n loop
              cmpacc := temp*a(i,k);
              Add(a(i,j),cmpacc);
              Clear(cmpacc);
            end loop;
          end loop;
          Clear(temp);
        end if;
        Clear(smax);
      end loop;
    end if;
    ipvt(n) := n;
    fltacc := AbsVal(a(n,n));
    if Equal(fltacc,0.0)
     then info := n;
    end if;
    Clear(fltacc);
  end lufac;

  procedure lufac ( a : in out Multprec_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

    kp1,ell,nm1 : integer32;
    ak,aj : Multprec_Complex_Vectors.Link_to_Vector;
    smax,fltacc : Floating_Number;
    temp,cmpacc : Complex_Number;

  begin
    info := 0;
    nm1 := n - 1;
    if nm1 >= 1 then
      for k in 1..nm1 loop
        kp1 := k + 1;
        ell := k;
        ak := a(k);
        smax := AbsVal(ak(k));                  -- find the pivot index ell
        for i in kp1..n loop
          fltacc := AbsVal(ak(i));
          if fltacc > smax then
            ell := i;
            Copy(fltacc,smax);
          end if;
          Clear(fltacc);
        end loop;
        ipvt(k) := ell;
        if Equal(smax,0.0) then     -- this column is already triangularized
          info := k;
        else
          if ell /= k then                       -- interchange if necessary
            Copy(ak(ell),temp);
            Copy(ak(k),ak(ell));
            Copy(temp,ak(k)); Clear(temp);
          end if;
          cmpacc := Create(integer(1));               -- compute multipliers
          Div(cmpacc,ak(k));
          Min(cmpacc);
          for i in kp1..n loop
            Mul(ak(i),cmpacc);
          end loop;
          Clear(cmpacc);
          for j in kp1..n loop                           -- row elimination
            aj := a(j);
            Copy(aj(ell),temp);
            if ell /= k then
              Copy(aj(k),aj(ell));
              Copy(temp,aj(k));
            end if;
            for i in kp1..n loop
              cmpacc := temp*ak(i);
              Add(aj(i),cmpacc);
              Clear(cmpacc);
            end loop;
          end loop;
          Clear(temp);
        end if;
        Clear(smax);
      end loop;
    end if;
    ipvt(n) := n;
    fltacc := AbsVal(a(n)(n));
    if Equal(fltacc,0.0)
     then info := n;
    end if;
    Clear(fltacc);
  end lufac;

  procedure estco ( a : in Multprec_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in Floating_Number;
                    rcond : out Floating_Number ) is

  -- NOTE :
  --   rcond = 1/(norm(a)*(estimate of norm(inverse(a))))
  --   estimate = norm(z)/norm(y) where a*z = y and ctrans(a)*y = e.
  --   ctrans(a) is the conjugate transpose of a.
  --   The components of e are chosen to cause maximum local
  --   growth in the elements of w where ctrans(u)*w = e.
  --   The vectors are frequently rescaled to avoid overflow.

    z : Multprec_Complex_Vectors.Vector(1..n);
    kb,kp1,ell : integer32;
    s,sm,sum,ynorm,fltacc1,fltacc2 : Floating_Number;
    ek,t,wk,wkm,cmpacc1,cmpacc2 : Complex_Number;

    use Multprec_Complex_Vectors;

  begin
    ek := Create(integer(1));                        -- solve ctrans(u)*w = e
    for j in 1..n loop
      z(j) := Create(integer(0));
    end loop;
    for k in 1..n loop
      fltacc1 := AbsVal(z(k));
      if not Equal(fltacc1,0.0) then
        cmpacc1 := -z(k);
        Copy(ek,cmpacc2);
        Clear(ek);
        ek := csign(cmpacc2,cmpacc1);
        Clear(cmpacc1);
        Clear(cmpacc2);
      end if;
      Clear(fltacc1);
      cmpacc1 := ek - z(k);
      fltacc1 := AbsVal(cmpacc1);
      Clear(cmpacc1);
      fltacc2 := AbsVal(a(k,k));
      if fltacc1 > fltacc2 then
        s := fltacc2/fltacc1;
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Mul(ek,cmpacc1);
        Clear(cmpacc1);
        Clear(s);
      end if;
      Clear(fltacc1); Clear(fltacc2);
      wk := ek - z(k);
      wkm := ek + z(k);
      Min(wkm);
      s := AbsVal(wk);
      sm := AbsVal(wkm);
      fltacc1 := AbsVal(a(k,k));
      if Equal(fltacc1,0.0) then
        Clear(wk);  wk := Create(integer(1));
        Clear(wkm); wkm := Create(integer(1));
      else
        cmpacc1 := dconjg(a(k,k));
        Div(wk,cmpacc1);
        Div(wkm,cmpacc1);
        Clear(cmpacc1);
      end if;
      Clear(fltacc1);
      kp1 := k + 1;
      if kp1 <= n then
        for j in kp1..n loop
          cmpacc2 := dconjg(a(k,j));
          cmpacc1 := wkm*cmpacc2;
          Add(cmpacc1,z(j));
          fltacc1 := AbsVal(cmpacc1);
          Add(sm,fltacc1);
          Clear(fltacc1); 
          Clear(cmpacc1);
          cmpacc1 := wk*cmpacc2;
          Add(z(j),cmpacc1);
          Clear(cmpacc1);
          Clear(cmpacc2);
          fltacc1 := AbsVal(z(j));
          Add(s,fltacc1);
          Clear(fltacc1);
        end loop;
        if s < sm then
          t := wkm - wk;
          Copy(wkm,wk);
          for j in kp1..n loop
            cmpacc2 := dconjg(a(k,j));
            cmpacc1 := t*cmpacc2;
            Add(z(j),cmpacc1);
            Clear(cmpacc1);
            Clear(cmpacc2);
          end loop;
          Clear(t);
        end if;
      end if;
      Copy(wk,z(k));
      Clear(wk); Clear(wkm);
      Clear(s); Clear(sm);
    end loop;
    Clear(ek);
    s := Inverse_Abs_Sum(z);
    cmpacc1 := Create(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    Clear(s);
    for k in 1..n loop                           -- solve ctrans(l)*y = w
      kb := n+1-k;
      if kb < n then
        t := Create(integer(0));
        for i in (kb+1)..n loop
          cmpacc2 := dconjg(a(i,kb));
          cmpacc1 := cmpacc2*z(i);
          Add(t,cmpacc1);
          Clear(cmpacc1);
          Clear(cmpacc2);
        end loop;
        Add(z(kb),t);
        Clear(t);
      end if;
      fltacc1 := AbsVal(z(kb));
      if fltacc1 > 1.0 then
        s := Create(integer(1));
        Div(s,fltacc1);
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Clear(cmpacc1);
        Clear(s);
      end if;
      Clear(fltacc1);
      ell := ipvt(kb);
      if ell /= kb then
        Copy(z(ell),t);
        Copy(z(kb),z(ell));
        Copy(t,z(kb));
        Clear(t);
      end if;
    end loop;
    s := Inverse_Abs_Sum(z);
    cmpacc1 := Create(s);
    Clear(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    ynorm := Create(integer(1));
    for k in 1..n loop                                    -- solve l*v = y
      ell := ipvt(k);
      if ell /= k then
        Copy(z(ell),t);
        Copy(z(k),z(ell));
        Copy(t,z(k));
      else
        Copy(z(ell),t);
      end if;
      if k < n then
        for i in (k+1)..n loop
          cmpacc1 := t*a(i,k);
          Add(z(i),cmpacc1);
          Clear(cmpacc1);
        end loop;
      end if;
      Clear(t);
      fltacc1 := AbsVal(z(k));
      if fltacc1 > 1.0 then
        s := Create(integer(1));
        Div(s,fltacc1);
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Clear(cmpacc1);
        Mul(ynorm,s);
        Clear(s);
      end if;
      Clear(fltacc1);
    end loop;
    s := Inverse_Abs_Sum(z);
    cmpacc1 := Create(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    Mul(ynorm,s);
    Clear(s);
    for k in 1..n loop                                    -- solve u*z = v
      kb := n+1-k;
      fltacc1 := AbsVal(z(kb));
      fltacc2 := AbsVal(a(kb,kb));
      if fltacc1 > fltacc2 then
        s := fltacc2/fltacc1;
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Clear(cmpacc1);
        Mul(ynorm,s);
        Clear(s);
      end if;
      if Equal(fltacc2,0.0) then
        Clear(z(kb));
        z(kb) := Create(integer(1));
      else
        Div(z(kb),a(kb,kb));
      end if;
      Clear(fltacc1);
      Clear(fltacc2);
      t := -z(kb);
      for i in 1..(kb-1) loop
        cmpacc1 := t*a(i,kb);
        Add(z(i),cmpacc1);
        Clear(cmpacc1);
      end loop;
      Clear(t);
    end loop;
    s := Inverse_Abs_Sum(z);                             -- make znorm = 1.0
    cmpacc1 := Create(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    Mul(ynorm,s);
    Clear(s);
    if Equal(anorm,0.0)
     then rcond := Create(integer(0));
     else rcond := ynorm/anorm;
    end if;
    Clear(ynorm); Clear(z);
  end estco;

  procedure estco ( a : in Multprec_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in Floating_Number;
                    rcond : out Floating_Number ) is

  -- NOTE :
  --   rcond = 1/(norm(a)*(estimate of norm(inverse(a))))
  --   estimate = norm(z)/norm(y) where a*z = y and ctrans(a)*y = e.
  --   ctrans(a) is the conjugate transpose of a.
  --   The components of e are chosen to cause maximum local
  --   growth in the elements of w where ctrans(u)*w = e.
  --   The vectors are frequently rescaled to avoid overflow.

    z : Multprec_Complex_Vectors.Vector(1..n);
    ak,aj : Multprec_Complex_Vectors.Link_to_Vector;
    kb,kp1,ell : integer32;
    s,sm,sum,ynorm,fltacc1,fltacc2 : Floating_Number;
    ek,t,wk,wkm,cmpacc1,cmpacc2 : Complex_Number;

    use Multprec_Complex_Vectors;

  begin
    ek := Create(integer(1));                        -- solve ctrans(u)*w = e
    for j in 1..n loop
      z(j) := Create(integer(0));
    end loop;
    for k in 1..n loop
      ak := a(k);
      fltacc1 := AbsVal(z(k));
      if not Equal(fltacc1,0.0) then
        cmpacc1 := -z(k);
        Copy(ek,cmpacc2);
        Clear(ek);
        ek := csign(cmpacc2,cmpacc1);
        Clear(cmpacc1);
        Clear(cmpacc2);
      end if;
      Clear(fltacc1);
      cmpacc1 := ek - z(k);
      fltacc1 := AbsVal(cmpacc1);
      Clear(cmpacc1);
      fltacc2 := AbsVal(ak(k));
      if fltacc1 > fltacc2 then
        s := fltacc2/fltacc1;
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Mul(ek,cmpacc1);
        Clear(cmpacc1);
        Clear(s);
      end if;
      Clear(fltacc1); Clear(fltacc2);
      wk := ek - z(k);
      wkm := ek + z(k);
      Min(wkm);
      s := AbsVal(wk);
      sm := AbsVal(wkm);
      fltacc1 := AbsVal(ak(k));
      if Equal(fltacc1,0.0) then
        Clear(wk);  wk := Create(integer(1));
        Clear(wkm); wkm := Create(integer(1));
      else
        cmpacc1 := dconjg(ak(k));
        Div(wk,cmpacc1);
        Div(wkm,cmpacc1);
        Clear(cmpacc1);
      end if;
      Clear(fltacc1);
      kp1 := k + 1;
      if kp1 <= n then
        for j in kp1..n loop
          aj := a(j);
          cmpacc2 := dconjg(aj(k));
          cmpacc1 := wkm*cmpacc2;
          Add(cmpacc1,z(j));
          fltacc1 := AbsVal(cmpacc1);
          Add(sm,fltacc1);
          Clear(fltacc1); 
          Clear(cmpacc1);
          cmpacc1 := wk*cmpacc2;
          Add(z(j),cmpacc1);
          Clear(cmpacc1);
          Clear(cmpacc2);
          fltacc1 := AbsVal(z(j));
          Add(s,fltacc1);
          Clear(fltacc1);
        end loop;
        if s < sm then
          t := wkm - wk;
          Copy(wkm,wk);
          for j in kp1..n loop
            aj := a(j);
            cmpacc2 := dconjg(aj(k));
            cmpacc1 := t*cmpacc2;
            Add(z(j),cmpacc1);
            Clear(cmpacc1);
            Clear(cmpacc2);
          end loop;
          Clear(t);
        end if;
      end if;
      Copy(wk,z(k));
      Clear(wk); Clear(wkm);
      Clear(s); Clear(sm);
    end loop;
    Clear(ek);
    s := Inverse_Abs_Sum(z);
    cmpacc1 := Create(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    Clear(s);
    for k in 1..n loop                           -- solve ctrans(l)*y = w
      kb := n+1-k;
      if kb < n then
        t := Create(integer(0));
        ak := a(kb);
        for i in (kb+1)..n loop
          cmpacc2 := dconjg(ak(i));
          cmpacc1 := cmpacc2*z(i);
          Add(t,cmpacc1);
          Clear(cmpacc1);
          Clear(cmpacc2);
        end loop;
        Add(z(kb),t);
        Clear(t);
      end if;
      fltacc1 := AbsVal(z(kb));
      if fltacc1 > 1.0 then
        s := Create(integer(1));
        Div(s,fltacc1);
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Clear(cmpacc1);
        Clear(s);
      end if;
      Clear(fltacc1);
      ell := ipvt(kb);
      if ell /= kb then
        Copy(z(ell),t);
        Copy(z(kb),z(ell));
        Copy(t,z(kb));
        Clear(t);
      end if;
    end loop;
    s := Inverse_Abs_Sum(z);
    cmpacc1 := Create(s);
    Clear(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    ynorm := Create(integer(1));
    for k in 1..n loop                                    -- solve l*v = y
      ell := ipvt(k);
      if ell /= k then
        Copy(z(ell),t);
        Copy(z(k),z(ell));
        Copy(t,z(k));
      else
        Copy(z(ell),t);
      end if;
      if k < n then
        ak := a(k);
        for i in (k+1)..n loop
          cmpacc1 := t*ak(i);
          Add(z(i),cmpacc1);
          Clear(cmpacc1);
        end loop;
      end if;
      Clear(t);
      fltacc1 := AbsVal(z(k));
      if fltacc1 > 1.0 then
        s := Create(integer(1));
        Div(s,fltacc1);
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Clear(cmpacc1);
        Mul(ynorm,s);
        Clear(s);
      end if;
      Clear(fltacc1);
    end loop;
    s := Inverse_Abs_Sum(z);
    cmpacc1 := Create(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    Mul(ynorm,s);
    Clear(s);
    for k in 1..n loop                                    -- solve u*z = v
      kb := n+1-k;
      ak := a(kb);
      fltacc1 := AbsVal(z(kb));
      fltacc2 := AbsVal(ak(kb));
      if fltacc1 > fltacc2 then
        s := fltacc2/fltacc1;
        cmpacc1 := Create(s);
        Mul(z,cmpacc1);
        Clear(cmpacc1);
        Mul(ynorm,s);
        Clear(s);
      end if;
      if Equal(fltacc2,0.0) then
        Clear(z(kb));
        z(kb) := Create(integer(1));
      else
        Div(z(kb),ak(kb));
      end if;
      Clear(fltacc1);
      Clear(fltacc2);
      t := -z(kb);
      for i in 1..(kb-1) loop
        cmpacc1 := t*ak(i);
        Add(z(i),cmpacc1);
        Clear(cmpacc1);
      end loop;
      Clear(t);
    end loop;
    s := Inverse_Abs_Sum(z);                             -- make znorm = 1.0
    cmpacc1 := Create(s);
    Mul(z,cmpacc1);
    Clear(cmpacc1);
    Mul(ynorm,s);
    Clear(s);
    if Equal(anorm,0.0)
     then rcond := Create(integer(0));
     else rcond := ynorm/anorm;
    end if;
    Clear(ynorm); Clear(z);
  end estco;

  procedure lufco ( a : in out Multprec_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out Floating_Number ) is

    anorm : Floating_Number := Norm1(a);
    info : integer32;

  begin
    lufac(a,n,ipvt,info);
    if info = 0
     then estco(a,n,ipvt,anorm,rcond);
     else rcond := Create(0.0);
    end if;
    Clear(anorm);
  end lufco;

  procedure lufco ( a : in out Multprec_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out Floating_Number ) is

    anorm : Floating_Number := Norm1(a);
    info : integer32;

  begin
    lufac(a,n,ipvt,info);
    if info = 0
     then estco(a,n,ipvt,anorm,rcond);
     else rcond := Create(0.0);
    end if;
    Clear(anorm);
  end lufco;

  procedure lusolve ( a : in Multprec_Complex_Matrices.Matrix;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Multprec_Complex_Vectors.Vector ) is

    ell,nm1,kb : integer32;
    temp,acc : Complex_Number;
 
  begin
    nm1 := n-1;
    if nm1 >= 1 then                                        -- solve L*y = b
      for k in 1..nm1 loop
        ell := ipvt(k);
        Copy(b(ell),temp);
        if ell /= k then
          Copy(b(k),b(ell));
          Copy(temp,b(k));
        end if;
        for i in (k+1)..n loop
          acc := temp*a(i,k);
          Add(b(i),acc);
          Clear(acc);
        end loop;
        Clear(temp);
      end loop;
    end if;
    for k in 1..n loop                                     -- solve U*x = y
      kb := n+1-k;
      Div(b(kb),a(kb,kb));
      temp := -b(kb);
      for j in 1..(kb-1) loop
        acc := temp*a(j,kb);
        Add(b(j),acc);
        Clear(acc);
      end loop;
      Clear(temp);
    end loop;
  end lusolve;

  procedure lusolve ( a : in Multprec_Complex_VecVecs.VecVec;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Multprec_Complex_Vectors.Vector ) is

    ell,nm1,kb : integer32;
    ak : Multprec_Complex_Vectors.Link_to_Vector;
    temp,acc : Complex_Number;
 
  begin
    nm1 := n-1;
    if nm1 >= 1 then                                        -- solve L*y = b
      for k in 1..nm1 loop
        ell := ipvt(k);
        Copy(b(ell),temp);
        if ell /= k then
          Copy(b(k),b(ell));
          Copy(temp,b(k));
        end if;
        ak := a(k);
        for i in (k+1)..n loop
          acc := temp*ak(i);
          Add(b(i),acc);
          Clear(acc);
        end loop;
        Clear(temp);
      end loop;
    end if;
    for k in 1..n loop                                     -- solve U*x = y
      kb := n+1-k;
      ak := a(kb);
      Div(b(kb),ak(kb));
      temp := -b(kb);
      for j in 1..(kb-1) loop
        acc := temp*ak(j);
        Add(b(j),acc);
        Clear(acc);
      end loop;
      Clear(temp);
    end loop;
  end lusolve;

  procedure Examine ( file : in file_type; tol : in double_float;
                      mpc : in Multprec_Complex_Numbers.Complex_Number;
                      stc : in Standard_Complex_Numbers.Complex_Number ) is

    rnd : Standard_Complex_Numbers.Complex_Number := Round(mpc);
    absdif : double_float;

    use Standard_Complex_Numbers;

  begin
    rnd := rnd - stc;
    absdif := AbsVal(rnd);
    if absdif > tol then
      put(file,"INCONSISTENTIE : diff = ");
      put(file,absdif); new_line(file);
    end if;
  end Examine;

  procedure Triangulate ( a : in out Multprec_Complex_Matrices.Matrix;
                          tol : in double_float;
                          size : in natural32; n,m : in integer32 ) is

    max,cbs : Floating_Number;
    temp : Complex_Number;
    pivot,k,kcolumn : integer32;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := Create(integer(0));                              -- find pivot
      pivot := 0;
      for l in k..n loop
        cbs := AbsVal(a(l,kcolumn));
        if ((cbs > tol) and then (cbs > max))
         then Copy(cbs,max); pivot := l;
        end if;
        Clear(cbs);
      end loop;
      Clear(max);
      if pivot = 0 then
        kcolumn := kcolumn + 1;
      else
        if pivot /= k then                      -- interchange if necessary
          for i in 1..m loop
            Copy(a(pivot,i),temp);
            Copy(a(k,i),a(pivot,i));
            Copy(temp,a(k,i));
            Clear(temp);
          end loop;
        end if;
        for j in (kcolumn+1)..m loop                      -- triangulate a
          Div(a(k,j),a(k,kcolumn));
        end loop;
        Clear(a(k,kcolumn));
        a(k,kcolumn) := Create(integer(1));
        Set_Size(a(k,kcolumn),size);
        for i in (k+1)..n loop
          for j in (kcolumn+1)..m loop
            temp := a(i,kcolumn)*a(k,j);
            Sub(a(i,j),temp);
            Clear(temp);
          end loop;
          Clear(a(i,kcolumn));
          a(i,kcolumn) := Create(integer(0));
          Set_Size(a(i,kcolumn),size);
        end loop;
        k := k + 1;
        kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Triangulate;

  procedure Triangulate ( file : in file_type;
                          a : in out Multprec_Complex_Matrices.Matrix;
                          tol : in double_float;
                          size : in natural32; n,m : in integer32 ) is

    max,cbs : Floating_Number;
    temp : Complex_Number;
    pivot,k,kcolumn : integer32;
    shadow : Standard_Complex_Matrices.Matrix(1..n,1..m);
    statmp : Standard_Complex_Numbers.Complex_Number;

    use Standard_Complex_Numbers;

  begin
    put_line(file,"Starting triangulate");
    put(file,"#rows : "); put(file,n,1);
    put(file,"  #columns : "); put(file,m,1); new_line(file);
    for i in 1..n loop
      for j in 1..m loop
        shadow(i,j) := Round(a(i,j));
      end loop;
    end loop;
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := Create(integer(0));                               -- find pivot
      pivot := 0;
      for l in k..n loop
        cbs := AbsVal(a(l,kcolumn));
        if ((cbs > tol) and then (cbs > max)) then
          Copy(cbs,max);
          pivot := l;
        end if;
        Clear(cbs);
      end loop;
      Clear(max);
      put(file,"pivot index : "); put(file,pivot,1); new_line(file);
      if pivot /= 0 then
        put(file,"pivot element : ");
        put(file,a(pivot,kcolumn)); new_line(file);
        put(file,"its shadow : ");
        put(file,shadow(pivot,kcolumn)); new_line(file);
        Examine(file,tol,a(pivot,kcolumn),shadow(pivot,kcolumn));
      end if;
      if pivot = 0 then
        kcolumn := kcolumn + 1;
      else
        if pivot /= k then                       -- interchange if necessary
          for i in 1..m loop
            Copy(a(pivot,i),temp);
            Copy(a(k,i),a(pivot,i));
            Copy(temp,a(k,i));
            Clear(temp);
            statmp := shadow(pivot,i);
            shadow(pivot,i) := shadow(k,i);
            shadow(k,i) := statmp;
          end loop;
        end if;
        put(file,"Dividing row "); put(file,k,1); put_line(file," :");
        for j in (kcolumn+1)..m loop                        -- triangulate a
          put(file,"before : "); put(file,a(k,j));
          new_line(file);
          put(file,"its shadow : "); put(file,shadow(k,j));
          new_line(file);
          Examine(file,tol,a(k,j),shadow(k,j));
          put(file,"the divisor : "); put(file,a(k,kcolumn));
          new_line(file);
          put(file,"its shadow : "); put(file,shadow(k,kcolumn));
          new_line(file);
          Examine(file,tol,a(k,kcolumn),shadow(k,kcolumn));
          Div(a(k,j),a(k,kcolumn));
          shadow(k,j) := shadow(k,j)/shadow(k,kcolumn);
          put(file,"after : "); put(file,a(k,j)); new_line(file);
          put(file,"its shadow element : "); put(file,shadow(k,j));
          new_line(file);
          Examine(file,tol,a(k,j),shadow(k,j));
        end loop;
        Clear(a(k,kcolumn));
        a(k,kcolumn) := Create(integer(1));
        Set_Size(a(k,kcolumn),size);
        shadow(k,kcolumn) := Create(1.0);
        for i in (k+1)..n loop
          for j in (kcolumn+1)..m loop
            temp := a(i,kcolumn)*a(k,j);
            Sub(a(i,j),temp);
            Clear(temp);
            shadow(i,j) := shadow(i,j) - shadow(i,kcolumn)*shadow(k,j);
          end loop;
          Clear(a(i,kcolumn));
          a(i,kcolumn) := Create(integer(0));
          Set_Size(a(i,kcolumn),size);
          shadow(i,kcolumn) := Create(0.0);
        end loop;
        put(file,"Row "); put(file,k+1,1);
        put_line(file," after update : ");
        for i in 1..m loop
          put(file,a(k+1,i)); new_line(file);
          put(file,"its shadow element : ");
          put(file,shadow(k+1,i)); new_line(file);
          Examine(file,tol,a(k+1,i),shadow(k+1,i));
        end loop;
        k := k + 1;
        kcolumn := kcolumn + 1;
      end if;
    end loop;
    put_line(file,"leaving triangulate");
  end Triangulate;

  procedure Diagonalize ( a : in out Multprec_Complex_Matrices.Matrix;
                          n,m : in integer32 ) is

    max : Floating_Number;
    temp : Complex_Number;
    pivot,k,kcolumn : integer32;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := Create(integer(0));                               -- find pivot
      for l in k..n loop
        if AbsVal(a(l,kcolumn)) > max then
          max := AbsVal(a(l,kcolumn));
          pivot := l;
        end if;
      end loop;
      if Equal(max,0.0) then
        kcolumn := kcolumn + 1;
      else
        if pivot /= k then                       -- interchange if necessary
          for i in 1..m loop
            temp := a(pivot,i);
            a(pivot,i) := a(k,i);
            a(k,i) := temp;
          end loop;
        end if;
        for j in (kcolumn+1)..m loop                        -- diagonalize a
          a(k,j) := a(k,j) / a(k,kcolumn);
        end loop;
        a(k,kcolumn) := Create(integer(1));
        for i in 1..(k-1) loop
          for j in (kcolumn+1)..m loop
            a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
          end loop;
        end loop;
        for i in (k+1)..n loop
          for j in (kcolumn+1)..m loop
            a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
          end loop;
        end loop;
        for j in 1..(k-1) loop
          a(j,kcolumn) := Create(integer(0));
        end loop;
        for j in (k+1)..n loop
          a(j,kcolumn) := Create(integer(0));
        end loop;
        k := k + 1;
        kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Diagonalize;

end Multprec_Complex_Linear_Solvers;
