with integer_io;                        use integer_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
--with Standard_Affine_Planes;            use Standard_Affine_Planes;
--with Standard_Affine_Solutions;         use Standard_Affine_Solutions;
with Standard_Plane_Operations;         use Standard_Plane_Operations;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Hypersurface_Samplers;             use Hypersurface_Samplers;
with P_Intrinsic_Diagonal_Continuation; use P_Intrinsic_Diagonal_Continuation;
with Affine_Sampling_Machine;           use Affine_Sampling_Machine;
--with Intrinsic_Sampling_Machine;        use Intrinsic_Sampling_Machine;
with Extrinsic_Diagonal_Homotopies;     use Extrinsic_Diagonal_Homotopies;

package body P_Intrinsic_Diagonal_Homotopies is

-- UTILITIES :

  function Product ( a,b : Vector ) return Solution_List is

    res,res_last : Solution_List;

  begin
    for i in a'range loop
      for j in b'range loop
        declare
          s : Solution(2);
        begin
          s.t := Create(0.0);
          s.m := 1;
          s.v(1) := a(i);
          s.v(2) := b(j);
          s.err := 0.0;
          s.rco := 1.0;
          s.res := 0.0;
          Append(res,res_last,s);
        end;
      end loop;
    end loop;
    return res;
  end Product;

  function Product ( sols : Solution_List; a : Vector ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      for i in a'range loop
        declare
          sa : Solution(ls.n+1);
        begin
          sa.t := Create(0.0);
          sa.m := 1;
          sa.v(ls.v'range) := ls.v;
          sa.v(sa.v'last) := a(i);
          sa.err := 0.0;
          sa.rco := 1.0;
          sa.res := 0.0;
          Append(res,res_last,sa);
        end;
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Product;

  function On_Hypersurface
             ( p : Poly; b,v : Vector; t : Complex_Number;
               tol : double_float ) return boolean is

    x : Vector(b'range) := b + t*v;
    y : Complex_Number := Eval(p,x);

  begin
    return (AbsVal(y) < tol);
  end On_Hypersurface;

  function On_Hypersurface
             ( p : Poly; b,v,t : Vector; tol : double_float )
             return Solution_List is

    res,res_last : Solution_List;
    common : Vector(1..t'length);
    cnt : natural := 0;

  begin
    for i in t'range loop
      if On_Hypersurface(p,b,v,t(i),tol)
       then cnt := cnt + 1;
            common(cnt) := t(i);   
      end if;
    end loop;
    if cnt > 0
     then declare
            s : Solution(cnt);
          begin
            s.t := Create(1.0);
            s.m := 1;
            s.v := common(1..cnt);
            s.err := 0.0;
            s.rco := 1.0;
            s.res := 0.0;
            Append(res,res_last,s);
          end;
    end if;
    return res;
  end On_Hypersurface;

  function Remove_on_Hypersurface
              ( p : Poly; b,v,t : Vector; tol : double_float )
              return Vector is

    res : Vector(t'range);
    ind : natural := t'first-1;

  begin
    for i in t'range loop
      if not On_Hypersurface(p,b,v,t(i),tol)
       then ind := ind + 1;
            res(ind) := t(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_on_Hypersurface;

-- WITNESS POINTS ON ONE HYPERSURFACE :

  procedure Hypersurface_Witness_Points
              ( n,dp : in natural; p : in Poly; ep : in Eval_Poly;
		b,v : in Vector; tp : out Vector; fail : out boolean ) is

    eps : constant double_float := 10.0**(-14 + dp);

  begin
    Generic_Points(p,ep,dp,b,v,eps,10*dp,tp,fail);
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( n,dp : in natural; p : in Poly; b,v : in Vector;
                tp : out Vector; fail : out boolean ) is

    ep : Eval_Poly := Create(p);

  begin
    Hypersurface_Witness_Points(n,dp,p,ep,b,v,tp,fail);
    Clear(ep);
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp : in natural; p : in Poly; ep : in Eval_Poly;
                b,v : in Vector; tp : out Vector; fail : out boolean ) is

    eps : constant double_float := 10.0**(-14 + dp);

  begin
    put(file,"The number of unknowns : "); put(file,n,1); new_line(file);
    put(file,"The degree of p : "); put(file,dp,1); new_line(file); 
    Generic_Points(file,p,ep,dp,b,v,eps,10*dp,tp,fail);
    if fail
     then put_line(file,"Calculation of witness points of p failed.");
     else put_line(file,"Calculation of witness points of p succeeded.");
    end if;
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp : in natural; p : in Poly; b,v : in Vector;
                tp : out Vector; fail : out boolean ) is

    ep : Eval_Poly := Create(p);

  begin
    Hypersurface_Witness_Points(file,n,dp,p,ep,b,v,tp,fail);
    Clear(ep);
  end Hypersurface_Witness_Points;

-- WITNESS POINTS ON TWO HYPERSURFACES :

  procedure Hypersurface_Witness_Points
              ( n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean ) is

    eps_p : constant double_float := 10.0**(-14 + dp);
    eps_q : constant double_float := 10.0**(-14 + dq);

  begin
    Generic_Points(p,ep,dp,bp,vp,eps_p,10*dp,tp,fail);
    if not fail 
     then Generic_Points(q,eq,dq,bq,vq,eps_q,10*dq,tq,fail);
    end if;
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( n,dp,dq : in natural; p,q : in Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean ) is

    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);

  begin
    Hypersurface_Witness_Points(n,dp,dq,p,q,ep,eq,bp,vp,bq,vq,tp,tq,fail);
    Clear(ep); Clear(eq);
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean ) is

    eps_p : constant double_float := 10.0**(-14 + dp);
    eps_q : constant double_float := 10.0**(-14 + dq);

  begin
    put(file,"The number of unknowns : "); put(file,n,1); new_line(file);
    put(file,"The degree of p : "); put(file,dp,1); new_line(file); 
    Generic_Points(file,p,ep,dp,bp,vp,eps_p,10*dp,tp,fail);
    if fail
     then put_line(file,"Calculation of witness points of p failed.");
     else put_line(file,"Calculation of witness points of p succeeded.");
    end if;
    put(file,"The degree of q : "); put(file,dq,1); new_line(file); 
    Generic_Points(file,q,eq,dq,bq,vq,eps_q,10*dq,tq,fail);
    if fail
     then put_line(file,"Calculation of witness points of q failed.");
     else put_line(file,"Calculation of witness points of q succeeded.");
    end if;
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp,dq : in natural; p,q : in Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean ) is

    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);

  begin
    Hypersurface_Witness_Points
      (file,n,dp,dq,p,q,ep,eq,bp,vp,bq,vq,tp,tq,fail);
    Clear(ep); Clear(eq);
  end Hypersurface_Witness_Points;

-- WITNESS POINTS ON ONE MULTIHOMOGENEOUS HYPERSURFACE :

  procedure Hypersurface_Witness_Points
              ( n,dp : in natural; p : in Poly; ep : in Eval_Poly;
                z : in Partition; b : in Vector; v : in VecVec;
                tp : out VecVec; fail : out boolean ) is

    eps : constant double_float := 10.0**(-14 + dp);
    maxit : constant natural := 10*dp;

  begin
    Generic_Points(p,ep,z,b,v,eps,maxit,tp,fail);
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( n : in natural; p : in Poly; z : in Partition;
                b : in Vector; v : in VecVec; tp : out VecVec;
                fail : out boolean ) is

    ep : Eval_Poly := Create(p);
    dp : constant natural := Degree(p);

  begin
    Hypersurface_Witness_Points(n,dp,p,ep,z,b,v,tp,fail);
    Clear(ep);
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp : in natural; p : in Poly; ep : in Eval_Poly;
                z : in Partition; b : in Vector; v : in VecVec; 
                tp : out VecVec; fail : out boolean ) is

    eps : constant double_float := 10.0**(-14 + dp);
    maxit : constant natural := 10*dp;

  begin
    put(file,"The number of unknowns : "); put(file,n,1); new_line(file);
    Generic_Points(file,p,ep,z,b,v,eps,maxit,tp,fail);
    if fail
     then put_line(file,"Calculation of witness points of p failed.");
     else put_line(file,"Calculation of witness points of p succeeded.");
    end if;
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n : in natural; p : in Poly; z : in Partition;
                b : in Vector; v : in VecVec; tp : out VecVec;
                fail : out boolean ) is

    ep : Eval_Poly := Create(p);
    dp : constant natural := Degree(p);

  begin
    Hypersurface_Witness_Points(file,n,dp,p,ep,z,b,v,tp,fail);
    Clear(ep);
  end Hypersurface_Witness_Points;

-- WITNESS POINTS ON TWO MULTIHOMOGENEOUS HYPERSURFACES :

  procedure Hypersurface_Witness_Points
              ( n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean ) is

    eps_p : constant double_float := 10.0**(-14 + dp);
    eps_q : constant double_float := 10.0**(-14 + dq);
    maxit_p : constant natural := 10*dp;
    maxit_q : constant natural := 10*dq;

  begin
    Generic_Points(p,ep,z,bp,vq,eps_p,maxit_p,tp,fail);
    if not fail
     then Generic_Points(q,eq,z,bq,vq,eps_q,maxit_q,tq,fail);
    end if;
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( n : in natural; p,q : in Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean ) is

    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);
    dp : constant natural := Degree(p);
    dq : constant natural := Degree(q);

  begin
    Hypersurface_Witness_Points
      (n,dp,dq,p,q,ep,eq,z,bp,bq,vp,vq,tp,tq,fail);
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean ) is

    eps_p : constant double_float := 10.0**(-14 + dp);
    eps_q : constant double_float := 10.0**(-14 + dq);
    maxit_p : constant natural := 10*dp;
    maxit_q : constant natural := 10*dq;

  begin
    put(file,"The number of unknowns : "); put(file,n,1); new_line(file);
    Generic_Points(file,p,ep,z,bp,vp,eps_p,maxit_p,tp,fail);
    if fail
     then put_line(file,"Calculation of witness points of p failed.");
     else put_line(file,"Calculation of witness points of p succeeded.");
          Generic_Points(file,q,eq,z,bq,vq,eps_q,maxit_q,tq,fail);
          if fail
           then put_line(file,"Calculation of witness points of q failed.");
           else put_line(file,"Calculation of witness points of q succeeded.");
          end if;
    end if;
  end Hypersurface_Witness_Points;

  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n : in natural; p,q : in Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean ) is

    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);
    dp : constant natural := Degree(p);
    dq : constant natural := Degree(q);

  begin
    Hypersurface_Witness_Points
      (file,n,dp,dq,p,q,ep,eq,z,bp,bq,vp,vq,tp,tq,fail);
  end Hypersurface_Witness_Points;

-- WITNESS POINTS ON THE INTERSECTION OF TWO HYPERSURFACES :

  function Special_Random ( v : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns a random vector of the same dimension and support as v.

    res : Vector(v'range) := Random_Vector(v'first,v'last);
    tol : constant double_float := 1.0E-10;

  begin
    for i in v'range loop
      if AbsVal(v(i)) < tol
       then res(i) := Create(0.0);
      end if;
    end loop;
    return res;
  end Special_Random;

  procedure Diagonal_Homotopy
              ( n : in natural; p,q : in Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec;
                sols : out Solution_List ) is

    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);

  begin
    Diagonal_Homotopy(n,p,q,ep,eq,bp,vp,bq,vq,tp,tq,b2,w2,sols);
    Clear(ep); Clear(eq);
  end Diagonal_Homotopy;

  procedure Diagonal_Homotopy
              ( n : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec;
                sols : out Solution_List ) is

   -- p2 : Poly := Append_Variables(n,p);
   -- q2 : Poly := Insert_Variables(n,q);
    n2 : constant natural := 2*n;
    w_1,w_2 : Vector(1..n2);
   -- p2q2 : Poly_Sys(1..2);
   -- ep2q2 : Eval_Poly_Sys(1..2);
   -- jm : Jaco_Mat(p2q2'range,1..n2);
   -- ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2));
    rv : Vector(1..n);
    start_b,target_b : Vector(1..n2);
    start_w,target_w : VecVec(1..2);
    tol : constant double_float := 1.0E-10;

  begin
    start_b(1..n) := bp; start_b(n+1..n2) := bq;
    target_b(1..n) := Random_Vector(1,n);
    target_b(n+1..n2) := target_b(1..n);
    w_1(1..n) := vp; w_1(n+1..n2) := (n+1..n2 => Create(0.0));
    w_2(1..n) := (1..n => Create(0.0)); w_2(n+1..n2) := vq;
    start_w(1) := new Standard_Complex_Vectors.Vector'(w_1);
    start_w(2) := new Standard_Complex_Vectors.Vector'(w_2);
   -- p2q2(1) := p2;          p2q2(2) := q2;
   -- ep2q2(1) := Create(p2); ep2q2(2) := Create(q2);
   -- jm := Create(p2q2);     ejm := Create(jm);
    target_w(1) := new Standard_Complex_Vectors.Vector(1..n2);
    rv := Special_Random(vp);
    for i in 1..n loop
      target_w(1)(i) := rv(i);
      target_w(1)(i+n) := rv(i);
    end loop;
    rv := Special_Random(vq);
    target_w(2) := new Standard_Complex_Vectors.Vector(1..n2);
    for i in 1..n loop
      target_w(2)(i) := rv(i);
      target_w(2)(i+n) := rv(i);
    end loop;
    sols := Product(tp,tq);
   -- version 0 :
   -- Silent_Affine_Sampler(ep2q2,ejm,start_b,target_b,start_w,target_w,sols);
   -- Silent_LU_Newton_Refiner
   --   (ep2q2,ejm,target_b,start_w,sols,1.0E-14,1.0E-14,4);
   -- version 1 :
   -- Silent_Diagonal_Continuation
   --   (n,p,q,start_b,target_b,start_w,target_w,sols);
   -- version 2 :
    declare
       ejm : Eval_Jaco_Mat(1..2,1..n) := Create(n,p,q);
    begin
      Silent_Diagonal_Continuation
        (n,ep,ep,ejm,start_b,target_b,start_w,target_w,sols);
      Clear(ejm);
    end;
    b2 := target_b; w2 := target_w;
  end Diagonal_Homotopy;

  procedure Diagonal_Homotopy
              ( file : in file_type;
                n : in natural; p,q : in Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec;
                sols : out Solution_List ) is

    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);

  begin
    Diagonal_Homotopy(file,n,p,q,ep,eq,bp,vp,bq,vq,tp,tq,b2,w2,sols);
    Clear(ep); Clear(eq);
  end Diagonal_Homotopy;

  procedure Diagonal_Homotopy
              ( file : in file_type;
                n : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec;
                sols : out Solution_List ) is

   -- p2 : Poly := Append_Variables(n,p);
   -- q2 : Poly := Insert_Variables(n,q);
    n2 : constant natural := 2*n;
    w_1,w_2 : Vector(1..n2);
   -- p2q2 : Poly_Sys(1..2);
   -- ep2q2 : Eval_Poly_Sys(1..2);
   -- jm : Jaco_Mat(p2q2'range,1..n2);
   -- ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2));
    rv : Vector(1..n);
    start_b,target_b : Vector(1..n2);
    start_w,target_w : VecVec(1..2);
    tol : constant double_float := 1.0E-10;
   -- x,y : Vector(1..2);

  begin
    start_b(1..n) := bp; start_b(n+1..n2) := bq;
    target_b(1..n) := Random_Vector(1,n);
    target_b(n+1..n2) := target_b(1..n);
    w_1(1..n) := vp; w_1(n+1..n2) := (n+1..n2 => Create(0.0));
    w_2(1..n) := (1..n => Create(0.0)); w_2(n+1..n2) := vq;
    start_w(1) := new Standard_Complex_Vectors.Vector'(w_1);
    start_w(2) := new Standard_Complex_Vectors.Vector'(w_2);
   -- put_line(file,"Evaluation at the product : ");
   -- for i in tp'range loop
   --   x(1) := tp(i);
   --   for j in tq'range loop
   --     x(2) := tq(j);
   --     y := Eval(ep2q2,start_b,start_w,x);
   --     put(file,"Evaluation at ("); put(file,i,1); put(file,",");
   --     put(file,j,1); put_line(file,") : "); put_line(file,y);
   --   end loop;
   -- end loop;
   -- p2q2(1) := p2;          p2q2(2) := q2;
   -- ep2q2(1) := Create(p2); ep2q2(2) := Create(q2);
   -- jm := Create(p2q2);     ejm := Create(jm);
    target_w(1) := new Standard_Complex_Vectors.Vector(1..n2);
    rv := Special_Random(vp);
    for i in 1..n loop
      target_w(1)(i) := rv(i);
      target_w(1)(i+n) := rv(i);
    end loop;
    rv := Special_Random(vq);
    target_w(2) := new Standard_Complex_Vectors.Vector(1..n2);
    for i in 1..n loop
      target_w(2)(i) := rv(i);
      target_w(2)(i+n) := rv(i);
    end loop;
    put_line(file,"Executing the intrinsic diagonal homotopy :");
    sols := Product(tp,tq);
   -- version 0 :
   -- Reporting_Affine_Sampler
   --   (file,ep2q2,ejm,start_b,target_b,start_w,target_w,sols);
   -- Reporting_LU_Newton_Refiner
   --   (file,ep2q2,ejm,target_b,target_w,sols,1.0E-14,1.0E-14,4);
   -- version 1 :
   -- Reporting_Diagonal_Continuation
   --   (file,n,p,q,start_b,target_b,start_w,target_w,sols);
   -- version 2 :
    declare
      ejm : Eval_Jaco_Mat(1..2,1..n) := Create(n,p,q);
    begin
      Reporting_Diagonal_Continuation
        (file,n,ep,eq,ejm,start_b,target_b,start_w,target_w,sols);
      Clear(ejm);
    end;
    w2 := target_w; b2 := target_b;
  end Diagonal_Homotopy;

  function Double_Start_Directions
              ( n : natural; v1 : VecVec; v2 : Vector ) return VecVec is

  -- DESCRIPTION :
  --   Doubles the directions of the affine slicing plane at the start.
  --   One direction is added, which is the direction of the line which
  --   defines the witness points of the new hypersurface we process 
  --   in the diagonal homotopy.  The direction of this line is in v2.

    n2 : constant natural := 2*n;
    res : VecVec(v1'first..v1'last+1);

  begin
    for i in v1'range loop     -- append zeros to given directions
      declare
        ww : Vector(1..n2);
      begin
        ww(v1(i)'range) := v1(i).all;
        ww(n+1..n2) := (n+1..n2 => Create(0.0));
        res(i) := new Vector'(ww);
      end;
    end loop;
    declare                   -- insert zeros for line cutting hypersurface
      ww : Vector(1..n2);
    begin
      ww(1..n) := (1..n => Create(0.0));
      ww(n+1..n2) := v2;
      res(res'last) := new Vector'(ww);
    end;
    return res;
  end Double_Start_Directions;

  function Double_Target_Directions
              ( n : natural; v1 : VecVec; v2 : Vector ) return VecVec is

  -- DESCRIPTION :
  --   Returns the directions of the target affine plane in the diagonal
  --   homotopy to intersection with a hypersurface.  The directions in v1
  --   are doubled and the doubled vector v2 occurs on the last position
  --   in the vector on return.

    n2 : constant natural := 2*n;
    res : VecVec(v1'first..v1'last+1);

  begin
    for i in v1'range loop     -- directions with same support as v1
      declare
        ww : Vector(1..n2);
        rv : Vector(1..n) := Special_Random(v1(i).all);
      begin
        ww(1..n) := rv;
        ww(n+1..n2) := rv;
        res(i) := new Vector'(ww);
      end;
    end loop;
    declare                    -- directions with same support as v2
      ww : Vector(1..n2);
      rv : Vector(1..n) := Special_Random(v2);
    begin
      ww(1..n) := rv;
      ww(n+1..n2) := rv;
      res(res'last) := new Vector'(ww);
    end;
    return res;
  end Double_Target_Directions;

  procedure Evaluate_Product
              ( file : in file_type; p : in Eval_Poly_Sys;
                b : in Vector; v : in VecVec; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the solutions (intrinsic format) in the system p.

    tmp : Solution_List := sols;

  begin
    put_line(file,"Evaluation of product ...");
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
        x : Vector(b'range) := b;
        y : Vector(p'range);
      begin
        for i in ls.v'range loop
          x := x + ls.v(i)*v(i).all;
        end loop;
        y := Eval(p,x);
        put_line(file,y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Evaluate_Product;

  procedure Diagonal_Homotopy
              ( n : in natural; p : in Poly_Sys; q : in Poly;
                bp : in Vector; vp : in VecVec; tp : in Solution_List;
                bq : in Vector; vq : in Vector; tq : in Vector;
                b2 : out Vector; w2 : out VecVec;
                sols : out Solution_List ) is

    n2 : constant natural := 2*n;
    p2 : Poly_Sys(p'range) := Append_Variables(n,p);
    q2 : Poly := Insert_Variables(n,q);
    p2q2 : Poly_Sys(p'first..p'last+1);
    ep2q2 : Eval_Poly_Sys(p2q2'range);
    jm : Jaco_Mat(p2q2'range,1..n2);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2));
    rv : Vector(1..n) := Random_Vector(1,n);
    start_b,target_b : Vector(1..n2);
    start_w : VecVec(vp'first..vp'last+1)
            := Double_Start_Directions(n,vp,vq);
    target_w : VecVec(vp'first..vp'last+1)
             := Double_Target_Directions(n,vp,vq);

  begin
    start_b(1..n) := bp;     target_b(1..n) := rv;
    start_b(n+1..n2) := bq;  target_b(n+1..n2) := rv;
    for i in p2'range loop
      p2q2(i) := p2(i);
    end loop;
    p2q2(p2q2'last) := q2;
    ep2q2 := Create(p2q2);
    jm := Create(p2q2);
    ejm := Create(jm);
    sols := Product(tp,tq);
   -- Evaluate_Product(Standard_Output,ep2q2,start_b,start_w,sols);
    Silent_Affine_Sampler
      (ep2q2,ejm,start_b,target_b,start_w,target_w,sols);
   -- Silent_LU_Newton_Refiner
    Silent_SV_Newton_Refiner
      (ep2q2,ejm,target_b,target_w,sols,1.0E-14,1.0E-14,4);
   -- Evaluate_Product(Standard_Output,ep2q2,target_b,target_w,sols);
    b2 := target_b; w2 := target_w;
  end Diagonal_Homotopy;

  procedure Diagonal_Homotopy
              ( file : in file_type;
                n : in natural; p : in Poly_Sys; q : in Poly;
                bp : in Vector; vp : in VecVec; tp : in Solution_List;
                bq : in Vector; vq : in Vector; tq : in Vector;
                b2 : out Vector; w2 : out VecVec;
                sols : out Solution_List ) is

    n2 : constant natural := 2*n;
    p2 : Poly_Sys(p'range) := Append_Variables(n,p);
    q2 : Poly := Insert_Variables(n,q);
    p2q2 : Poly_Sys(p'first..p'last+1);
    ep2q2 : Eval_Poly_Sys(p2q2'range);
    jm : Jaco_Mat(p2q2'range,1..n2);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2));
    rv : Vector(1..n) := Random_Vector(1,n);
    start_b,target_b : Vector(1..n2);
    start_w : VecVec(vp'first..vp'last+1)
            := Double_Start_Directions(n,vp,vq);
    target_w : VecVec(vp'first..vp'last+1)
             := Double_Target_Directions(n,vp,vq);

  begin
    start_b(1..n) := bp;     target_b(1..n) := rv;
    start_b(n+1..n2) := bq;  target_b(n+1..n2) := rv;
    for i in p2'range loop
      p2q2(i) := p2(i);
    end loop;
    p2q2(p2q2'last) := q2;
    ep2q2 := Create(p2q2);
    jm := Create(p2q2);
    ejm := Create(jm);
    sols := Product(tp,tq);
   -- Evaluate_Product(file,ep2q2,start_b,start_w,sols);
    Reporting_Affine_Sampler
      (file,ep2q2,ejm,start_b,target_b,start_w,target_w,sols);
   -- Reporting_LU_Newton_Refiner
    Reporting_SV_Newton_Refiner
      (file,ep2q2,ejm,target_b,target_w,sols,1.0E-14,1.0E-14,4);
   -- Evaluate_Product(file,ep2q2,target_b,target_w,sols);
    b2 := target_b; w2 := target_w;
  end Diagonal_Homotopy;

  procedure Combine_Solutions
              ( file : in file_type; n : in natural; sols : in Solution_List;
                b2 : in Vector; w2 : in VecVec ) is

    trunc_sols : Solution_List := Truncate(sols,n);
    w : VecVec(w2'range) := Truncate(w2,n);
   -- combi_sols : Solution_List := Combine(trunc_sols,b2(1..n),w);
    combi_sols : Solution_List := Expand(trunc_sols,b2(1..n),w);

  begin
    put_line(file,"The solutions in extrinsic format : ");
    put(file,Length_Of(combi_sols),Head_Of(combi_sols).n,combi_sols);
  end Combine_Solutions;

  procedure Combine_Solutions
              ( file : in file_type; n : in natural;
                sols : in Array_of_Solution_Lists;
                b2 : in Vector; w2 : in Array_of_VecVecs ) is

    w : Array_of_VecVecs(w2'range) := Truncate(w2,n);

  begin
    for i in sols'range loop
      put(file,"Solution list "); put(file,i,1);
      if Is_Null(sols(i))
       then put_line(file," is empty.");
       else put_line(file," in extrinsic format : ");
            declare
              trunc_sols : Solution_List := Truncate(sols(i),n);
              combi_sols : Solution_List
                        -- := Combine(trunc_sols,b2(1..n),w(i).all);
                         := Expand(trunc_sols,b2(1..n),w(i).all);
            begin
              put(file,Length_Of(combi_sols),Head_Of(combi_sols).n,combi_sols);
            end;
      end if;
    end loop;
  end Combine_Solutions;

end P_Intrinsic_Diagonal_Homotopies;
