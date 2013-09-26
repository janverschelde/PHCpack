with integer_io;                        use integer_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;         use Standard_Solution_Filters;
with Hypersurface_Samplers;             use Hypersurface_Samplers;
--with Standard_Affine_Planes;            use Standard_Affine_Planes;
--with Standard_Affine_Solutions;         use Standard_Affine_Solutions;
with Standard_Plane_Operations;         use Standard_Plane_Operations;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Witness_Sets_Formats;              use Witness_Sets_Formats;
with Extrinsic_Diagonal_Homotopies;     use Extrinsic_Diagonal_Homotopies;
with P_Intrinsic_Diagonal_Homotopies;   use P_Intrinsic_Diagonal_Homotopies;
with Multihomogeneous_Solutions;        use Multihomogeneous_Solutions;

package body P_Intrinsic_Diagonal_Solvers is

  procedure Witness_Points
              ( n : in natural; p,q : in Poly; b2 : out Vector;
                w2 : out VecVec; sols : out Array_of_Solution_Lists;
                fail : out boolean ) is

    dp : constant natural := Degree(p);
    dq : constant natural := Degree(q);
    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);
    bp : Vector(1..n) := Random_Vector(1,n);
    vp : Vector(1..n) := Random_Vector(1,n);
    bq : Vector(1..n) := Random_Vector(1,n);
    vq : Vector(1..n) := Random_Vector(1,n);
    tp : Vector(1..dp);
    tq : Vector(1..dq);
    tol : constant double_float := 1.0E-8;
    nbcom : natural;
    paths,end_paths : Solution_List;

  begin
   -- Hypersurface_Witness_Points
   --   (Standard_Output,n,dp,dq,p,q,bp,vp,bq,vq,tp,tq,fail);
    Hypersurface_Witness_Points
      (Standard_Output,n,dp,dq,p,q,ep,eq,bp,vp,bq,vq,tp,tq,fail);
    if fail
     then put_line("Computation of witness points failed.");
     else put_line("Witness points of p : "); put_line(tp);
          put_line("Witness points of q : "); put_line(tq);
          sols(1) :=  On_Hypersurface(q,bp,vp,tp,tol);
          nbcom := Length_Of(sols(1));
          put("Number of common factors : "); put(nbcom,1); new_line;
          if nbcom = 0
           then Diagonal_Homotopy
                  (Standard_Output,n,p,q,bp,vp,bq,vq,tp,tq,b2,w2,paths);
           else declare
                  rtp : constant Vector
                      := Remove_on_Hypersurface(q,bp,vp,tp,tol);
                  rtq : constant Vector
                      := Remove_on_Hypersurface(p,bq,vq,tq,tol);
                begin
                  Diagonal_Homotopy
                    (Standard_Output,n,p,q,bp,vp,bq,vq,rtp,rtq,b2,w2,paths);
                end;
          end if;
          end_paths := On_Target_Filter(paths,Create(1.0),tol);
          sols(2) := Vanishing_Filter(end_paths,tol);
    end if;
  end Witness_Points;

  procedure Witness_Points 
              ( file : in file_type;
                n : in natural; p,q : in Poly; b2 : out Vector;
                w2 : out VecVec; sols : out Array_of_Solution_Lists;
                fail : out boolean ) is

    dp : constant natural := Degree(p);
    dq : constant natural := Degree(q);
    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);
    bp : Vector(1..n) := Random_Vector(1,n);
    vp : Vector(1..n) := Random_Vector(1,n);
    bq : Vector(1..n) := Random_Vector(1,n);
    vq : Vector(1..n) := Random_Vector(1,n);
    tp : Vector(1..dp);
    tq : Vector(1..dq);
    tol : constant double_float := 1.0E-8;
    nbcom : natural;
    paths,end_paths : Solution_List;

  begin
   -- Hypersurface_Witness_Points(file,n,dp,dq,p,q,bp,vp,bq,vq,tp,tq,fail);
    Hypersurface_Witness_Points
      (file,n,dp,dq,p,q,ep,eq,bp,vp,bq,vq,tp,tq,fail);
    if fail
     then put_line(file,"Computation of witness points failed.");
     else put_line(file,"Witness points of p : "); put_line(file,tp);
          put_line(file,"Witness points of q : "); put_line(file,tq);
          sols(1) := On_Hypersurface(q,bp,vp,tp,tol);
          nbcom := Length_Of(sols(1));
          put(file,"Number of common factors : ");
          put(file,nbcom,1); new_line(file);
          if nbcom = 0
           then Diagonal_Homotopy(file,n,p,q,bp,vp,bq,vq,tp,tq,b2,w2,paths);
           else declare
                  rtp : constant Vector
                      := Remove_on_Hypersurface(q,bp,vp,tp,tol);
                  rtq : constant Vector
                      := Remove_on_Hypersurface(q,bq,vq,tq,tol);
                begin
                  Diagonal_Homotopy
                    (file,n,p,q,bp,vp,bq,vq,rtp,rtq,b2,w2,paths);
                end;
          end if;
          end_paths := On_Target_Filter(paths,Create(1.0),tol);
          sols(2) := Vanishing_Filter(end_paths,tol);
    end if;
  end Witness_Points;

  procedure Write_Witness_Points
              ( file : in file_type; v,t : in VecVec ) is

  -- DESCRIPTION :
  --   Writes the directions v(i) and values for t which define
  --   witness points on the affine lines b + t*v(i).

  begin
    for i in v'range loop
      put_line(file,"special direction :");
      put_line(file,v(i).all);
      if t(i) = null
       then put_line(file,"no witness points for this direction.");
       else put_line(file,"the values for t are :");
            put_line(file,t(i).all);
      end if;
    end loop;
  end Write_Witness_Points;

  function Admissible ( v1,v2 : Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the special directions v1 and v2 will admit
  --   solutions on the diagonal for a multi-homogeneous structure.

    cnt : natural := 0;

  begin
    for i in v1'range loop
      if v1(i) /= Create(0.0)
       then cnt := cnt + 1;
       elsif v2(i) /= Create(0.0)
           then cnt := cnt + 1;
      end if;
    end loop;
    return (cnt > 1);
  end Admissible;

  function Admissible ( v1 : VecVec; v2 : Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the special directions v1 and v2 will admit 
  --   solutions on the diagonal for a multi-homogeneous structure.

  -- REQUIRED : the vectors in v1 span an orthogonal basis.

    tol : constant double_float := 1.0E-8;

  begin
    return not In_Span(v1,v2,tol);
  end Admissible;

  procedure Witness_Points 
              ( n : in natural; p,q : in Poly; z : in Partition;
                b2 : out Vector; w2 : out Array_of_VecVecs;
                sols : out Array_of_Solution_Lists ) is

    bp : Vector(1..n) := Random_Vector(1,n);
    bq : Vector(1..n) := Random_Vector(1,n);
    vp : VecVec(z'range) := Random_Multihomogeneous_Directions(n,z);
    vq : VecVec(z'range) := Random_Multihomogeneous_Directions(n,z);
    tp,tq : VecVec(z'range);
    fail : boolean;
    ind : natural;

  begin
    Hypersurface_Witness_Points(n,p,q,z,bp,bq,vp,vq,tp,tq,fail);
    if not fail
     then ind := 0;
          for i in vp'range loop
            for j in vq'range loop
              ind := ind + 1;
              if Admissible(vp(i).all,vq(j).all)
               then if tp(i) /= null and tq(j) /= null
                     then declare
                            ww : VecVec(1..2);
                          begin
                            Diagonal_Homotopy
                              (n,p,q,bp,vp(i).all,bq,vq(j).all,
                               tp(i).all,tq(j).all,b2,ww,sols(ind));
                            w2(ind) := new VecVec'(ww);
                          end;
                    end if;
              end if;
            end loop;
          end loop;
    end if;
  end Witness_Points;

  procedure Witness_Points 
              ( file : in file_type;
                n : in natural; p,q : in Poly; z : in Partition;
                b2 : out Vector; w2 : out Array_of_VecVecs;
                sols : out Array_of_Solution_Lists ) is

    bp : Vector(1..n) := Random_Vector(1,n);
    bq : Vector(1..n) := Random_Vector(1,n);
    vp : VecVec(z'range) := Random_Multihomogeneous_Directions(n,z);
    vq : VecVec(z'range) := Random_Multihomogeneous_Directions(n,z);
    tp,tq : VecVec(z'range);
    fail : boolean;
    ind : natural;

  begin
    Hypersurface_Witness_Points(file,n,p,q,z,bp,bq,vp,vq,tp,tq,fail);
    if fail
     then put_line(file,"Computation of witness points failed.");
     else put_line(file,"Witness points of p : ");
          Write_Witness_Points(file,vp,tp);
          put_line(file,"Witness points of q : ");
          Write_Witness_Points(file,vq,tq);
          ind := 0;
          for i in vp'range loop
            for j in vq'range loop
              ind := ind + 1;
              put(file,"Combining direction "); put(file,i,1);
              put(file," with direction "); put(file,j,1);
              if not Admissible(vp(i).all,vq(j).all)
               then put_line(file," will not lead to a solution.");
               else if tp(i) /= null and tq(j) /= null
                     then put_line(file,"...");
                          declare
                            ww : VecVec(1..2);
                          begin
                            Diagonal_Homotopy
                              (file,n,p,q,bp,vp(i).all,bq,vq(j).all,
                               tp(i).all,tq(j).all,b2,ww,sols(ind));
                            w2(ind) := new VecVec'(ww);
                          end;
                     else put_line(file," has no start solutions.");
                    end if;
              end if;
            end loop;
          end loop;
    end if;
  end Witness_Points;

  function Remove_Common_Factors
              ( file : file_type; p : Poly_Sys;
                b,v,t : Vector; tol : double_float ) return Vector is

  -- DESCRIPTION :
  --   Removes points tp(k) from tp for which b + t(k)*v satisfies
  --   one of the polynomials in p.

    res : Vector(t'range);
    cnt : natural := t'first-1;
    retain : boolean;

  begin
    for i in t'range loop
      retain := true;
      declare
        x : constant Vector := b + t(i)*v;
        y : Complex_Number;
      begin
        for j in p'range loop
          y := Eval(p(j),x);
          put(file,"  p("); put(file,j,1);
          put(file,") = "); put(file,y); new_line(file);
          if AbsVal(y) < tol
           then retain := false;
          end if;
          exit when not retain;
        end loop;
      end;
      if retain
       then put(file,"witness point "); put(file,i,1);
            put(file," does not satisfy any previous equation.");
            put_line(file,"  Retained");
            cnt := cnt + 1;
            res(cnt) := t(i);
       else put(file,"witness point "); put(file,i,1);
            put(file," satisfies some previous equation.");
	    put_line(file,"  Discarded");
      end if;
    end loop;
    return res(res'first..cnt);
  end Remove_Common_Factors;

  procedure Filter_by_Evaluation
              ( file : in file_type; n : in natural; p : in Poly;
                b : in Vector; v : in VecVec; tol : in double_float;
		sols : in Solution_List; on_p,off_p : out Solution_List ) is

    on_p_last,off_p_last : Solution_List;
    newsols : Solution_List := Expand(sols,b,v); -- Combine(sols,b,v);
    tmp : Solution_List := newsols;
    restmp : Solution_List := sols;
    ep : Eval_Poly := Create(p);
    cnt : natural := 0;

  begin
    for i in 1..Length_Of(newsols) loop
      put(file,"Value at witness point "); put(file,i,1); put(file," : ");
      declare
        ls : Link_to_Solution := Head_Of(tmp);
        resls : Link_to_Solution := Head_Of(restmp);
        eva : Complex_Number := Eval(ep,ls.v(1..n));
        abseva : double_float := AbsVal(eva);
      begin
        put(file,abseva); new_line(file);
        if abseva > tol
         then Append(off_p,off_p_last,resls.all);
         else cnt := cnt + 1;
              Append(on_p,on_p_last,resls.all);
        end if;
      end;
      tmp := Tail_Of(tmp);
      restmp := Tail_Of(restmp);
    end loop;
    Clear(ep);
    put(file,"#solutions on hypersurface : "); put(file,cnt,1);
    new_line(file);
  end Filter_by_Evaluation;

  procedure Shift_to_End ( p : in out Poly_Sys; ind : in natural ) is

    backup : Poly := p(ind);

  begin
    for i in ind+1..p'last loop
      p(i-1) := p(i);
    end loop;
    p(p'last) := backup;
  end Shift_to_End;

  procedure Sanity_Check ( file : in file_type; p : in Poly_Sys;
                           n,k,ind : in natural; b2 : in Vector;
                           w : in VecVec; sols : in Solution_List;
                           tol : in double_float ) is

  -- DESCRIPTION :
  --   Checks whether the shifted polynomials are still superfluous.

    on_p,off_p : Solution_List;

  begin
    if ind = p'last
     then put_line(file,
                   "Checking if shifted polynomials are superfluous...");
          for i in ind..p'last loop
            Filter_by_Evaluation
              (file,n,p(i),b2,w(1..k),tol,sols,on_p,off_p);
            if Is_Null(off_p)
             then put(file,"Polynomial "); put(file,i,1);
                  put_line(file," is still superfluous.");
             else put(file,"Polynomial "); put(file,i,1);
                  put_line(file," is no longer superfluous!");
            end if;
            Clear(on_p); Clear(off_p);
          end loop;
    end if;
  end Sanity_Check;

  procedure Witness_Points
              ( file : in file_type; n : in natural;
                p : in out Poly_Sys; k,ind : in out natural;
                b2 : in out Vector; w : in out VecVec;
		sols : in out Solution_List; fail : out boolean ) is

    p3 : Poly := p(ind);
    dp : constant natural := Degree(p3);
    bp : Vector(1..n) := Random_Vector(1,n);
    vp : Vector(1..n) := Random_Vector(1,n);
    tp : Vector(1..dp);
    tol_vanish : constant double_float := 1.0E-8;
    tol_singular : constant double_float := 1.0E-8;
    w3 : VecVec(1..k+1);
    sols3,tar_sols,van_sols,reg_sols,onp3,offp3 : Solution_List;

  begin
    Filter_by_Evaluation(file,n,p3,b2,w(1..k),tol_vanish,sols,onp3,offp3);
    if Is_Null(offp3)
     then put(file,"Equation "); put(file,ind,1);
          put_line(file," is superfluous, will skip it.");
          Shift_to_End(p,ind);
     else if not Is_Null(onp3)
           then Write_Embedding(file,p,n,k,b2(1..n),w(1..k),onp3);
          end if;
          Copy(offp3,sols); Clear(offp3);
          Hypersurface_Witness_Points(file,n,dp,p3,bp,vp,tp,fail);
          if fail
           then put(file,"Failed to compute witness points at equation ");
                put(file,ind,1); new_line(file);
           else declare
                  ntp : constant Vector
                      := Remove_Common_Factors
                           (file,p(p'first..p'first+ind-2),
                            bp,vp,tp,tol_vanish);
                begin
                  Diagonal_Homotopy
                    (file,n,p(p'first..p'first+ind-2),p3,
                          b2(1..n),w(1..k),sols,bp,vp,ntp,b2,w3,sols3);
                end;
                tar_sols := On_Target_Filter(sols3,Create(1.0),tol_vanish);
                van_sols := Vanishing_Filter(tar_sols,tol_vanish);
                reg_sols := Regular_Filter(van_sols,tol_singular);
                put_line(file,"The regular vanishing solutions : ");
                put(file,Length_Of(reg_sols),Head_Of(reg_sols).n,reg_sols);
                Copy(reg_sols,sols);
                Clear(sols3);
                Clear(tar_sols); Clear(van_sols); Clear(reg_sols);
                k := k+1;
                for i in 1..k loop
                  Clear(w(i));
                  w(i) := new Vector'(w3(i).all);
                  Clear(w3(i));
                end loop;
                ind := ind+1;
          end if;
    end if;
  end Witness_Points;

  procedure Total_Degree_Hypersurface_Solver
              ( file : in file_type; n : in natural; p : in out Poly_Sys;
                b : out Vector; w : out VecVec;
                sols : out Array_of_Solution_lists) is

    n2 : constant natural := 2*n;
    b2 : Vector(1..n2);
    work : Solution_List;
    fail : boolean;
    k,ind : natural;
    tol : constant double_float := 1.0E-8;

  begin
    Witness_Points(file,n,p(p'first),p(p'first+1),b2,w(1..2),sols,fail);
    if not Is_Null(sols(1))
     then Write_Embedding(file,p,n,b2(1..n),w(1).all,sols(1));
    end if;
    work := sols(2);
    if not fail and (p'length > 2)
     then k := 2;
          ind := p'first+2;
          for i in p'first+2..p'last loop
            Witness_Points(file,n,p,k,ind,b2,w,work,fail);
          end loop;
          Sanity_Check(file,p,n,k,ind,b2,w,work,tol);
    end if;
    Combine_Solutions(file,n,work,b2,w);
    b := b2(1..n);
  end Total_Degree_Hypersurface_Solver;

  procedure Multihomogeneous_Intersection
              ( file : in file_type;
                mhs : in out Link_to_Multihomogeneous_Solution;
                p : in Poly_Sys; z : in Partition;
                q_v,q_t : in VecVec ) is

  -- DESCRIPTION :
  --   Intersection of a k-dimensional multi-homogeneous solution component
  --   with a multi-homogeneous hypersurface q(x) = 0.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   mhs      multi-homogeneous solution set of co-dimension k = mhs.k;
  --   p        the first k equations define the solution set in mhs;
  --   z        partition of the set of variables, defines the type
  --            of multi-homogenization;
  --   q_v      directions of the affine special lines cutting p(k+1);
  --   q_t      values for the the witness points on the special lines.

  -- ON RETURN :
  --   mhs      multi-homogeneous solution set of co-dimension k-1.

    m : constant natural := z'length;
    nbm : constant natural := mhs.nb*m;
    sols : Array_of_Solution_Lists(1..nbm);
    n2 : constant natural :=  2*mhs.n;
    b2 : Vector(1..n2);
    w2 : Array_of_VecVecs(1..nbm);
    ind : natural := 0;
    new_mhs : Link_to_Multihomogeneous_Solution;

  begin
    b2(1..mhs.n) := mhs.b;
    b2(mhs.n+1..b2'last) := mhs.b;
    for i in mhs.w'range loop
      for j in q_v'range loop
        ind := ind + 1;
        put(file,"Vectors "); put(file,i,1);
        put(file," with vector "); put(file,j,1);
        if q_t(j) = null
         then
           put_line(file," has no start solution.");
         else
           if not Admissible(mhs.w(i).all,q_v(j).all)
            then put_line(file," will not lead to any solution.");
            else put_line(file," invoking a diagonal homotopy...");
                 declare
                   ww : VecVec(1..mhs.k+1);
                 begin
                   Diagonal_Homotopy
                     (file,mhs.n,p(1..mhs.k),p(mhs.k+1),
                           mhs.b,mhs.w(i).all,mhs.s(i),
                           mhs.b,q_v(j).all,q_t(j).all,b2,ww,sols(ind));
                   w2(ind) := new VecVec'(ww);
                 end;
           end if;
        end if;
      end loop;
    end loop;
    new_mhs := Create(mhs.n,mhs.b,w2,sols);
    Clear(mhs);
    mhs := new_mhs;
  end Multihomogeneous_Intersection;

  procedure Multihomogeneous_Hypersurface_Solver
              ( file : in file_type; n : in natural; z : in Partition;
                p : in out Poly_Sys ) is

    n2 : constant natural := 2*n;
    b2 : Vector(1..n2);
    m : constant natural := z'length;
    m2 : constant natural := m*m;
    w2 : Array_of_VecVecs(1..m2);
    sols : Array_of_Solution_Lists(1..m2);
    mhs : Link_to_Multihomogeneous_Solution;

  begin
    Witness_Points(file,n,p(p'first),p(p'first+1),z,b2,w2,sols);
    mhs := Create(n,b2(1..n),w2,sols);
    for i in 3..p'last loop
      put(file,"Intersection with polynomial ");
      put(file,i,1); put(file," ...");
      declare
        pi_v : VecVec(z'range) := Random_Multihomogeneous_Directions(n,z);
        pi_t : VecVec(1..m);
        fail : boolean;
      begin
        Hypersurface_Witness_Points
          (file,n,p(i),z,b2(1..n),pi_v,pi_t,fail);
        Multihomogeneous_Intersection(file,mhs,p,z,pi_v,pi_t);
      end;
    end loop;
    Combine_Solutions(file,mhs.n,mhs.s,mhs.b,mhs.w);
  end Multihomogeneous_Hypersurface_Solver;

end P_Intrinsic_Diagonal_Solvers;
