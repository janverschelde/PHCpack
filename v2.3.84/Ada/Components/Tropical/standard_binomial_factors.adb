with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Power_Transformations;
with Standard_Durand_Kerner;             use Standard_Durand_Kerner;
with Standard_Lattice_Polygons;          use Standard_Lattice_Polygons;
with Standard_Initial_Forms;             use Standard_Initial_Forms;

package body Standard_Binomial_Factors is

-- I. compute tropisms

  function Common_Normals
               ( N,M : Standard_Integer64_Matrices.Matrix ) return List is

    res,res_last : List;
    v : Standard_Integer_Vectors.Vector(1..2);

  begin
    for i in N'range(2) loop
      for j in M'range(2) loop
        if (N(1,i) = M(1,j)) and (N(2,i) = M(2,j)) then
          v(1) := integer32(N(1,i)); v(2) := integer32(N(2,i));
          Append(res,res_last,v);
        end if;
      end loop;
    end loop;
    return res;
  end Common_Normals;

  function Common_Inner_Normals
               ( A,B : Standard_Integer64_Matrices.Matrix ) return List is

  -- DESCRIPTION :
  --   Returns the common inner normals for the points with coordinates
  --   in the columns of A and B.

    V : constant Standard_Integer64_Matrices.Matrix := Convex_Hull_2D(A);
    N : constant Standard_Integer64_Matrices.Matrix := Inner_Normals(V);
    W : constant Standard_Integer64_Matrices.Matrix := Convex_Hull_2D(B);
    M : constant Standard_Integer64_Matrices.Matrix := Inner_Normals(W);

  begin
    return Common_Normals(N,M);
  end Common_Inner_Normals;

  procedure Common_Inner_Normals
               ( A,B : in Standard_Integer64_Matrices.Matrix;
                 output : in boolean;
                 V,W,N,M,T : out List; fail : out boolean ) is

    Vm : constant Standard_Integer64_Matrices.Matrix := Convex_Hull_2D(A);
    Nm : constant Standard_Integer64_Matrices.Matrix := Inner_Normals(Vm);
    Wm : constant Standard_Integer64_Matrices.Matrix := Convex_Hull_2D(B);
    Mm : constant Standard_Integer64_Matrices.Matrix := Inner_Normals(Wm);
    bug1,bug2 : boolean;

  begin
    Check(A,Vm,Nm,false,bug1);
    if output then
      if bug1 then
        Check(A,Vm,Nm,true,bug1);
      else
        put_line("check of first computation of inner normals passed");
        put_line("vertices of the first Newton polygon :"); put(Vm);
        put_line("inner normals to edges of the Newton polygon :"); put(Nm);
      end if;
    end if;
    Convert(Vm,V); Convert(Nm,N);
    Check(B,Wm,Mm,false,bug2);
    if output then
      if bug2 then
        Check(B,Wm,Mm,true,bug2);
      else
        put_line("check of second computation of inner normals passed");
        put_line("vertices of the Newton polygon of :"); put(Wm);
        put_line("inner normals to edges of the Newton polygon :"); put(Mm);
      end if;
    end if;
    Convert(Wm,W); Convert(Mm,M);
    fail := bug1 or bug2;
    if not bug1 and not bug2 then
      T := Common_Normals(Nm,Mm);
      if output then
        if Is_Null(T)
         then put_line("no common inner normals");
         else put_line("the candidate tropisms :"); put(T);
        end if;
      end if;
    end if;
  end Common_Inner_Normals;

  function Common_Inner_Normals ( f,g : Poly ) return List is

    p : List := Create(f);
    q : List := Create(g);
    A : Standard_Integer64_Matrices.Matrix(1..2,1..integer32(Length_Of(p)));
    B : Standard_Integer64_Matrices.Matrix(1..2,1..integer32(Length_Of(q)));
    
  begin
    Convert(p,A); Lexicographic_Decreasing_Sort(A); Clear(p);
    Convert(q,B); Lexicographic_Decreasing_Sort(B); Clear(q);
    return Common_Inner_Normals(A,B);
  end Common_Inner_Normals;

  procedure Common_Inner_Normals
              ( f,g : in Poly; output : in boolean;
                A,B,V,W,N,M,T : out List; fail : out boolean ) is

    p : List := Create(f);
    q : List := Create(g);
    Ap : Standard_Integer64_Matrices.Matrix(1..2,1..integer32(Length_Of(p)));
    Bq : Standard_Integer64_Matrices.Matrix(1..2,1..integer32(Length_Of(q)));

  begin
    Convert(p,Ap); Lexicographic_Decreasing_Sort(Ap); Clear(p);
    Convert(q,Bq); Lexicographic_Decreasing_Sort(Bq); Clear(q);
    Convert(Ap,A); Convert(Bq,B);
    Common_Inner_Normals(Ap,Bq,output,V,W,N,M,T,fail);
  end Common_Inner_Normals;

-- II. compute initial roots

  function Coefficients ( p : Poly ) return Standard_Complex_Vectors.Vector is

    min : constant integer32 := Minimal_Degree(p,1);
    max : constant integer32 := Maximal_Degree(p,1);
    res : Standard_Complex_Vectors.Vector(min..max);

    procedure Coeff_Term ( t : in Term; continue : out boolean ) is
    begin
      res(t.dg(1)) := t.cf;
      continue := true;
    end Coeff_Term;
    procedure Coeff_Terms is new Visiting_Iterator(Coeff_Term);

  begin
    for i in res'range loop
      res(i) := Create(0.0);
    end loop;
    Coeff_Terms(p);
    return res;
  end Coefficients;

  function Normalize ( c : Standard_Complex_Vectors.Vector )
                     return Standard_Complex_Vectors.Vector is

    d : constant integer32 := c'first;

  begin
    if d = 0 then
      return c;
    else
      declare
        res : Standard_Complex_Vectors.Vector(0..c'last-d);
      begin
        for i in c'range loop
          res(i-d) := c(i);
        end loop;
        return res;
      end;
    end if;
  end Normalize;

  function Roots ( c : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector is

    tol : constant double_float := 1.0E-12;
    max : constant natural32 := 20;

  begin
    return Roots(c,max,tol);
  end Roots;

  function Roots ( c : Standard_Complex_Vectors.Vector;
                   max : natural32; tol : double_float )
                 return Standard_Complex_Vectors.Vector is
 
    res : Standard_Complex_Vectors.Vector(1..c'last)
        := Standard_Random_Vectors.Random_Vector(1,c'last);
    r : Standard_Complex_Vectors.Vector(res'range);
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,res,r,max,tol,nb,fail);
    return res;
  end Roots;

  procedure Roots ( c : in Standard_Complex_Vectors.Vector;
                    max : in natural32; tol : in double_float;
                    z : out Standard_Complex_Vectors.Vector;
                    nb : out natural32; res : out double_float;
                    fail : out boolean ) is
 
    r : Standard_Complex_Vectors.Vector(z'range);

  begin
    z := Standard_Random_Vectors.Random_Vector(1,c'last);
    Silent_Durand_Kerner(c,z,r,max,tol,nb,fail);
    res := Max_Norm(r);
  end Roots;

  function Match ( v : Standard_Integer_Vectors.Vector;
                   z1,z2 : Standard_Complex_Vectors.Vector;
                   tol : double_float ) return List_of_Terms is

    res,res_last : List_of_Terms;

  begin
    for i in z1'range loop
      for j in z2'range loop
        if AbsVal(z1(i) - z2(j)) < tol then
          declare
            t : Term;
          begin      
            t.cf := (z1(i) + z2(j))/2.0;
            t.dg := new Standard_Integer_Vectors.Vector'(v);
            Append(res,res_last,t);
          end;
        end if;
      end loop;
    end loop;
    return res;
  end Match;

  function Common_Initial_Roots
              ( f,g : Poly; v : Standard_Integer_Vectors.Vector;
                tol : double_float ) return List_of_Terms is

    c1 : constant Standard_Complex_Vectors.Vector := Coefficients(f);
    n1 : constant Standard_Complex_Vectors.Vector := Normalize(c1);
    c2 : constant Standard_Complex_Vectors.Vector := Coefficients(g);
    n2 : constant Standard_Complex_Vectors.Vector := Normalize(c2);
    z1 : constant Standard_Complex_Vectors.Vector := Roots(n1);
    z2 : constant Standard_Complex_Vectors.Vector := Roots(n2);

  begin
    return Match(v,z1,z2,tol);
  end Common_Initial_Roots;

  procedure Common_Initial_Roots
              ( f,g : in Poly; v : in Standard_Integer_Vectors.Vector;
                output : in boolean;
                max : in natural32; eps,tol : in double_float;
                t : out List_of_Terms; fail : out boolean ) is

    c1 : constant Standard_Complex_Vectors.Vector := Coefficients(f);
    n1 : constant Standard_Complex_Vectors.Vector := Normalize(c1);
    c2 : constant Standard_Complex_Vectors.Vector := Coefficients(g);
    n2 : constant Standard_Complex_Vectors.Vector := Normalize(c2);
    z1 : Standard_Complex_Vectors.Vector(1..n1'last);
    z2 : Standard_Complex_Vectors.Vector(1..n2'last);
    nb1,nb2 : natural32;
    res1,res2 : double_float;
    fail1,fail2 : boolean;

  begin
    Roots(n1,max,eps,z1,nb1,res1,fail1);
    Roots(n2,max,eps,z2,nb2,res2,fail2);
    fail := fail1 or fail2;
    t := Match(v,z1,z2,tol);
    if output then
      put_line("Roots of the first polynomial :"); put_line(z1);
      put("#steps : "); put(nb1,1);
      put("  Residual : "); put(res1,3);
      if fail1
       then put_line("  failure");
       else put_line("  success");
      end if;
      put_line("Roots of the second polynomial :"); put_line(z2);
      put("#steps : "); put(nb2,1);
      put("  Residual : "); put(res2,3);
      if fail2
       then put_line("  failure");
       else put_line("  success");
      end if;
    end if;
  end Common_Initial_Roots;

  function Initial_Terms
              ( f,g : Poly; v : Standard_Integer_Vectors.Vector )
              return List_of_Terms is

    res : List_of_Terms;
    tol : constant double_float := 1.0E-8;
    f1 : Poly := Initial(f,v);
    g1 : Poly := Initial(g,v);
    k : constant integer32 := Standard_Power_Transformations.Pivot(v,v'first);
    m : constant Standard_Integer_Matrices.Matrix(1..2,1..2)
      := Standard_Power_Transformations.Eliminate(v,k);
    tf1 : Poly := Transform(f1,m);
    tg1 : Poly := Transform(g1,m);
    ef1 : Poly := Eliminate(tf1,k);
    eg1 : Poly := Eliminate(tg1,k);

  begin
    res := Common_Initial_Roots(ef1,eg1,v,tol);
    Clear(f1); Clear(g1); Clear(tf1); Clear(tg1); Clear(ef1); Clear(eg1);
    return res;
  end Initial_Terms;

  procedure Initial_Terms
              ( f,g : in Poly; v : in Standard_Integer_Vectors.Vector;
                output : in boolean;
                max : in natural32; eps,tol : in double_float;
                t : out List_of_Terms; fail : out boolean ) is

    f1 : Poly := Initial(f,v);
    g1 : Poly := Initial(g,v);
    k : constant integer32 := Standard_Power_Transformations.Pivot(v,v'first);
    m : constant Standard_Integer_Matrices.Matrix(1..2,1..2)
      := Standard_Power_Transformations.Eliminate(v,k);
    tf1 : Poly := Transform(f1,m);
    tg1 : Poly := Transform(g1,m);
    ef1 : Poly := Eliminate(tf1,k);
    eg1 : Poly := Eliminate(tg1,k);

  begin
    if output then
      put("The coordinate transformation for v = "); put(v);
      put(", pivot = "); put(k,1);
      put_line(" :"); put(m);
    end if;
    Common_Initial_Roots(ef1,eg1,v,output,max,eps,tol,t,fail);
    Clear(f1); Clear(g1); Clear(tf1); Clear(tg1); Clear(ef1); Clear(eg1);
  end Initial_Terms;

  function Initial_Terms ( f,g : Poly; t : List ) return List_of_Terms is

    res,res_last : List_of_Terms;
    p : List := t;
    h : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(p) loop
      h := Head_Of(p);
      declare
        s : constant List_of_Terms := Initial_Terms(f,g,h.all);
      begin
        Concat(res,res_last,s);
      end;
      p := Tail_Of(p);
    end loop;
    return res;
  end Initial_Terms;

  procedure Initial_Terms
              ( f,g : in Poly; v : in List; output : in boolean;
                max : in natural32; eps,tol : in double_float;
                t : out List_of_Terms; fail : out boolean ) is

    last : List_of_Terms;
    p : List := v;
    h : Standard_Integer_Vectors.Link_to_Vector;

  begin
    fail := false;
    while not Is_Null(p) loop
      h := Head_Of(p);
      declare
        s : List_of_Terms;
        fl : boolean;
      begin
        Initial_Terms(f,g,h.all,output,max,eps,tol,s,fl);
        fail := fail or fl;
        Concat(t,last,s);
      end;
      p := Tail_Of(p);
    end loop;
  end Initial_Terms;

-- III. evaluate

  function Evaluate ( t,x : Term ) return Term is

    res : Term;
    k : constant integer32 := Standard_Power_Transformations.Pivot(x.dg.all);
    m : constant Standard_Integer_Matrices.Matrix(1..2,1..2)
      := Standard_Power_Transformations.Eliminate(x.dg.all,k);
    exp : integer32;

  begin
    res.dg := new Standard_Integer_Vectors.Vector(1..1);
    if k = 2 then
      res.dg(1) := t.dg(2)*x.dg(2);
      exp := t.dg(1);
    else
      res.dg(1) := t.dg(1)*x.dg(1) + t.dg(2)*x.dg(2);
      exp := t.dg(1)*m(2,1) + t.dg(2)*m(2,2);
    end if;
    res.cf := t.cf*((x.cf)**integer(exp));
    return res;
  end Evaluate;

  function Evaluate ( p : Poly; x : Term ) return Poly is

    res : Poly := Null_Poly;

    procedure Eval_Term ( t : Term; continue : out boolean ) is
 
      v : Term := Evaluate(t,x);

    begin
      Add(res,v);
      Clear(v);
      continue := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Evaluate;

  function Residual ( p : Poly ) return double_float is

    res : double_float := 0.0;

    procedure Scan_Term ( t : Term; continue : out boolean ) is

      v : constant double_float := AbsVal(t.cf);

    begin
      if v > res
       then res := v;
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Residual;

  procedure Split ( f,g : in Poly; t : in List_of_Terms;
                    tol : in double_float; b,r : out List_of_Terms ) is

    tmp : List_of_Terms := t;
    blast,rlast : List_of_Terms;
    vf,vg : Poly;
    fres,gres : double_float;
    trm : Term;

  begin
    while not Is_Null(tmp) loop
      trm := Head_Of(tmp);
      vf := Evaluate(f,trm); fres := Residual(vf);
      vg := Evaluate(g,trm); gres := Residual(vg);
      if fres + gres < tol
       then Append(b,blast,trm);
       else Append(r,rlast,trm);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Split;

end Standard_Binomial_Factors;
