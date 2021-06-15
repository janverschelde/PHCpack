with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Natural_Vectors;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Embed_Polynomials;        use Standard_Embed_Polynomials;
with Standard_Numerical_Rank;           use Standard_Numerical_Rank;
with Standard_Deflation_Trees_io;       use Standard_Deflation_Trees_io;
with Standard_Deflation_Methods;        use Standard_Deflation_Methods;
with Standard_Nullity_Matrices;         use Standard_Nullity_Matrices;
with Standard_Probe_Kernel;             use Standard_Probe_Kernel;

package body Standard_Multiple_Deflation is

  function Symbolic_Deflate 
              ( nq,nv : natural32;
                a : Standard_Complex_Poly_Matrices.Matrix;
                h : Standard_Complex_Matrices.Matrix ) return Poly_Sys is

    len : constant integer32 := integer32(nq) + a'last(1) + h'last(1);
    res : Poly_Sys(1..len);
    nm : constant integer32 := h'last(2)-1;  -- #multipliers
    mt : Term;
    acc : Poly;
    offset : integer32;

  begin
    for i in 1..integer32(nq) loop
      if a(i,1) /= Null_Poly
       then res(i) := Add_Variables(a(i,1),natural32(nm));
      end if;
    end loop;
    mt.cf := Create(1.0);
    mt.dg := new Standard_Natural_Vectors.Vector'(1..integer32(nv)+nm => 0);
    for i in a'range(1) loop
      res(integer32(nq)+i) := Null_Poly;
      for j in 1..nm loop
        if a(i,j+1) /= Null_Poly then
          mt.dg(integer32(nv)+j) := 1;
          acc := Add_Variables(a(i,j+1),natural32(nm));
          Mul(acc,mt);
          Add(res(integer32(nq)+i),acc);
          Clear(acc);
          mt.dg(integer32(nv)+j) := 0;
        end if;
      end loop;
    end loop;
    offset := integer32(nq)+a'last(1);
    for i in h'range(1) loop
      mt.cf := h(i,1);
      res(offset+i) := Create(mt); 
      for j in 1..nm loop
        mt.dg(integer32(nv)+j) := 1;
        mt.cf := h(i,j+1);
        Add(res(offset+i),mt);
        mt.dg(integer32(nv)+j) := 0;
      end loop;
    end loop;
    Clear(mt);
    return res;
  end Symbolic_Deflate;

  function Symbolic_Deflate
              ( nq,nv,r : natural32; 
                a : Standard_Complex_Poly_Matrices.Matrix ) return Poly_Sys is

    h : constant Standard_Complex_Matrices.Matrix(1..integer32(r),a'range(2))
      := Random_Matrix(r,natural32(a'last(2)));

  begin
    return Symbolic_Deflate(nq,nv,a,h);
  end Symbolic_Deflate;

  function Symbolic_Deflate
              ( p : Poly_Sys; d,r : natural32 ) return Poly_Sys is

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    nr,nc : integer32;

  begin
    Dimensions_of_Nullity_Matrix(nq,nv,d,natural32(nr),natural32(nc));
    declare
      a : constant Standard_Complex_Poly_Matrices.Matrix(1..nr,1..nc)
        := Create_Nullity_Matrix(nq,nv,natural32(nr),natural32(nc),d,p);
      h : constant Standard_Complex_Matrices.Matrix(1..integer32(r),1..nc)
        := Random_Matrix(r,natural32(nc));
    begin
      return Symbolic_Deflate(nq,nv,a,h);
    end;
  end Symbolic_Deflate;

  procedure Predict_Order
              ( p : in Eval_Poly_Sys;
                A : in out Standard_Complex_Matrices.Matrix;
                z : in Vector; tol : in double_float; max_d : in integer32;
                corank : out natural32; d : out natural32 ) is

    U : Standard_Complex_Matrices.Matrix(A'range(1),A'range(1));
    V : Standard_Complex_Matrices.Matrix(A'range(2),A'range(2));
    m : constant integer32 := Min0(A'last(1)+1,A'last(2));
    S : Standard_Complex_Vectors.Vector(1..m);
    rco : double_float;
    rank : natural32;
    w : Vector(z'range);
    t,y,c : Vector(0..max_d);

  begin
    Numerical_Rank(A,tol,S,U,V,rco,rank);
    corank := natural32(A'last(2)) - rank;
    if corank = 0 then
      d := 0;
    else
      w := Random_Vector_in_Kernel(V,corank);
      t := Random_Vector(0,max_d);
      y := Sample_Sum_on_Line(p,z,w,t);
      c := Interpolation_Coefficients(t,y);
      d := Numerical_Order(c,tol);
      if d > 0
       then d := d - 1;
      end if;
    end if;
  end Predict_Order;

  procedure Predict_Order
              ( file : in file_type; p : in Eval_Poly_Sys;
                A : in out Standard_Complex_Matrices.Matrix;
                z : in Vector; tol : in double_float; max_d : in integer32;
                corank : out natural32; d : out natural32 ) is

    use Standard_Complex_Matrices;

    AA : Matrix(A'range(1),A'range(2)) := A;
    U : Matrix(A'range(1),A'range(1));
    V : Matrix(A'range(2),A'range(2));
    m : constant integer32 := Min0(A'last(1)+1,A'last(2));
    S : Vector(1..m);
    rco : double_float;
    rank : natural32;
    w : Vector(z'range);
    r : Vector(A'range(1));
    t,y,c : Vector(0..max_d);

  begin
   -- put_line(file,"The approximate zero : "); put_line(file,z);
    Numerical_Rank(AA,tol,S,U,V,rco,rank);
    put_line(file,"The singular values :");
    put_line(file,S);
    put(file,"The numerical rank : "); put(file,rank,1);
    if integer32(rank) = A'last(2)
     then put(file," = ");
     else put(file," < ");
    end if;
    put(file,A'last(2),1);
    corank := natural32(A'last(2)) - rank;
    put(file,"  Corank is ");
    put(file,corank,1); put(file," w.r.t. tol :");
    put(file,tol,3); put_line(file,".");
    if corank = 0 then
      put_line(file,"Corank is zero, thus also d = 0."); d := 0;
    else
      w := Random_Vector_in_Kernel(V,corank);
     -- put_line(file,"A random vector in the kernel :");
     -- put_line(file,w);
      r := A*w;
      put_line(file,"The residual of the kernel vector : ");
      put_line(file,r);
      t := Random_Vector(0,max_d);
      t(0) := Create(0.0);
      y := Sample_Sum_on_Line(p,z,w,t);
     -- put_line(file,"Values of the sampled sum on the line :"); 
     -- put_line(file,y);
      c := Interpolation_Coefficients(t,y);
      put_line(file,"The coefficients of the interpolating polynomial :");
      put_line(file,c);
      d := Numerical_Order(c,tol);
      if d > 0
       then d := d - 1;
      end if;
      put(file,"The numerical order is ");
      put(file,d,1); put(file," w.r.t. tol :");
      put(file,tol,3); put_line(file,".");
    end if;
  end Predict_Order;

  procedure Predict_Order
              ( p : in Poly_Sys; z : in Vector; tol : in double_float;
                corank : out natural32; d : out natural32 ) is

    use Standard_Complex_Matrices;

    jm : Jaco_Mat(p'range,z'range) := Create(p);
    A : Matrix(p'range,z'range) := Eval(jm,z);
    f : Eval_Poly_Sys(p'range) := Create(p);
    max_d : constant integer32 := Standard_Probe_Kernel.Maximal_Degree(p)+1;

  begin
    Predict_Order(f,A,z,tol,max_d,corank,d);
    Clear(jm); Clear(f);
  end Predict_Order;

  procedure Predict_Order
              ( file : in file_type;
                p : in Poly_Sys; z : in Vector; tol : in double_float;
                corank : out natural32; d : out natural32 ) is

    use Standard_Complex_Matrices;

    jm : Jaco_Mat(p'range,z'range) := Create(p);
    A : Matrix(p'range,z'range) := Eval(jm,z);
    f : Eval_Poly_Sys(p'range) := Create(p);
    max_d : constant integer32 := Standard_Probe_Kernel.Maximal_Degree(p)+1;

  begin
    put(file,"Bound for the maximal order : "); put(file,max_d,1);
    put_line(file,".");
    Predict_Order(file,f,A,z,tol,max_d,corank,d);
    Clear(jm); Clear(f);
  end Predict_Order;

  procedure Numeric_Deflate
              ( p : in Poly_Sys; z : in Vector;
                d : in natural32; tol : in double_float;
                r : out natural32; dp : out Link_to_Poly_Sys ) is

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    nr,nc : integer32;

  begin
    Dimensions_of_Nullity_Matrix(nq,nv,d,natural32(nr),natural32(nc));
    declare
      a : constant Standard_Complex_Poly_Matrices.Matrix(1..nr,1..nc)
        := Create_Nullity_Matrix(nq,nv,natural32(nr),natural32(nc),d,p);
      y : constant Standard_Complex_Matrices.Matrix(1..nr,1..nc-1)
        := Eval1(a,z);
      rank : constant natural32 := Numerical_Rank(y,tol);
    begin
      r := natural32(nc) - rank;
      dp := new Poly_Sys'(Symbolic_Deflate(nq,nv,r,a));
    end;
  end Numeric_Deflate;

  procedure Numeric_Deflate
              ( file : in file_type; p : in Poly_Sys; z : in Vector;
                d : in natural32; tol : in double_float;
                r : out natural32; dp : out Link_to_Poly_Sys ) is

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    nr,nc : natural32;

  begin
    Dimensions_of_Nullity_Matrix(nq,nv,d,nr,nc);
    put(file,"Creating a "); put(file,nr,1);
    put(file,"-by-"); put(file,nc,1); put_line(file," nullity matrix...");
    declare
      a : constant Standard_Complex_Poly_Matrices.Matrix
            (1..integer32(nr),1..integer32(nc))
        := Create_Nullity_Matrix(nq,nv,nr,nc,d,p);
      y : Standard_Complex_Matrices.Matrix
            (1..integer32(nr),1..integer32(nc)-1) := Eval1(a,z);
      u : Standard_Complex_Matrices.Matrix(y'range(1),y'range(1));
      v : Standard_Complex_Matrices.Matrix(y'range(2),y'range(2));
      m : constant integer32 := Min0(integer32(nr)+1,integer32(nc)-1);
      s : Standard_Complex_Vectors.Vector(1..m);
      rco : double_float;
      rank : natural32;
    begin
      Numerical_Rank(y,tol,s,u,v,rco,rank);
      put_line(file,"The singular values :");
      put_line(file,s);
      put(file,"The numerical rank is "); put(file,rank,1); new_line(file);
      r := nc - rank;
      put(file,"The corank is "); put(file,r,1); new_line(file);
      dp := new Poly_Sys'(Symbolic_Deflate(nq,nv,r,a));
    end;
  end Numeric_Deflate;

  procedure Add_Multipliers 
              ( wz : in out Link_to_Vector; nm : in natural32 ) is

    nz : Vector(1..wz'last+integer32(nm));

  begin
    nz(wz'range) := wz.all;
    for i in wz'last+1..nz'last loop
      nz(i) := Create(0.0);
    end loop;
    Clear(wz);
    wz := new Vector'(nz);
  end Add_Multipliers;

  procedure Interactive_Symbolic_Deflation
              ( file : in file_type; p : in Poly_Sys;
                sol : in out Solution; tol : in double_float ) is

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    wz : Link_to_Vector := new Vector'(sol.v);
    err,rco,res : double_float;
    wp,dp : Link_to_Poly_Sys;
    rank,d,nm,corank,k : natural32 := 0;
    ans : character;

  begin
    wp := new Poly_Sys(p'range);
    Copy(p,wp.all);
    loop
      Interactive_Symbolic_Newton(file,wp.all,wz.all,err,rco,res,tol,rank);
      sol.v := wz(sol.v'range);
      sol.err := err;
      sol.rco := rco;
      sol.res := res;
      Predict_Order(standard_output,wp.all,wz.all,tol,corank,d);
      put("Do you want to deflate ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the order of the deflation : "); get(d);
        Numeric_Deflate(standard_output,wp.all,wz.all,d,tol,rank,dp);
        Clear(wp); wp := new Poly_Sys(dp'range);
        wp.all := dp.all;
        nm := Number_of_Unknowns(wp(wp'first)) - n;
        Add_Multipliers(wz,nm);
        k := k + 1;
        Add_Multiplier_Symbols(k,nm);
      else
        put("Do you want to continue with this root ? (y/n) ");
        Ask_Yes_or_No(ans);
      end if;
      exit when (ans /= 'y');
    end loop;
  end Interactive_Symbolic_Deflation;

  procedure Interactive_Symbolic_Deflation
              ( file : in file_type; p : in Poly_Sys;
                sols : in out Solution_List; tol : in double_float ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Interactive_Symbolic_Deflation(file,p,ls.all,tol);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Interactive_Symbolic_Deflation;

end Standard_Multiple_Deflation;
