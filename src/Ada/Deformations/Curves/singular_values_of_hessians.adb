with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Floating_Two_Norms;
with Double_Double_Two_Norms;
with Quad_Double_Two_Norms;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Homotopy;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Homotopy;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Homotopy;

package body Singular_Values_of_Hessians is

  procedure Singular_Values
             ( A : in out Standard_Complex_Matrices.Matrix;
               s : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Standard_Complex_Matrices.Matrix(1..n,1..n);
    v : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Standard_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,job,info);
  end Singular_Values;

  procedure Singular_Values
             ( A : in out DoblDobl_Complex_Matrices.Matrix;
               s : out DoblDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    e : DoblDobl_Complex_Vectors.Vector(1..p);
    u : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : DoblDobl_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    DoblDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,job,info);
  end Singular_Values;

  procedure Singular_Values
             ( A : in out QuadDobl_Complex_Matrices.Matrix;
               s : out QuadDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    u : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : QuadDobl_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    QuadDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,job,info);
  end Singular_Values;

  function Standard_Singular_Values
             ( h : Standard_Complex_Hessians.Link_to_Hessian;
               x : Standard_Complex_Vectors.Vector )
             return  Standard_Floating_Vectors.Vector is

    A : Standard_Complex_Matrices.Matrix(h'range(1),h'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    d : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Standard_Complex_Vectors.Vector(1..d);
    res : Standard_Floating_Vectors.Vector(1..d);

  begin
    A := Standard_Complex_Hessians.Eval(h,x);
    Singular_Values(A,s);
    for i in s'range loop
      res(i) := Standard_Complex_Numbers.Real_Part(s(i));
    end loop;
    return res;
  end Standard_Singular_Values;

  function DoblDobl_Singular_Values
             ( h : DoblDobl_Complex_Hessians.Link_to_Hessian;
               x : DoblDobl_Complex_Vectors.Vector )
             return  Double_Double_Vectors.Vector is

    A : DoblDobl_Complex_Matrices.Matrix(h'range(1),h'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    d : constant integer32 := DoblDobl_Complex_Singular_Values.Min0(n+1,p);
    s : DoblDobl_Complex_Vectors.Vector(1..d);
    res : Double_Double_Vectors.Vector(1..d);

  begin
    A := DoblDobl_Complex_Hessians.Eval(h,x);
    Singular_Values(A,s);
    for i in s'range loop
      res(i) := DoblDobl_Complex_Numbers.Real_Part(s(i));
    end loop;
    return res;
  end DoblDobl_Singular_Values;

  function QuadDobl_Singular_Values
             ( h : QuadDobl_Complex_Hessians.Link_to_Hessian;
               x : QuadDobl_Complex_Vectors.Vector )
             return  Quad_Double_Vectors.Vector is

    A : QuadDobl_Complex_Matrices.Matrix(h'range(1),h'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    d : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    s : QuadDobl_Complex_Vectors.Vector(1..d);
    res : Quad_Double_Vectors.Vector(1..d);

  begin
    A := QuadDobl_Complex_Hessians.Eval(h,x);
    Singular_Values(A,s);
    for i in s'range loop
      res(i) := QuadDobl_Complex_Numbers.Real_Part(s(i));
    end loop;
    return res;
  end QuadDobl_Singular_Values;

   function Standard_Singular_Values
              ( h : Standard_Complex_Hessians.Array_of_Hessians;
                x : Standard_Complex_Vectors.Vector )
              return Standard_Floating_Vectors.Vector is

     res : Standard_Floating_Vectors.Vector(h'range);

   begin
     for i in h'range loop
       declare
         wrk : constant Standard_Floating_Vectors.Vector
             := Standard_Singular_Values(h(i),x);
       begin
         res(i) := wrk(wrk'first);
       end;
     end loop;
     return res;
   end Standard_Singular_Values;

   function DoblDobl_Singular_Values
              ( h : DoblDobl_Complex_Hessians.Array_of_Hessians;
                x : DoblDobl_Complex_Vectors.Vector )
              return Double_Double_Vectors.Vector is

     res : Double_Double_Vectors.Vector(h'range);

   begin
     for i in h'range loop
       declare
         wrk : constant Double_Double_Vectors.Vector
             := DoblDobl_Singular_Values(h(i),x);
       begin
         res(i) := wrk(wrk'first);
       end;
     end loop;
     return res;
   end DoblDobl_Singular_Values;

   function QuadDobl_Singular_Values
              ( h : QuadDobl_Complex_Hessians.Array_of_Hessians;
                x : QuadDobl_Complex_Vectors.Vector )
              return Quad_Double_Vectors.Vector is

     res : Quad_Double_Vectors.Vector(h'range);

   begin
     for i in h'range loop
       declare
         wrk : constant Quad_Double_Vectors.Vector
             := QuadDobl_Singular_Values(h(i),x);
       begin
         res(i) := wrk(wrk'first);
       end;
     end loop;
     return res;
   end QuadDobl_Singular_Values;

  function Standard_Distance
             ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in Standard_Complex_Hessians.Array_of_Hessians;
               xt : in Standard_Complex_Vectors.Vector )
             return double_float is

    jmx : Standard_Complex_Matrices.Matrix(jm'range(1),jm'range(2));
    n : constant integer32 := jm'last(1);
    p : constant integer32 := jm'last(2);
    dim : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    jsv : Standard_Complex_Vectors.Vector(1..dim);
    nrm,sigma1 : double_float;
    sv : constant Standard_Floating_Vectors.Vector
       := Standard_Singular_Values(hs,xt);

  begin
    jmx := Standard_Complex_Jaco_Matrices.Eval(jm,xt);
    Singular_Values(jmx,jsv);
    sigma1 := Standard_Complex_Numbers.REAL_PART(jsv(jsv'last));
    nrm := Standard_Floating_Two_Norms.Norm2(sv);
    return (2.0*sigma1)/nrm;
  end Standard_Distance;

  function DoblDobl_Distance
             ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in DoblDobl_Complex_Hessians.Array_of_Hessians;
               xt : in DoblDobl_Complex_Vectors.Vector )
             return double_double is

    jmx : DoblDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2));
    n : constant integer32 := jm'last(1);
    p : constant integer32 := jm'last(2);
    dim : constant integer32 := DoblDobl_Complex_Singular_Values.Min0(n+1,p);
    jsv : DoblDobl_Complex_Vectors.Vector(1..dim);
    nrm,sigma1 : double_double;
    sv : constant Double_Double_Vectors.Vector
       := DoblDobl_Singular_Values(hs,xt);

  begin
    jmx := DoblDobl_Complex_Jaco_Matrices.Eval(jm,xt);
    Singular_Values(jmx,jsv);
    sigma1 := DoblDobl_Complex_Numbers.REAL_PART(jsv(jsv'last));
    nrm := Double_Double_Two_Norms.Norm2(sv);
    return (2.0*sigma1)/nrm;
  end DoblDobl_Distance;

  function QuadDobl_Distance
             ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in QuadDobl_Complex_Hessians.Array_of_Hessians;
               xt : in QuadDobl_Complex_Vectors.Vector )
             return quad_double is

    jmx : QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2));
    n : constant integer32 := jm'last(1);
    p : constant integer32 := jm'last(2);
    dim : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    jsv : QuadDobl_Complex_Vectors.Vector(1..dim);
    nrm,sigma1 : quad_double;
    sv : constant Quad_Double_Vectors.Vector
       := QuadDobl_Singular_Values(hs,xt);

  begin
    jmx := QuadDobl_Complex_Jaco_Matrices.Eval(jm,xt);
    Singular_Values(jmx,jsv);
    sigma1 := QuadDobl_Complex_Numbers.REAL_PART(jsv(jsv'last));
    nrm := Quad_Double_Two_Norms.Norm2(sv);
    return (2.0*sigma1)/nrm;
  end QuadDobl_Distance;

  function Standard_Distance
             ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in Standard_Complex_Hessians.Array_of_Hessians;
               sol : in Standard_Complex_Solutions.Solution )
             return double_float is

    xt : Standard_Complex_Vectors.Vector(1..sol.n+1);

  begin
    xt(sol.v'range) := sol.v;
    xt(xt'last) := sol.t;
    return Standard_Distance(jm,hs,xt);
  end Standard_Distance;

  function DoblDobl_Distance
             ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in DoblDobl_Complex_Hessians.Array_of_Hessians;
               sol : in DoblDobl_Complex_Solutions.Solution )
             return double_double is

    xt : DoblDobl_Complex_Vectors.Vector(1..sol.n+1);

  begin
    xt(sol.v'range) := sol.v;
    xt(xt'last) := sol.t;
    return DoblDobl_Distance(jm,hs,xt);
  end DoblDobl_Distance;

  function QuadDobl_Distance
             ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in QuadDobl_Complex_Hessians.Array_of_Hessians;
               sol : in QuadDobl_Complex_Solutions.Solution )
             return quad_double is

    xt : QuadDobl_Complex_Vectors.Vector(1..sol.n+1);

  begin
    xt(sol.v'range) := sol.v;
    xt(xt'last) := sol.t;
    return QuadDobl_Distance(jm,hs,xt);
  end QuadDobl_Distance;

  procedure Standard_Jacobian_Hessians_of_Homotopy
              ( jm : out Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out Standard_Complex_Hessians.Link_to_Array_of_Hessians ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys
      := Standard_Homotopy.Homotopy_System;
    nbeq : constant integer32 := p'last;
    n : constant integer32
      := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    hs : constant Standard_Complex_Hessians.Array_of_Hessians(p'range)
       := Standard_Complex_Hessians.Create(p,n);
    cm : Standard_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n1 loop
        cm(i,j) := Standard_Complex_Polynomials.Diff(p(i),j);
      end loop;
    end loop;
    jm := new Standard_Complex_Jaco_Matrices.Jaco_Mat'(cm);
    h := new Standard_Complex_Hessians.Array_of_Hessians'(hs);
  end Standard_Jacobian_Hessians_of_Homotopy;

  procedure DoblDobl_Jacobian_Hessians_of_Homotopy
              ( jm : out DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians ) is

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
      := DoblDobl_Homotopy.Homotopy_System;
    nbeq : constant integer32 := p'last;
    n : constant integer32
      := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    hs : constant DoblDobl_Complex_Hessians.Array_of_Hessians(p'range)
       := DoblDobl_Complex_Hessians.Create(p,n);
    cm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n1 loop
        cm(i,j) := DoblDobl_Complex_Polynomials.Diff(p(i),j);
      end loop;
    end loop;
    jm := new DoblDobl_Complex_Jaco_Matrices.Jaco_Mat'(cm);
    h := new DoblDobl_Complex_Hessians.Array_of_Hessians'(hs);
  end DoblDobl_Jacobian_Hessians_of_Homotopy;

  procedure QuadDobl_Jacobian_Hessians_of_Homotopy
              ( jm : out QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians ) is

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
      := QuadDobl_Homotopy.Homotopy_System;
    nbeq : constant integer32 := p'last;
    n : constant integer32
      := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    hs : constant QuadDobl_Complex_Hessians.Array_of_Hessians(p'range)
       := QuadDobl_Complex_Hessians.Create(p,n);
    cm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n1 loop
        cm(i,j) := QuadDobl_Complex_Polynomials.Diff(p(i),j);
      end loop;
    end loop;
    jm := new QuadDobl_Complex_Jaco_Matrices.Jaco_Mat'(cm);
    h := new QuadDobl_Complex_Hessians.Array_of_Hessians'(hs);
  end QuadDobl_Jacobian_Hessians_of_Homotopy;

  procedure Standard_Jacobian_Hessians_of_Homotopy
              ( k : in integer32;
                jm : out Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out Standard_Complex_Hessians.Link_to_Array_of_Hessians ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys
      := Standard_Homotopy.Homotopy_System;
    nbeq : constant integer32 := p'last;
    n : constant integer32
      := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    hs : constant Standard_Complex_Hessians.Array_of_Hessians(p'range)
       := Standard_Complex_Hessians.Create(p,k);
    cm : Standard_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n loop
        if j < k then
          cm(i,j) := Standard_Complex_Polynomials.Diff(p(i),j);
        elsif j > k then
          cm(i,j-1) := Standard_Complex_Polynomials.Diff(p(i),j);
        end if;
      end loop;
    end loop;
    jm := new Standard_Complex_Jaco_Matrices.Jaco_Mat'(cm);
    h := new Standard_Complex_Hessians.Array_of_Hessians'(hs);
  end Standard_Jacobian_Hessians_of_Homotopy;

  procedure DoblDobl_Jacobian_Hessians_of_Homotopy
              ( k : in integer32;
                jm : out DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians ) is

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
      := DoblDobl_Homotopy.Homotopy_System;
    nbeq : constant integer32 := p'last;
    n : constant integer32
      := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    hs : constant DoblDobl_Complex_Hessians.Array_of_Hessians(p'range)
       := DoblDobl_Complex_Hessians.Create(p,k);
    cm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n loop
        if j < k then
          cm(i,j) := DoblDobl_Complex_Polynomials.Diff(p(i),j);
        elsif i > k then
          cm(i,j-1) := DoblDobl_Complex_Polynomials.Diff(p(i),j);
        end if;
      end loop;
    end loop;
    jm := new DoblDobl_Complex_Jaco_Matrices.Jaco_Mat'(cm);
    h := new DoblDobl_Complex_Hessians.Array_of_Hessians'(hs);
  end DoblDobl_Jacobian_Hessians_of_Homotopy;

  procedure QuadDobl_Jacobian_Hessians_of_Homotopy
              ( k : in integer32;
                jm : out QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians ) is

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
      := QuadDobl_Homotopy.Homotopy_System;
    nbeq : constant integer32 := p'last;
    n : constant integer32
      := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    hs : constant QuadDobl_Complex_Hessians.Array_of_Hessians(p'range)
       := QuadDobl_Complex_Hessians.Create(p,k);
    cm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n loop
        if j < k then
          cm(i,j) := QuadDobl_Complex_Polynomials.Diff(p(i),j);
        elsif j > k then
          cm(i,j-1) := QuadDobl_Complex_Polynomials.Diff(p(i),j);
        end if;
      end loop;
    end loop;
    jm := new QuadDobl_Complex_Jaco_Matrices.Jaco_Mat'(cm);
    h := new QuadDobl_Complex_Hessians.Array_of_Hessians'(hs);
  end QuadDobl_Jacobian_Hessians_of_Homotopy;

end Singular_Values_of_Hessians;
