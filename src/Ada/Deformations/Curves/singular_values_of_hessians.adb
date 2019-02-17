with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Floating_Two_Norms;
with Double_Double_Two_Norms;
with Quad_Double_Two_Norms;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;

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
    return sigma1/nrm;
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
    return sigma1/nrm;
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
    return sigma1/nrm;
  end QuadDobl_Distance;

end Singular_Values_of_Hessians;
