with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;

package body Singular_Values_of_Hessians is

-- AUXILIARY WRAPPERS :

  procedure Singular_Values
             ( A : in out Standard_Complex_Matrices.Matrix;
               s : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the singular values of A and returns the result in s,
  --   computed in standard double precision.

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

  -- DESCRIPTION :
  --   Computes the singular values of A and returns the result in s,
  --   computed in double double precision.

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

  -- DESCRIPTION :
  --   Computes the singular values of A and returns the result in s,
  --   computed in quad double precision.

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

-- TARGET FUNCTIONS :

  function Standard_Singular_Values
             ( h : Standard_Complex_Hessians.Link_to_Hessian;
               x : Standard_Complex_Vectors.Vector )
             return  Standard_Complex_Vectors.Vector is

    A : Standard_Complex_Matrices.Matrix(h'range(1),h'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    d : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    res : Standard_Complex_Vectors.Vector(1..d);

  begin
    A := Standard_Complex_Hessians.Eval(h,x);
    Singular_Values(A,res);
    return res;
  end Standard_Singular_Values;

  function DoblDobl_Singular_Values
             ( h : DoblDobl_Complex_Hessians.Link_to_Hessian;
               x : DoblDobl_Complex_Vectors.Vector )
             return  DoblDobl_Complex_Vectors.Vector is

    A : DoblDobl_Complex_Matrices.Matrix(h'range(1),h'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    d : constant integer32 := DoblDobl_Complex_Singular_Values.Min0(n+1,p);
    res : DoblDobl_Complex_Vectors.Vector(1..d);

  begin
    A := DoblDobl_Complex_Hessians.Eval(h,x);
    Singular_Values(A,res);
    return res;
  end DoblDobl_Singular_Values;

  function QuadDobl_Singular_Values
             ( h : QuadDobl_Complex_Hessians.Link_to_Hessian;
               x : QuadDobl_Complex_Vectors.Vector )
             return  QuadDobl_Complex_Vectors.Vector is

    A : QuadDobl_Complex_Matrices.Matrix(h'range(1),h'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    d : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    res : QuadDobl_Complex_Vectors.Vector(1..d);

  begin
    A := QuadDobl_Complex_Hessians.Eval(h,x);
    Singular_Values(A,res);
    return res;
  end QuadDobl_Singular_Values;

end Singular_Values_of_Hessians;
