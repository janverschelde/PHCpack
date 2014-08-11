with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;    use QuadDobl_Complex_Numbers_Polar;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Linear_Solvers;   use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Random_Vectors;           use QuadDobl_Random_Vectors;
with QuadDobl_Random_Matrices;          use QuadDobl_Random_Matrices;

package body QuadDobl_Plane_Representations is

-- AUXILIARY OPERATION FOR COMPUTING GENERATORS :

  function Pivot ( v : Vector; a,b : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of the largest entry in v, between a and b.

    max : constant quad_double := Radius(v(a));
    res : integer32 := a;

  begin
    for i in a+1..b loop
      if Radius(v(i)) > max
       then res := i;
      end if;
    end loop;
    return res;
  end Pivot;

-- DATA STRUCTURE CONVERSIONS :

  function Equations_to_Matrix ( c : VecVec; n : integer32 ) return Matrix is

    res : Matrix(c'range,0..n);

  begin
    for i in c'range loop
      for j in 0..n loop
        res(i,j) := c(i)(j);
      end loop;
    end loop;
    return res;
  end Equations_to_Matrix;

  function Equations_to_Matrix ( c : VecVec ) return Matrix is

    cf : constant Link_to_Vector := c(c'first);

  begin
    return Equations_to_Matrix(c,cf'last);
  end Equations_to_Matrix;

  function Equations_to_VecVec ( c : Matrix ) return VecVec is

    res : VecVec(c'range(1));
    n : constant integer32 := c'last(2);

  begin
    for i in c'range(1) loop
      res(i) := new QuadDobl_Complex_Vectors.Vector(0..n);
      for j in c'range(2) loop
        res(i)(j) := c(i,j);
      end loop;
    end loop;
    return res;
  end Equations_to_VecVec;

  function Generators_to_Matrix ( b : Vector; v : VecVec ) return Matrix is

    res : Matrix(b'range,0..v'last);

  begin
    for i in b'range loop
      res(i,0) := b(i);
      for j in v'range loop
        res(i,j) := v(j)(i);
      end loop;
    end loop;
    return res;
  end Generators_to_Matrix;

  procedure Generators_to_VecVec
              ( g : in Matrix; b : out Vector; v : out VecVec ) is
  begin
    for i in b'range loop
      b(i) := g(i,0);
    end loop;
    for j in v'range loop
      v(j) := new QuadDobl_Complex_Vectors.Vector(b'range);
      for i in g'range(1) loop
        v(j)(i) := g(i,j);
      end loop;
    end loop;
  end Generators_to_VecVec;

-- FROM EQUATIONS TO GENERATORS :

  procedure Generators1 ( hyp : in Vector;
                          basis : out Vector; directions : out VecVec ) is

    n : constant integer32 := hyp'last;
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    sum : Complex_Number := Create(zero);
    piv,cnt : integer32;

  begin
    basis := Random_Vector(1,n);
    piv := Pivot(hyp,1,n);
    sum := hyp(0);
    for i in 1..n loop
      if i /= piv
       then sum := sum + hyp(i)*basis(i);
      end if;
    end loop;
    basis(piv) := -sum/hyp(piv);
    cnt := 0;
    for i in 1..n loop
      if i /= piv then
        cnt := cnt + 1;
        directions(cnt)
          := new QuadDobl_Complex_Vectors.Vector'(1..n => Create(zero));
        directions(cnt)(i) := Create(one);
        directions(cnt)(piv) := -hyp(i)/hyp(piv);
      end if;
    end loop;
  end Generators1;

  procedure Generators ( n,k : in integer32; hyp : in VecVec;
                         basis : out Vector; directions : out VecVec ) is

    cff : Matrix(1..k,1..n+1);
    mat : Matrix(1..k,1..k);
    rhs : QuadDobl_Complex_Vectors.Vector(1..k);
    ipvt : Standard_Integer_Vectors.Vector(1..k);
    info,cnt : integer32;
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    for i in 1..k loop
      for j in 1..n loop
        cff(i,j) := hyp(i)(j);
      end loop;
      cff(i,n+1) := hyp(i)(0);
      for j in 1..k loop
        mat(i,j) := cff(i,j);
      end loop;
    end loop;
    lufac(mat,k,ipvt,info);
    for i in 1..k loop
      rhs(i) := -cff(i,n+1);
    end loop;
    lusolve(mat,k,ipvt,rhs);
    basis(rhs'range) := rhs;
    basis(k+1..n) := (k+1..n => Create(zero));
    cnt := 0;
    for i in k+1..n loop
      for j in 1..k loop
        rhs(j) := -cff(j,i);
      end loop;
      lusolve(mat,k,ipvt,rhs);
      cnt := cnt+1;
      directions(cnt)
        := new QuadDobl_Complex_Vectors.Vector'(1..n => Create(zero));
      for j in 1..k loop
        directions(cnt)(j) := rhs(j);
      end loop;
      directions(cnt)(i) := Create(one);
    end loop;
  end Generators;

  function Generators ( hyp : Matrix ) return Matrix is

    k : constant integer32 := hyp'last(1);
    n : constant integer32 := hyp'last(2);
    mat : Matrix(1..k,1..k);
    res : Matrix(1..n,0..n-k);
    ipvt : Standard_Integer_Vectors.Vector(1..k);
    info : integer32;
    b : Vector(1..k);
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in 1..k loop
      for j in 1..k loop
        mat(i,j) := hyp(i,j);
      end loop;
      b(i) := -hyp(i,0);
    end loop;
    lufac(mat,k,ipvt,info);
    lusolve(mat,k,ipvt,b);
    for i in b'range loop
      res(i,0) := b(i);
    end loop;
    for i in k+1..n loop
      res(i,0) := Create(zero);
    end loop;
    for j in k+1..n loop
      for i in b'range loop
        b(i) := -hyp(i,j);
      end loop;
      lusolve(mat,k,ipvt,b);
      for i in b'range loop
        res(i,j-k) := b(i); 
      end loop;
      for i in k+1..n loop
        res(i,j-k) := Create(zero);
      end loop;
      res(j,j-k) := Create(one);
    end loop;
    return res;
  end Generators;

  function Orthogonalize ( v : Matrix ) return Matrix is

    res : Matrix(v'range(1),v'range(2));
    m1,m2 : Matrix(v'range(1),1..v'last(2));

  begin
    for i in v'range(1) loop
      for j in 1..v'last(2) loop
        m1(i,j) := v(i,j);
      end loop;
    end loop;
    m2 := QuadDobl_Random_Matrices.Orthogonalize(m1);
    for i in v'range(1) loop
      res(i,0) := v(i,0);
      for j in 1..v'last(2) loop
        res(i,j) := m2(i,j);
      end loop;
    end loop;
    return res;
  end Orthogonalize;

  procedure Orthogonalize ( v : in out Matrix ) is

    m1,m2 : Matrix(v'range(1),1..v'last(2));

  begin
    for i in v'range(1) loop
      for j in 1..v'last(2) loop
        m1(i,j) := v(i,j);
      end loop;
    end loop;
    m2 := QuadDobl_Random_Matrices.Orthogonalize(m1);
    for i in v'range(1) loop
      for j in 1..v'last(2) loop
        v(i,j) := m2(i,j);
      end loop;
    end loop;
  end Orthogonalize;

  function Orthogonalize ( v : VecVec ) return VecVec is

    res : VecVec(v'range);
    n : constant integer32 := v(v'first)'length;
    m1,m2 : Matrix(1..n,v'range);

  begin
    for j in v'range loop
      for i in v(j)'range loop
        m1(i,j) := v(j)(i);
      end loop;
    end loop;
    m2 := QuadDobl_Random_Matrices.Orthogonalize(m1);
    for j in res'range loop
      res(j) := new QuadDobl_Complex_vectors.Vector(1..n);
      for i in 1..n loop
        res(j)(i) := m2(i,j);
      end loop;
    end loop;
    return res;
  end Orthogonalize;

-- FROM GENERATORS TO EQUATIONS :

  function Equations1 ( basis,direction : Vector ) return VecVec is

    n : constant integer32 := basis'last;
    res : VecVec(1..n-1);
    piv : constant integer32 := Pivot(direction,1,n);
    cnt : integer32 := 0;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in 1..n-1 loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'(0..n => Create(zero));
    end loop;
    for i in 1..n loop
      if i /= piv then
        cnt := cnt+1;
        res(cnt)(i) := Create(one);
        res(cnt)(piv) := -direction(i)/direction(piv);
        res(cnt)(0) := -basis(i) - basis(piv)*res(cnt)(piv);
      end if;
    end loop;
    return res;
  end Equations1;

  function Equations ( basis : Vector; directions : VecVec ) return VecVec is

    n : constant integer32 := basis'last;
    k : constant integer32 := directions'last;
    res : VecVec(1..n-k);
    cff : Matrix(1..k,1..n);
    mat : Matrix(1..k,1..k);
    rhs : QuadDobl_Complex_Vectors.Vector(1..k);
    ipvt : Standard_Integer_Vectors.Vector(1..k);
    info,cnt : integer32;
    sum : Complex_Number;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in 1..n-k loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'(0..n => Create(zero));
    end loop;
    for i in 1..k loop
      for j in 1..n loop
        cff(i,j) := directions(i)(j);
      end loop;
    end loop;
    for i in 1..k loop
      for j in 1..k loop
        mat(i,j) := cff(i,j);
      end loop;
    end loop;
    lufac(mat,k,ipvt,info);
    cnt := 0;
    for i in k+1..n loop
      for j in 1..k loop
        rhs(j) := -cff(j,i); 
      end loop;
      lusolve(mat,k,ipvt,rhs);
      cnt := cnt+1;
      for j in 1..k loop
        res(cnt)(j) := rhs(j);
      end loop;
      res(cnt)(i) := Create(one);
      sum := basis(i);
      for j in 1..k loop
        sum := sum + rhs(j)*basis(j);
      end loop;
      res(cnt)(0) := -sum;
    end loop;
    return res;
  end Equations;

  function Equations ( g : Matrix ) return Matrix is

    n : constant integer32 := g'last(1);
    k : constant integer32 := g'last(2);
    res : Matrix(1..n-k,0..n);
    mat : Matrix(1..k,1..k);
    b : QuadDobl_Complex_Vectors.Vector(1..k);
    ipvt : Standard_Integer_Vectors.Vector(1..k);
    info : integer32;
    sum : Complex_Number;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in 1..k loop
      for j in 1..k loop
        mat(i,j) := g(j,i);
      end loop;
    end loop;
    lufac(mat,k,ipvt,info);
    for j in k+1..n loop
      for i in 1..k loop
        b(i) := -g(j,i); 
      end loop;
      lusolve(mat,k,ipvt,b);
      for i in 1..k loop
        res(j-k,i) := b(i);
      end loop;
      for i in k+1..n loop
        res(j-k,i) := Create(zero);
      end loop;
      res(j-k,j) := Create(one);
      sum := g(j,0);
      for i in 1..k loop
        sum := sum + b(i)*g(i,0);
      end loop;
      res(j-k,0) := -sum;
    end loop;
    return res;
  end Equations;

end QuadDobl_Plane_Representations;
