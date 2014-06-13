with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_QR_Least_Squares;  use Standard_Complex_QR_Least_Squares;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;

package body Standard_Polynomial_Interpolators is

-- AUXILIARIES TO SAMPLER :

  function Evaluate ( p : Poly; v : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector with the coefficients for the polynomial
  --   obtained after replacing the first n-1 variables by the values in v.

    res : Standard_Complex_Vectors.Vector(0..Degree(p))
        := (0..Degree(p) => Create(0.0));

    procedure Eval_Term ( t : Term; continue : out boolean ) is

      ind : constant integer32 := integer32(t.dg(t.dg'last));
      acc : Complex_Number := t.cf;

    begin
      for i in v'range loop
        for j in 1..t.dg(i) loop
          acc := acc*v(i);
        end loop;
      end loop;
      res(ind) := res(ind) + acc;
      continue := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Evaluate;

  function Eval ( p : Standard_Complex_Vectors.Vector;
                  x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the function value of the polynomial in x.

    res : Complex_Number := p(p'last);

  begin
    for i in reverse 0..p'last-1 loop
      res := res*x + p(i);
    end loop;
    return res;
  end Eval;

  function Diff ( p : Standard_Complex_Vectors.Vector;
                  x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the derivative of p at x.

    exp : Complex_Number := Create(p'last);
    res : Complex_Number := p(p'last)*exp;

  begin
    for i in reverse 1..p'last-1 loop
      exp := Create(i);
      res := res*x + p(i)*exp;
    end loop;
    return res;
  end Diff;

  procedure Newton ( p : in Standard_Complex_Vectors.Vector;
                     x : in out Complex_Number ) is

  -- DESCRIPTION :
  --   Applies Newton's method to refine the solution of p(x) = 0.
  --   The coefficient of x^i in p is given by p(i).

  -- NOTE : this Newton's method assumes generic inputs.

    eps : constant double_float := 10.0**(-14);
    y,dx : Complex_Number;

  begin
    loop
      y := Eval(p,x);
      exit when (AbsVal(y) < eps);
      dx := y/Diff(p,x);
      exit when (AbsVal(dx) < eps);
      x := x - dx;
    end loop;
  end Newton;

-- AUXILIARIES TO INTERPOLATOR :

  function Minimum ( n1,n2 : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the minimum of n1 and n2.

  begin
    if n1 <= n2
     then return n1;
     else return n2;
    end if;
  end Minimum;

  procedure Coefficient_System
              ( p : in Poly; n : in natural32; v : in VecVec;
                mat : in out Standard_Complex_Matrices.Matrix;
                rhs : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Builds the matrix and right-hand-side vector of the coefficient
  --   system for the polynomial p, interpolating at the points in v.
  --   The dimension of the system is n.

    row,col : integer32;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      col := col+1;
      if col <= integer32(n)
       then mat(row,col) := Eval(t.dg,Create(1.0),v(row).all);
       else rhs(row) := rhs(row) - Eval(t,v(row).all);
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    for i in 1..n loop
      row := integer32(i); col := 0;
      Scan_Terms(p);
    end loop;
  end Coefficient_System;

  function Maximal_Column_Entries
             ( n : natural32; mat : Standard_Complex_Matrices.Matrix )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the maximal absolute value of the elements for every column
  --   in the given matrix.

    res : Standard_Floating_Vectors.Vector(1..integer32(n))
        := (1..integer32(n) => 0.0);

  begin
    for j in mat'range(2) loop
      for i in mat'range(1) loop
        if AbsVal(mat(i,j)) > res(j)
         then res(j) := AbsVal(mat(i,j));
        end if;
      end loop;
    end loop;
    return res;
  end Maximal_Column_Entries;

  procedure Scale_Coefficient_Matrix
              ( mat : in out Standard_Complex_Matrices.Matrix;
                scl : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Every column in the matrix is divided by the corresponding
  --   entry in scl.

  begin
    for j in mat'range(2) loop
      for i in mat'range(1) loop
        mat(i,j) := mat(i,j)/scl(j);
      end loop;
    end loop;
  end Scale_Coefficient_Matrix;

  procedure Scale_Coefficients
              ( cff : in out Standard_Complex_Vectors.Vector;
                scl : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Multiplies every entry in cff with the corresponding entry in scl.

  begin
    for i in cff'range loop
      cff(i) := cff(i)/scl(i);
    end loop;
  end Scale_Coefficients;

  function Assign_Coefficients
             ( p : Poly; c : Standard_Complex_Vectors.Vector ) return Poly is

  -- DESCRIPTION :
  --   The polynomial on return has the same exponents as p,
  --   but with coefficients as in c.

    res : Poly := Null_Poly;
    ind : integer32 := c'first-1;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      nt : Term;

    begin
      ind := ind + 1;
      if ind > c'last then
        Add(res,t);
      else
        nt.cf := c(ind);
        nt.dg := t.dg;
        Add(res,nt);
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Assign_Coefficients;

-- TARGET ROUTINES :

  function Number_of_Terms ( d,n : natural32 ) return natural32 is

    sum : natural32;

  begin
    if n = 2 then
      if d = 0
       then return 1;
       else return ((d+1)*(d+2)/2);
      end if;
    else
      sum := 0;
      for i in 0..d loop
        sum := sum + Number_of_Terms(i,n-1);
      end loop;
      return sum;
    end if;
  end Number_of_Terms;

  function Create ( d,n,cff : natural32 ) return Poly is

    res : Poly := Null_Poly;
    t : Term;

    procedure Generate_Monomials
                ( accu : in out Term; k,sum : in natural32 ) is

    -- DESCRIPTION :
    --   Accumulating procedure to generate all monomials up to degree d.

    -- ON ENTRY :
    --   accu     accumulator contains the current exponent vector;
    --   k        current component;
    --   sum      sum of current entries in the accumulator.

    -- ON RETURN :
    --   accu     accumulator determined up to k component, k included.

      ranflt : double_float;

    begin
      if k > n then
        if cff = 1 then
          accu.cf := Create(1.0);
        elsif cff = 2 then
          ranflt := Random;
          accu.cf := Create(ranflt);
        else
          accu.cf := Random1;
        end if;
        Add(res,accu);
      else
        for i in 0..d loop
          if sum + i <= d then
            accu.dg(integer32(k)) := i;
            Generate_Monomials(accu,k+1,sum+i);
            accu.dg(integer32(k)) := 0;
          end if;
        end loop;
      end if;
    end Generate_Monomials;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    Generate_Monomials(t,1,0);
    return res;
  end Create;

  function Sample ( p : Poly; m : natural32 ) return VecVec is

    res : VecVec(1..integer32(m));
    d : constant integer32 := Degree(p);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    point : Standard_Complex_Vectors.Vector(1..n);
    p1 : Standard_Complex_Vectors.Vector(0..d);

  begin
    for i in 1..integer32(m) loop
      point := Random_Vector(1,n);
      p1 := Evaluate(p,point(1..n-1));
      Newton(p1,point(n));
      res(i) := new Standard_Complex_Vectors.Vector'(point);
    end loop;
    return res;
  end Sample;

  function Sample ( p : Poly; m : natural32;
                    scl : double_float ) return VecVec is

    res : VecVec(1..integer32(m));
    d : constant integer32 := Degree(p);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    point : Standard_Complex_Vectors.Vector(1..n);
    p1 : Standard_Complex_Vectors.Vector(0..d);

  begin
    for i in 1..m loop
      for i in point'range loop
        point(i) := Random1*scl;
      end loop;
      p1 := Evaluate(p,point(1..n-1));
      Newton(p1,point(n));
      res(integer32(i)) := new Standard_Complex_Vectors.Vector'(point);
    end loop;
    return res;
  end Sample;

  procedure Interpolate ( p : in Poly; v : in VecVec;
                          ip : out Poly; rcond : out double_float ) is

    res : Poly;
    n : constant integer32
      := integer32(Minimum(Number_of_Terms(p)-1,natural32(v'length)));
    mat : Standard_Complex_Matrices.Matrix(1..n,1..n);
    rhs : Standard_Complex_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
   -- info : integer;         -- needed for lufac and QRD

  begin
    for i in 1..n loop
      ipvt(i) := i;
      rhs(i) := Create(0.0);
      for j in 1..n loop
        mat(i,j) := Create(0.0);
      end loop;
    end loop;
    Coefficient_System(p,natural32(n),v,mat,rhs);
   -- Gaussian elimination on the coefficient matrix :
   -- lufac(mat,n,ipvt,info);  -- eiter lufac
    lufco(mat,n,ipvt,rcond); -- or lufco -> rcond
    lusolve(mat,n,ipvt,rhs);
   -- put("rcond of coeffmat : "); put(rcond,3); new_line;
    res := Assign_Coefficients(p,rhs);
    ip := res;
  end Interpolate;

  function Interpolate ( p : Poly; v : VecVec ) return Poly is

    res : Poly;
    n : constant integer32
      := integer32(Minimum(Number_of_Terms(p)-1,natural32(v'length)));
    mat : Standard_Complex_Matrices.Matrix(1..n,1..n);
    rhs : Standard_Complex_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    qraux : Standard_Complex_Vectors.Vector(1..n) := (1..n => Create(0.0));
    dum : Standard_Complex_Vectors.Vector(1..n);
    info : integer32;         -- needed for lufac and QRD
   -- maxcols : Standard_Floating_Vectors.Vector(1..n);
   -- rcond : double_float;   -- needed for lufco

  begin
    for i in 1..n loop
      ipvt(i) := i;
      rhs(i) := Create(0.0);
      for j in 1..n loop
        mat(i,j) := Create(0.0);
      end loop;
    end loop;
    Coefficient_System(p,natural32(n),v,mat,rhs);
   -- maxcols := Maximal_Column_Entries(n,mat);
   -- put_line("The maximal column entries :");
   -- put_line(maxcols);
   -- Scale_Coefficient_Matrix(mat,maxcols);
   -- Gaussian elimination on the coefficient matrix :
   -- lufac(mat,n,ipvt,info);  -- eiter lufac
   -- lufco(mat,n,ipvt,rcond); -- or lufco -> rcond
   -- lusolve(mat,n,ipvt,rhs);
   -- put("rcond of coeffmat : "); put(rcond,3); new_line;
   -- QR decomposition followed by Least Squares :
    QRD(mat,qraux,ipvt,false);                         -- QR without pivoting
    QRLS(mat,n,n,qraux,rhs,dum,dum,rhs,dum,dum,100,info);  -- least squares
   -- Scale_Coefficients(rhs,maxcols);
    res := Assign_Coefficients(p,rhs);
    return res;
  end Interpolate;

  function Norm ( p : Poly ) return double_float is

    cff : constant Standard_Complex_Vectors.Vector := Coeff(p);

  begin
    return Max_Norm(cff);
  end Norm;

  function Distance ( p,q : Poly ) return double_float is

    cp : constant Standard_Complex_Vectors.Vector := Coeff(p);
    cq : constant Standard_Complex_Vectors.Vector := Coeff(q);
    use Standard_Complex_Vectors;

  begin
    if cp'last > cq'last
     then return Max_Norm(cp(cq'range) - cq);
     else return Max_Norm(cp - cq(cp'range));
    end if;
  end Distance;

end Standard_Polynomial_Interpolators;
