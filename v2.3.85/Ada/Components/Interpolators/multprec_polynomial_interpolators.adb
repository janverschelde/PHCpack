with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Complex_Matrices;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
with Multprec_Complex_QR_Least_Squares;  use Multprec_Complex_QR_Least_Squares;
with Multprec_Complex_Poly_Functions;    use Multprec_Complex_Poly_Functions;

package body Multprec_Polynomial_Interpolators is

-- AUXILIARIES TO SAMPLER :

  function Evaluate ( p : Poly; v : Multprec_Complex_Vectors.Vector )
                    return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector with the coefficients for the polynomial
  --   obtained after replacing the first n-1 variables by the values in v.

    res : Multprec_Complex_Vectors.Vector(0..Degree(p))
        := (0..Degree(p) => Create(integer(0)));

    procedure Eval_Term ( t : Term; continue : out boolean ) is

      ind : constant integer32 := integer32(t.dg(t.dg'last));
      acc : Complex_Number;

    begin
      Copy(t.cf,acc);
      for i in v'range loop
        for j in 1..t.dg(i) loop
         -- acc := acc*v(i);
          Mul(acc,v(i));
        end loop;
      end loop;
     -- res(ind) := res(ind) + acc;
      Add(res(ind),acc);
      continue := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Evaluate;

  function Eval ( p : Multprec_Complex_Vectors.Vector;
                  x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the function value of the polynomial in x.

    res : Complex_Number;

  begin
    Copy(p(p'last),res);
    for i in reverse 0..p'last-1 loop
     -- res := res*x + p(i);
      Mul(res,x);
      Add(res,p(i));
    end loop;
    return res;
  end Eval;

  function Diff ( p : Multprec_Complex_Vectors.Vector;
                  x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the derivative of p at x.

    exp : Complex_Number := Create(p'last);
    res : Complex_Number := p(p'last)*exp;
    fac : Complex_Number;

  begin
    for i in reverse 1..p'last-1 loop
      exp := Create(i);
     -- res := res*x + p(i)*exp;
      Mul(res,x);
      fac := p(i)*exp;
      Add(res,fac);
      Clear(fac);
    end loop;
    return res;
  end Diff;

  procedure Newton ( p : in Multprec_Complex_Vectors.Vector;
                     x : in out Complex_Number ) is

  -- DESCRIPTION :
  --   Applies Newton's method to refine the solution of p(x) = 0.
  --   The coefficient of x^i in p is given by p(i).

  -- NOTE : this Newton's method assumes generic inputs.

    eps : constant Floating_Number := Create(1.0E-100);
    y,dx,df : Complex_Number;

  begin
    for i in 1..50 loop
      y := Eval(p,x);
      exit when (AbsVal(y) < eps);
      df := Diff(p,x);
      dx := y/df;
      exit when (AbsVal(dx) < eps);
      Sub(x,dx);
      Clear(y); Clear(dx); Clear(df);
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

  procedure Coefficient_System1
              ( p : in Poly; n : in natural32; v : in VecVec;
                mat : in out Multprec_Complex_Matrices.Matrix;
                rhs : in out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Builds the matrix and right-hand-side vector of the coefficient
  --   system for the polynomial p, interpolating at the points in v.
  --   The dimension of the system is n.
  --   Here the constant term of p is assumed to be nonzero.

    row,col : integer32;
    one : Complex_Number := Create(integer(1));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      eva : Complex_Number;

    begin
      col := col + 1;
      if col <= integer32(n) then
        mat(row,col) := Eval(t.dg,one,v(row).all);
      else
        eva := Eval(t,v(row).all);
        Sub(rhs(row),eva);
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    for i in 1..n loop
      row := integer32(i); col := 0;
      Scan_Terms(p);
    end loop;
    Clear(one);
  end Coefficient_System1;

  procedure Coefficient_System
              ( p : in Poly; n : in natural32; v : in VecVec;
                mat : in out Multprec_Complex_Matrices.Matrix;
                rhs : in out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Builds the matrix and right-hand-side vector of the coefficient
  --   system for the polynomial p, interpolating at the points in v.
  --   The dimension of the system is n.

    row,col : integer32;
    one : Complex_Number := Create(integer(1));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      eva : Complex_Number;

    begin
      col := col+1;
      if col <= integer32(n) then
        mat(row,col) := Eval(t.dg,one,v(row-1).all);
      else
        eva := Eval(t,v(row-1).all);
        Sub(rhs(row),eva);
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    for j in 1..integer32(n) loop       -- first equation is sum of coeff = 1
      Copy(one,mat(1,j));
    end loop;
    Copy(one,rhs(1));
    for i in 2..integer32(n) loop       -- run over the other equations
      row := i; col := 0;
      rhs(i) := Create(integer(0));
      Scan_Terms(p);
    end loop;
    Clear(one);
  end Coefficient_System;

  function Assign_Coefficients
             ( p : Poly; c : Multprec_Complex_Vectors.Vector ) return Poly is

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
        Copy(c(ind),nt.cf);
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
      rancmp : Standard_Complex_Numbers.Complex_Number;

    begin
      if k > n then
        if cff = 1 then
          accu.cf := Create(integer(1));
        elsif cff = 1 then
          ranflt := Random;
          rancmp := Standard_Complex_Numbers.Create(ranflt);
          accu.cf := Create(rancmp);
        else
          accu.cf := Create(Random1);
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

  function Sample ( p : in Poly; m,sz : natural32 ) return VecVec is

    res : VecVec(1..integer32(m));
    d : constant integer32 := Degree(p);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    point : Multprec_Complex_Vectors.Vector(1..n);
    p1 : Multprec_Complex_Vectors.Vector(0..d);

  begin
    for i in 1..m loop
      point := Random_Vector(1,n,sz);
      p1 := Evaluate(p,point(1..n-1));
      Newton(p1,point(n));
      res(integer32(i)) := new Multprec_Complex_Vectors.Vector'(point);
    end loop;
    return res;
  end Sample;

  procedure Interpolate1 ( p : in Poly; v : in VecVec; ip : out Poly;
                           invcnd : out Floating_Number ) is

  -- NOTE : this interpolate takes the constant term equal to one.

    n : constant integer32
      := integer32(Minimum(Number_of_Terms(p)-1,natural32(v'length)));
    mat : Multprec_Complex_Matrices.Matrix(1..n,1..n);
    rhs : Multprec_Complex_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
   -- qraux : Standard_Complex_Vectors.Vector(1..n) := (1..n => Create(0.0));
   -- dum : Standard_Complex_Vectors.Vector(1..n);
   -- info : integer;         -- needed for lufac and QRD

  begin
    for i in 1..n loop
      ipvt(i) := i;
      rhs(i) := Create(integer(0));
      for j in 1..n loop
        mat(i,j) := Create(integer(0));
      end loop;
    end loop;
    Coefficient_System1(p,natural32(n),v,mat,rhs);
   -- Gaussian elimination on the coefficient matrix :
   -- lufac(mat,n,ipvt,info);
   -- invcnd := Create(integer(1));
    lufco(mat,n,ipvt,invcnd);
    lusolve(mat,n,ipvt,rhs);
   -- QR decomposition followed by Least Squares :
   -- QRD(mat,qraux,ipvt,false);                         -- QR without pivoting
   -- QRLS(mat,n,n,n,qraux,rhs,dum,dum,rhs,dum,dum,100,info);  -- least squares
   -- Scale_Coefficients(rhs,maxcols);
    ip := Assign_Coefficients(p,rhs);
    Multprec_Complex_Matrices.Clear(mat);
    Multprec_Complex_Vectors.Clear(rhs);
  end Interpolate1;

  procedure Interpolate ( p : in Poly; v : in VecVec; ip : out Poly;
                          invcnd : out Floating_Number ) is

  -- NOTE : with this interpolate, all coefficients have to add up to one.

    n : constant integer32 := integer32(Number_of_Terms(p));
      --Minimum(Number_of_Terms(p),v'length);
    mat : Multprec_Complex_Matrices.Matrix(1..n,1..n);
    rhs : Multprec_Complex_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
   -- qraux : Standard_Complex_Vectors.Vector(1..n) := (1..n => Create(0.0));
   -- dum : Standard_Complex_Vectors.Vector(1..n);
   -- info : integer;         -- needed for lufac and QRD

  begin
    for i in 1..n loop
      ipvt(i) := i;
      rhs(i) := Create(integer(0));
      for j in 1..n loop
        mat(i,j) := Create(integer(0));
      end loop;
    end loop;
    Coefficient_System(p,natural32(n),v,mat,rhs);
   -- Gaussian elimination on the coefficient matrix :
   -- lufac(mat,n,ipvt,info);
   -- invcnd := Create(integer(1));
    lufco(mat,n,ipvt,invcnd);
    lusolve(mat,n,ipvt,rhs);
   -- QR decomposition followed by Least Squares :
   -- QRD(mat,qraux,ipvt,false);                         -- QR without pivoting
   -- QRLS(mat,n,n,n,qraux,rhs,dum,dum,rhs,dum,dum,100,info);  -- least squares
   -- Scale_Coefficients(rhs,maxcols);
    ip := Assign_Coefficients(p,rhs);
    Multprec_Complex_Matrices.Clear(mat);
    Multprec_Complex_Vectors.Clear(rhs);
  end Interpolate;

  function Interpolate ( p : Poly; v : VecVec ) return Poly is

    res : Poly;
    n : constant integer32 := integer32(Number_of_Terms(p));
      --Minimum(Number_of_Terms(p),v'length);
    mat : Multprec_Complex_Matrices.Matrix(1..n,1..n);
    rhs : Multprec_Complex_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    qraux : Multprec_Complex_Vectors.Vector(1..n) := (1..n => Create(integer(0)));
    dum : Multprec_Complex_Vectors.Vector(1..n);
    info : integer32;         -- needed for lufac and QRD

  begin
    for i in 1..n loop
      ipvt(i) := i;
      rhs(i) := Create(integer(0));
      for j in 1..n loop
        mat(i,j) := Create(integer(0));
      end loop;
    end loop;
    Coefficient_System(p,natural32(n),v,mat,rhs);
   -- Gaussian elimination on the coefficient matrix :
   -- lufac(mat,n,ipvt,info);
   -- invcnd := Create(integer(1));
   -- lufco(mat,n,ipvt,invcnd);
   -- lusolve(mat,n,ipvt,rhs);
   -- QR decomposition followed by Least Squares :
    QRD(mat,qraux,ipvt,false);                         -- QR without pivoting
    QRLS(mat,n,n,n,qraux,rhs,dum,dum,rhs,dum,dum,100,info);  -- least squares
   -- Scale_Coefficients(rhs,maxcols);
    res := Assign_Coefficients(p,rhs);
    Multprec_Complex_Matrices.Clear(mat);
    Multprec_Complex_Vectors.Clear(rhs);
    return res;
  end Interpolate;

  function Norm ( p : Poly ) return Floating_Number is

    cff : constant Multprec_Complex_Vectors.Vector := Coeff(p);

  begin
    return Max_Norm(cff);
  end Norm;

  function Distance ( p,q : Poly ) return Floating_Number is

    cp : constant Multprec_Complex_Vectors.Vector := Coeff(p);
    cq : constant Multprec_Complex_Vectors.Vector := Coeff(q);
    use Multprec_Complex_Vectors;

  begin
    if cp'last > cq'last
     then return Max_Norm(cp(cq'range) - cq);
     else return Max_Norm(cp - cq(cp'range));
    end if;
  end Distance;

end Multprec_Polynomial_Interpolators;
