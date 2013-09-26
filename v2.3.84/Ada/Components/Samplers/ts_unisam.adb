with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;

procedure ts_unisam is

-- DESCRIPTION :
--   Interactive environment to develop sampling of hypersurfaces.
--   The idea is to apply the method of Durand-Kerner, either in
--   isolation, or in a homotopy continuation setting.

-- THE METHOD OF WEIERSTRASS (aka Durand-Kerner) :

  function Compute_q ( i : integer32; a : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   Computes the quotient needed in the Durand-Kerner step.

    res : Complex_Number;

  begin
    res := Create(1.0);
    for j in a'range loop
      if j /= i
       then res := res*(a(i)-a(j));
      end if;
    end loop;
    return res;
  end Compute_q;

  function Eval ( p : Eval_Poly; basis,direc : Vector; t : Complex_Number )
                return Complex_Number is

  -- DESCRIPTION :
  --   Returns the function value of p at the point basis + t*direc.

    point : constant Vector(basis'range) := basis + t*direc;

  begin
    return Eval(p,point);
  end Eval;

  procedure DK ( p : in Eval_Poly; basis,direc : in Vector;
                 lc : Complex_Number; z,res : in out Vector ) is

  -- DESCRIPTION :
  --   Computes one step in the Durand-Kerner iteration with
  --   the polynomial p/lc, where lc is the leading coefficient
  --   for the parameteric representation using direc.

  begin
    for i in z'range loop
      declare
      begin
        z(i) := z(i) - Eval(p,basis,direc,z(i))/(lc*Compute_q(i,z));
        res(i) := Eval(p,basis,direc,z(i))/lc;
      exception
        when others => --put("Exception occurred at component ");
                       --put(i,1); put_line(".");
                       z(i) := Random1;
      end;
    end loop;
  end DK;

  procedure Reporting_Durand_Kerner
              ( file : in file_type;
                p : in Eval_Poly; basis,direc : in Vector;
                lc : in Complex_Number; z,res : in out Vector;
                maxsteps : in natural32; eps : in double_float;
                nb : out natural32; fail : out boolean ) is

    nrm,absres : double_float;

  begin
    fail := true;
    for k in 1..maxsteps loop
      DK(p,basis,direc,lc,z,res);
      nrm := Max_Norm(res);
      put(file,"Max norm at step "); put(file,k,1);
      put(file," : "); put(file,nrm); new_line(file);
      for i in z'range loop
        put(file,z(i)); put(file," : ");
        absres := AbsVal(res(i));
        put(file,absres); new_line(file);
      end loop;
      if nrm <= eps then
        fail := false;
        nb := k;
        return;
      end if;
    end loop;
    nb := maxsteps;
  end Reporting_Durand_Kerner;

-- SETTING UP THE SYSTEM :

  function Determine_Slices ( n : integer32 ) return VecVec is

  -- DESCRIPTION :
  --   Returns n-1 hyperplanes given by the user or randomly generated.

    res : VecVec(1..n-1);
    cff : Vector(0..n);
    ans : character;

  begin
    new_line;
    put("MENU to generate "); put(n-1,1); put_line(" general slices :");
    put_line("  1. give complex coefficients of the slices; or");
    put_line("  2. let random number generator determine the coefficients.");
    put("Type 1 or 2 to select : "); Ask_Alternative(ans,"12");
    if ans = '1' then
      for j in 1..n-1 loop
        put("Reading the coefficients of slice "); put(j,1);
        put_line(" ...");
        for i in 0..n loop
          put("  cff("); put(j,1); put(",");
          put(i,1); put(") : "); get(cff(i));
        end loop;
        res(j) := new Vector'(cff);
      end loop;
    else
      for j in 1..n-1 loop
        res(j) := new Vector'(Random_Vector(0,n));
      end loop;
    end if;
    put_line("The coefficients of the slicing hyperplanes : ");
    for i in res'range loop
      put("Coefficients of slice "); put(i,1); put_line(" :");
      put_line(res(i));
    end loop;
    return res;
  end Determine_Slices;

  procedure Parametric_Representation
              ( n : in integer32; hyps : in VecVec;
                basis : out Vector; direc : out Vector ) is

  -- DESCRIPTION :
  --   Returns the parametric representation of the line defined
  --   by the n-1 hyperplanes in n-space.

  -- REQUIRED : the normals to the hyperplanes span n-1 space.

  -- ON ENTRY :
  --   n        dimension of the ambient space;
  --   hyps     n-1 vectors of range 0..n representing coefficients
  --            of the hyperplanes c(0) + c(1)*x(1) + .. + c(n)*x(n) = 0.

  -- ON RETURN :
  --   basis    point that satisfies all hyperplanes;
  --   direc    direction of the line.

  -- ALGORITHM : add twice a random slice to it and solve the system.
  --   The first time we find the basis vector, next time we make the
  --   difference with first and second solution to get the direction.

    mat,wrk : Matrix(1..n,1..n);
    rhs : Vector(1..n);
    pvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    for i in hyps'range loop
      for j in 1..n loop
        mat(i,j) := hyps(i)(j);
      end loop;
      rhs(i) := -hyps(i)(0);
    end loop;
    for j in 1..n loop
      mat(n,j) := Random1;
    end loop;
    rhs(n) := Random1;
    wrk := mat; 
    lufac(wrk,n,pvt,info);
    basis := rhs;
    lusolve(wrk,n,pvt,basis);
    for j in 1..n loop
      mat(n,j) := Random1;
    end loop;
    rhs(n) := Random1; wrk := mat;
    lufac(wrk,n,pvt,info);
    direc := rhs;
    lusolve(wrk,n,pvt,direc);
    for i in direc'range loop
      direc(i) := direc(i) - basis(i);
    end loop;
  end Parametric_Representation;

  function Lead_Coefficient ( t : Term; direc : Vector )
                            return Complex_Number is

  -- DESCRIPTION :
  --   For x(i) = a + l*direc(i), with l some parameter, the function
  --   return the coefficients of the highest degree term in l when
  --   expanding the x(i) in terms of a + l*direc(i).

    res : Complex_Number := t.cf;

  begin
    for i in t.dg'range loop
      for j in 1..t.dg(i) loop
        res := res*direc(i);
      end loop;
    end loop;
    return res;
  end Lead_Coefficient;

  function Lead_Coefficient ( p : Poly; d : natural32; direc : Vector )
                            return Complex_Number is

  -- DESCRIPTION :
  --   Applies Lead_Coefficient to every term of degree d in the
  --   polynomials and returns the sum.

    res : Complex_Number := Create(0.0);

    procedure Lead_by_Term ( t : in Term; cont : out boolean ) is

    begin
      if Standard_Natural_Vectors.Sum(t.dg.all) = d
       then res := res + Lead_Coefficient(t,direc);
      end if;
      cont := true;
    end Lead_by_Term;
    procedure Lead_by_Terms is new Visiting_Iterator(Lead_by_Term);

  begin
    Lead_by_Terms(p);
    return res;
  end Lead_Coefficient;

  procedure Slice ( n : integer32; p : in Poly ) is

  -- DESCRIPTION :
  --   This procedure reads the coefficients for a hyperplane and
  --   uses the parametric representation of the plane to sample
  --   points from the hypersurface defined by the polynomial.

    basis,direc,point : Vector(1..n);
    val,t : Complex_Number;
    ans : character;
    hyps : VecVec(1..n-1) := Determine_Slices(n);
    ep : Eval_Poly;
    z,res : Vector(1..Degree(p));
    max,nit : natural32 := 0;
    eps : double_float := 1.0E-8;
    fail : boolean;
    lc : Complex_Number;

  begin
    Parametric_Representation(n,hyps,basis,direc);
    put_line("Basis vector :"); put_line(basis);
    put_line("Its value at the slices : ");
    for i in 1..n-1 loop
      val := hyps(i)(0) + hyps(i)(basis'range)*basis;
      put(val); new_line;
    end loop;
    loop
      t := Random1;
      point := basis + t*direc;
      put_line("A random point on the slice : ");
      put_line(point);
      put_line("Its value at the slices : ");
      for i in 1..n-1 loop
        val := hyps(i)(0) + hyps(i)(basis'range)*basis;
        put(val); new_line;
      end loop;
      put("Do you want more samples ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
    put_line("Finding the roots on the slices ...");
    put("Give the maximal number of steps : "); get(max);
    put("Give the desired accuracy : "); get(eps);
    ep := Create(p);
    z := Random_Vector(1,z'last);
    res := z;
    lc := Lead_Coefficient(p,natural32(Degree(p)),direc);
    Reporting_Durand_Kerner
      (Standard_Output,ep,basis,direc,lc,z,res,max,eps,nit,fail);
    Standard_Complex_VecVecs.Clear(hyps);
  end Slice;

  procedure Main is

    n : natural32 := 0;
    p : Poly;
 
  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give your polynomial : "); get(p);
    put("-> Your polynomial : "); put(p); new_line;
    Slice(integer32(n),p);
  end Main;

begin
  new_line;
  put_line("Sampling a multivariate polynomial with complex coefficients.");
  Main;
end ts_unisam;
