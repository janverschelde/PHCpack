with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Multprec_Complex_Poly_Functions;   use Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_SysFun;      use Multprec_Complex_Poly_SysFun;
with Monomial_Hashing;                  use Monomial_Hashing;
 
package body Multprec_Nullity_Polynomials is

  function Derivative 
              ( p : Poly; m : Standard_Natural_Vectors.Vector ) return Poly is

    res : Poly;

  begin
    Copy(p,res);
    for i in m'range loop
      if m(i) > 0 then
        for j in 1..m(i) loop
          Diff(res,i);
          exit when (res = Null_Poly);
        end loop; 
      end if;
    end loop;
    return res;
  end Derivative;

  function Factorial ( m : Standard_Natural_Vectors.Vector )
                     return natural32 is

    res : natural32 := 1;

  begin
    for i in m'range loop
      if m(i) > 1 then
        for j in 2..m(i) loop
          res := res*j;
        end loop;
      end if;
    end loop;
    return res;
  end Factorial;

  function Monomial_Multiple
              ( m : Standard_Natural_Vectors.Vector;
                f : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(f'range);
    t : Term;

  begin
    t.cf := Create(integer32(1));
    t.dg := new Standard_Natural_Vectors.Vector'(m);
    for i in f'range loop
      res(i) := t*f(i);
    end loop;
    Clear(t);
    return res;
  end Monomial_Multiple;

  procedure Evaluate_Derivatives
              ( a : in out Multprec_Complex_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in Multprec_Complex_Vectors.Vector ) is

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is

      dp : Poly;
      row : natural32 := r;
      fac : constant natural32 := Factorial(m);
      mpf_fac : Floating_Number := Create(fac);

    begin
      for i in 1..integer32(nq) loop
        dp := Derivative(f(i),m);
        if dp = Null_Poly then
          a(integer32(row),integer32(c)) := Create(integer32(0));
        else
          a(integer32(row),integer32(c)) := Eval(dp,z);
          Div(a(integer32(row),integer32(c)),mpf_fac);
          Clear(dp);
        end if;
        row := row + 1;
      end loop;       
      Clear(mpf_fac);
      c := c + 1;
    end Monomial;
    procedure Enumerate is new Enumerate_Monomials(Monomial);

  begin
    Enumerate(k,nv);
  end Evaluate_Derivatives;

  procedure Evaluate_Derivatives
              ( file : in file_type;
                a : in out Multprec_Complex_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in Multprec_Complex_Vectors.Vector ) is

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is

      dp : Poly;
      row : natural32 := r;
      fac : constant natural32 := Factorial(m);
      mpf_fac : Floating_Number := Create(fac);

    begin
      for i in 1..integer32(nq) loop
        dp := Derivative(f(i),m);
        if dp = Null_Poly then
          a(integer32(row),integer32(c)) := Create(integer32(0));
        else
          a(integer32(row),integer32(c)) := Eval(dp,z);
          Div(a(integer32(row),integer32(c)),mpf_fac);
          Clear(dp);
        end if;
        put(file,"row = "); put(file,row,1);
        put(file,"  column = "); put(file,c,1); new_line(file);
        row := row + 1;
      end loop;
      Clear(mpf_fac);
      c := c + 1;
    end Monomial;
    procedure Enumerate is new Enumerate_Monomials(Monomial);

  begin
    Enumerate(k,nv);
  end Evaluate_Derivatives;

  procedure Compute_Derivatives
              ( a : in out Multprec_Complex_Poly_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys ) is

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is

      dp : Poly;
      row : natural32 := r;
      fac : constant natural32 := Factorial(m);
      mpf_fac : Floating_Number := Create(fac);
      mpc_fac : Complex_Number := Create(1.0/mpf_fac);

    begin
      for i in 1..integer32(nq) loop
        dp := Derivative(f(i),m);
        if dp = Null_Poly then
          a(integer32(row),integer32(c)) := Null_Poly;
        else
          Mul(dp,mpc_fac);
          a(integer32(row),integer32(c)) := dp;
        end if;
        row := row + 1;
      end loop;
      Clear(mpf_fac);
      Clear(mpc_fac);
      c := c + 1;
    end Monomial;
    procedure Enumerate is new Enumerate_Monomials(Monomial);

  begin
    Enumerate(k,nv);
  end Compute_Derivatives;

  procedure Compute_Derivatives
              ( file : in file_type;
                a : in out Multprec_Complex_Poly_Matrices.Matrix;
                r : in natural32; c : in out natural32; nq,nv,k : in natural32;
                f : in Poly_Sys ) is

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is

      dp : Poly;
      row : natural32 := r;
      fac : constant natural32 := Factorial(m);
      mpf_fac : Floating_Number := Create(fac);
      mpc_fac : Complex_Number := Create(1.0/mpf_fac);

    begin
      for i in 1..integer32(nq) loop
        dp := Derivative(f(i),m);
        if dp = Null_Poly then
          a(integer32(row),integer32(c)) := Null_Poly;
        else
          Mul(dp,mpc_fac);
          a(integer32(row),integer32(c)) := dp;
        end if;
        put(file,"row = "); put(file,row,1);
        put(file,"  column = "); put(file,c,1); new_line(file);
        row := row + 1;
      end loop;
      Clear(mpf_fac);
      Clear(mpc_fac);
      c := c + 1;
    end Monomial;
    procedure Enumerate is new Enumerate_Monomials(Monomial);

  begin
    Enumerate(k,nv);
  end Compute_Derivatives;

  procedure Evaluate_All_Derivatives
              ( a : in out Multprec_Complex_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in Multprec_Complex_Vectors.Vector ) is

    wr : natural32 := r;

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is
 
    -- DESCRIPTION :
    --   Evaluates all derivatives up to order k of m*f.

      mf : Poly_Sys(f'range) := Monomial_Multiple(m,f);
      row : natural32 := wr;
      wc : natural32;
 
    begin
      for i in 1..integer32(nq) loop
        a(integer32(row),integer32(c)) := Eval(mf(i),z);
        row := row + 1;
      end loop;
      wc := c+1;
      for i in 1..k loop
        Evaluate_Derivatives(a,wr,wc,nq,nv,i,mf,z);
      end loop;
      wr := wr + nq;
      Clear(mf);
    end Monomial;
    procedure Enumerate is new Enumerate_Monomials(Monomial);

  begin
    Enumerate(k-1,nv);
  end Evaluate_All_Derivatives;

  procedure Compute_All_Derivatives
              ( a : in out Multprec_Complex_Poly_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32;
                f : in Poly_Sys ) is

    wr : natural32 := r;

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is
 
    -- DESCRIPTION :
    --   Computes all derivatives up to order k of m*f.

      mf : constant Poly_Sys(f'range) := Monomial_Multiple(m,f);
      row : natural32 := wr;
      wc : natural32;
 
    begin
      for i in 1..integer32(nq) loop
        a(integer32(row),integer32(c)) := mf(i);
        row := row + 1;
      end loop;
      wc := c+1;
      for i in 1..k loop
        Compute_Derivatives(a,wr,wc,nq,nv,i,mf);
      end loop;
      wr := wr + nq;
    end Monomial;
    procedure Enumerate is new Enumerate_Monomials(Monomial);

  begin
    Enumerate(k-1,nv);
  end Compute_All_Derivatives;

  procedure Evaluate_Highest_Order
              ( a : in out Multprec_Complex_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32;
                f : in Poly_Sys; z : in Multprec_Complex_Vectors.Vector ) is

    wr : natural32 := r;

    procedure Highest_Order ( m : in Standard_Natural_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Evaluates the highest order derivative of m*f at z.

      mf : constant Poly_Sys(f'range) := Monomial_Multiple(m,f);
      wc : natural32 := c;

    begin
      Evaluate_Derivatives(a,wr,wc,nq,nv,k,mf,z);
      wr := wr+nq;
    end Highest_Order;
    procedure Enumerate is new Enumerate_Monomials(Highest_Order);

  begin
    for i in 1..k-1 loop
      Enumerate(i,nv);
    end loop;
  end Evaluate_Highest_Order;

  procedure Compute_Highest_Order
              ( a : in out Multprec_Complex_Poly_Matrices.Matrix;
                r,c : in natural32; nq,nv,k : in natural32;
                f : in Poly_Sys ) is

    wr : natural32 := r;

    procedure Highest_Order ( m : in Standard_Natural_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Computes the highest order derivative of m*f at z.

      mf : constant Poly_Sys(f'range) := Monomial_Multiple(m,f);
      wc : natural32 := c;

    begin
      Compute_Derivatives(a,wr,wc,nq,nv,k,mf);
      wr := wr+nq;
    end Highest_Order;
    procedure Enumerate is new Enumerate_Monomials(Highest_Order);

  begin
    for i in 1..k-1 loop
      Enumerate(i,nv);
    end loop;
  end Compute_Highest_Order;

  procedure Evaluate_Monomial_Multiples
              ( a : in out Multprec_Complex_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32;
                f : in Poly_Sys; z : in Multprec_Complex_Vectors.Vector ) is

    wr,wc : natural32;

  begin
    wc := nc1+1;
    Evaluate_Derivatives(a,1,wc,nq,nv,k,f,z);
    wr := nq+1;
    wc := nc1+1;
    Evaluate_Highest_Order(a,wr,wc,nq,nv,k,f,z);
    wr := r;
    Evaluate_All_Derivatives(a,r,1,nq,nv,k,f,z);
  end Evaluate_Monomial_Multiples;

  procedure Evaluate_Monomial_Multiples
              ( file : in file_type;
                a : in out Multprec_Complex_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32;
                f : in Poly_Sys; z : in Multprec_Complex_Vectors.Vector ) is

    wr,wc : natural32;

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is
 
    -- DESCRIPTION :
    --   Evaluates all derivatives up to order k of m*f.

      mf : Poly_Sys(f'range) := Monomial_Multiple(m,f);
      row : natural32 := wr;
 
    begin
      put(file,"multiply f with monomial "); put(file,m); new_line(file);
      for i in 1..integer32(nq) loop
        put(file,"row = "); put(file,row,1);
        put(file,"  column = "); put(file,c,1); new_line(file);
        a(integer32(row),integer32(c)) := Eval(mf(i),z);
        row := row + 1;
      end loop;
      wc := c+1;
      for i in 1..k loop
       -- Evaluate_Derivatives(file,a,wr,wc,nq,nv,i,mf,z);
        Evaluate_Derivatives(a,wr,wc,nq,nv,i,mf,z);
      end loop;
      wr := wr + nq;
      Clear(mf);
    end Monomial;
    procedure Evaluate_All_Derivatives is new Enumerate_Monomials(Monomial);

    procedure Highest_Order ( m : in Standard_Natural_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Evaluates the highest order derivative of m*f at z.

      mf : constant Poly_Sys(f'range) := Monomial_Multiple(m,f);

    begin
      put(file,"multiply f with monomial "); put(file,m); new_line(file);
      wc := nc1+1;
     -- Evaluate_Derivatives(file,a,wr,wc,nq,nv,k,mf,z);
      Evaluate_Derivatives(a,wr,wc,nq,nv,k,mf,z);
      wr := wr+nq;
    end Highest_Order;
    procedure Evaluate_Highest_Order is new Enumerate_Monomials(Highest_Order);

  begin
    wc := nc1+1;
    put_line(file,"evaluating highest order for i = 0");
   -- Evaluate_Derivatives(file,a,1,wc,nq,nv,k,f,z);
    Evaluate_Derivatives(a,1,wc,nq,nv,k,f,z);
    wr := nq+1;
    for i in 1..k-1 loop
      put(file,"evaluating highest order for i = ");
      put(file,i,1); new_line(file);
      Evaluate_Highest_Order(i,nv);
    end loop;
    wr := r;
    Evaluate_All_Derivatives(k-1,nv);
  end Evaluate_Monomial_Multiples;

  procedure Compute_Monomial_Multiples
              ( a : in out Multprec_Complex_Poly_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32; f : in Poly_Sys ) is

    wr,wc : natural32;

  begin
    wc := nc1+1;
    Compute_Derivatives(a,1,wc,nq,nv,k,f);
    wr := nq+1;
    wc := nc1+1;
    Compute_Highest_Order(a,wr,wc,nq,nv,k,f);
    wr := r;
    Compute_All_Derivatives(a,r,1,nq,nv,k,f);
  end Compute_Monomial_Multiples;

  procedure Compute_Monomial_Multiples
              ( file : in file_type;
                a : in out Multprec_Complex_Poly_Matrices.Matrix;
                r,c,nq,nv,k,nc1 : in natural32; f : in Poly_Sys ) is

    wr,wc : natural32;

    procedure Monomial ( m : in Standard_Natural_Vectors.Vector ) is
 
    -- DESCRIPTION :
    --   Computes all derivatives up to order k of m*f.

      mf : constant Poly_Sys(f'range) := Monomial_Multiple(m,f);
      row : natural32 := wr;
 
    begin
      put(file,"multiply f with monomial "); put(file,m); new_line(file);
      for i in 1..integer32(nq) loop
        put(file,"row = "); put(file,row,1);
        put(file,"  column = "); put(file,c,1); new_line(file);
        a(integer32(row),integer32(c)) := mf(i);
        row := row + 1;
      end loop;
      wc := c+1;
      for i in 1..k loop
       -- Compute_Derivatives(file,a,wr,wc,nq,nv,i,mf);
        Compute_Derivatives(a,wr,wc,nq,nv,i,mf);
      end loop;
      wr := wr + nq;
    end Monomial;
    procedure Compute_All_Derivatives is new Enumerate_Monomials(Monomial);

    procedure Highest_Order ( m : in Standard_Natural_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Computes the highest order derivative of m*f.

      mf : constant Poly_Sys(f'range) := Monomial_Multiple(m,f);

    begin
      put(file,"multiply f with monomial "); put(file,m); new_line(file);
      wc := nc1+1;
     -- Compute_Derivatives(file,a,wr,wc,nq,nv,k,mf);
      Compute_Derivatives(a,wr,wc,nq,nv,k,mf);
      wr := wr+nq;
    end Highest_Order;
    procedure Compute_Highest_Order is new Enumerate_Monomials(Highest_Order);

  begin
    wc := nc1+1;
    put_line(file,"computing highest order for i = 0");
   -- Compute_Derivatives(file,a,1,wc,nq,nv,k,f);
    Compute_Derivatives(a,1,wc,nq,nv,k,f);
    wr := nq+1;
    for i in 1..k-1 loop
      put(file,"computing highest order for i = ");
      put(file,i,1); new_line(file);
      Compute_Highest_Order(i,nv);
    end loop;
    wr := r;
    Compute_All_Derivatives(k-1,nv);
  end Compute_Monomial_Multiples;

end Multprec_Nullity_Polynomials;
