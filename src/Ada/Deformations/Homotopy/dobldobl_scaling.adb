with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;

package body DoblDobl_Scaling is

-- SIMPLE EQUATION SCALING :

  procedure Scale ( p : in out Poly ) is

    sum : double_double := create(0.0);
    number : natural32 := 0;
    dd_fac : double_double;
    factor : Complex_Number;

    procedure Add_To_Sum ( t : in Term; continue : out boolean ) is
    begin
      sum := sum + AbsVal(t.cf);
      number := number + 1;
      continue := true;
    end Add_To_Sum;
    procedure Compute_Sum is new Visiting_Iterator(Add_To_Sum);

  begin
    Compute_Sum(p);
    dd_fac := double_double_numbers.create(number)/sum;
    factor := Create(dd_fac);
    Mul(p,factor);
  end Scale;

  procedure Scale ( s : in out Poly_Sys ) is
  begin
    for i in s'range loop
      scale(s(i));
    end loop;
  end Scale;

-- AUXILIARIES TO IMPLEMENT VARIABLE SCALING :
 
  function Center_Coefficients
             ( s : Poly_Sys; n,m,npm : integer32; bas : natural32 )
             return Poly is

  -- DESCRIPTION :
  --   Returns the polynomial to be minimized for centering the coefficients
  --   about unity, which is denoted by "r1" in Alexander Morgan's book. 

  -- ON ENTRY :
  --   s       polynomial system in need of scaling;
  --   n       number of equations;
  --   m       number of unknowns;
  --   npm     is n+m;
  --   bas     basis of the logarithm.

    r1,r1i : Poly;

    procedure Scan ( p : in Poly; i : in integer32 ) is

      init : Poly;
      t_init : Term;
   
      procedure Scan_Term ( t : in Term; continue : out boolean ) is

        tt : Term;
        temp : Poly;

      begin
        Copy(init,temp);
        continue := true;
        for k in t.dg'range loop
          if t.dg(k) /= 0 then
            tt.cf := Create(double_double_numbers.create(t.dg(k)));
            tt.dg := new Standard_Natural_Vectors.Vector'(1..npm => 0);
            tt.dg(k) := 1;
            Add(temp,tt);
            Clear(tt);
          end if;
        end loop;
        tt.cf := Create(LOG10(AbsVal(t.cf)));
        tt.dg := new Standard_Natural_Vectors.Vector'(1..npm => 0);
        Add(temp,tt);
        Clear(tt);
        Mul(temp,temp);
        Add(r1i,temp);
        Clear(temp);
      end Scan_Term;
      procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

    begin
      t_init.cf := Create(double_double_numbers.create(1.0));
      t_init.dg := new Standard_Natural_Vectors.Vector'(1..npm => 0);
      t_init.dg(m+i) := 1;
      init := Create(t_init);
      Clear(t_init);
      Scan_Terms(p);
    end Scan;

  begin
    r1 := Null_Poly;
    for i in s'range loop
      Scan(s(i),i);
      Add(r1,r1i);
      Clear(r1i);
    end loop;
    return r1;
  end Center_Coefficients;

  function Reduce_Variability
             ( s : Poly_Sys; n,m,npm : integer32; bas : natural32 )
             return Poly is

  -- DESCRIPTION :
  --   Returns the polynomial to be minimized for reducing the variability of
  --   the coefficients, which is denoted by "r2" in Alexander Morgan's book. 

  -- ON ENTRY :
  --   s       polynomial system in need of scaling;
  --   n       number of equations;
  --   m       number of unknowns;
  --   npm     is n+m;
  --   bas     basis of the logarithm.

    r2,r2i : Poly;

    procedure Scan2 ( p : in Poly; t : in Term; nr : in natural32 ) is

      count : natural32 := 0;

      procedure Scan2_Term ( t2 : in Term; continue : out boolean ) is

        tt : Term;
        temp : Poly := Null_Poly;

      begin
        continue := true;
        count := count + 1;
        if count > nr then
          for i in t2.dg'range loop
            if t.dg(i) /= t2.dg(i) then
              tt.dg := new Standard_Natural_Vectors.Vector'(1..npm => 0);
              tt.dg(i) := 1;
              tt.cf := Create(double_double_numbers.create(t.dg(i)-t2.dg(i)));
              Add(temp,tt);
              Clear(tt);
            end if;
          end loop;
        end if;
        tt.dg := new Standard_Natural_Vectors.Vector'(1..npm => 0);
        tt.cf := Create(LOG10(AbsVal(t.cf)/AbsVal(t2.cf)));
        Add(temp,tt);
        Clear(tt);
        Mul(temp,temp);
        Add(r2i,temp);
        Clear(temp);
      end Scan2_Term;
      procedure Scan2_Terms is new Visiting_Iterator(Scan2_Term);

    begin
      Scan2_Terms(p);
    end Scan2;

    procedure Scan1 ( p : in Poly ) is

      nr : natural32 := 0;

      procedure Scan_Term ( t : in Term; continue : out boolean ) is
      begin
        nr := nr + 1;
        continue := true;
        Scan2(p,t,nr);
      end Scan_Term;
      procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

    begin
      Scan_Terms(p);
    end Scan1;

  begin
    r2 := Null_Poly;
    for i in s'range loop
      Scan1(s(i));
      Add(r2,r2i);
      Clear(r2i);
    end loop;
    return r2;
  end Reduce_Variability;

  procedure Init_Linear_System ( m : out Matrix; r : out Vector ) is

  -- DESCRIPTION :
  --   Initializes the matrix and right-hand-side vector of a linear
  --   system to zero.

  begin
    for i in m'range(1) loop
      r(i) := Create(double_double_numbers.create(0.0));
      for j in m'range(2) loop
        m(i,j) := Create(double_double_numbers.create(0.0));
      end loop;
    end loop;
  end Init_Linear_System;

  procedure Make_Linear_System ( r : in Poly; mat : out Matrix;
                                 right : out Vector ) is

    drj : Poly;

    procedure Scan ( p : in Poly; j : in integer32 ) is

      procedure Scan_Term ( t : in Term; continue : out boolean ) is
      begin
        continue := true;
        for i in t.dg'range loop
          if t.dg(i) = 1
           then mat(j,i) := t.cf; return;
          end if;
        end loop;
        right(j) := -t.cf;
      end Scan_Term;
      procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

    begin
      Scan_Terms(p);
    end Scan;

  begin
    Init_Linear_System(mat,right);
    for j in 1..integer32(Number_Of_Unknowns(r)) loop
      drj := DoblDobl_Complex_Polynomials.Diff(r,j);
      Scan(drj,j);
    end loop;
  end Make_Linear_System;

  procedure Scale ( s : in out Poly_Sys; bas : in natural32 := 10;
                    diff : in boolean; cond : out double_double;
                    sccff : out vector ) is
   
    r1,r2,target : Poly;
    n : constant integer32 := s'length;
    m : constant integer32 := integer32(Number_of_Unknowns(s(s'first)));
    npm : constant integer32 := n+m;
    mat : Matrix(1..npm,1..npm);
    right : Vector(1..npm);

    procedure Scale ( s : in out Poly_Sys; bas : in natural32;
                      mat : in out Matrix; right : in out Vector;
                      cond : out double_double ) is
 
      ipvt : Standard_Integer_Vectors.Vector(mat'range(2));
 
      procedure Scale ( p : in out Poly; ip : in integer32;
                        scalingcoeff : in Vector; bas : in natural32 ) is

        procedure Scale_Term ( t : in out Term; continue : out boolean ) is

          exp : double_double := create(0.0);
          ddb : double_double := create(bas);
          tdgk : double_double;

        begin
          exp := REAL_PART(scalingcoeff(m+ip));
          for k in t.dg'range loop
            tdgk := create(t.dg(k));
            exp := exp + tdgk*REAL_PART(scalingcoeff(k));
          end loop;
          t.cf := t.cf * Create(ddb**exp);
          continue := true;
        end Scale_Term;
        procedure Scale_Terms is new Changing_Iterator(Scale_Term);

      begin
        Scale_Terms(p);
      end Scale;
 
    begin
      lufco(mat,npm,ipvt,cond);
      lusolve(mat,npm,ipvt,right);
      for i in s'range loop
        Scale(s(i),i,right,bas);
      end loop;
    end Scale;
        
  begin
    r1 := Center_Coefficients(s,n,m,npm,bas);
    if diff then
      r2 := Reduce_Variability(s,n,m,npm,bas);
      target := r1 + r2;
      Clear(r1); Clear(r2);
    else
      Copy(r1,target); Clear(r1);
    end if;
    Make_Linear_System(target,mat,right);
    Clear(target);
    Scale(s,bas,mat,right,cond);
    sccff := right;
  end Scale;

  procedure Scale ( basis : in natural32; sccff : in Vector;
                    s : in out Solution ) is

    dd_basis : double_double := create(basis);
 
  begin
    for i in s.v'range loop
      s.v(i) := Create(dd_basis**REAL_PART(sccff(i))) * s.v(i);
    end loop;
  end Scale;

  procedure Scale ( basis : in natural32; sccff : in Vector;
                    sols : in out Solution_List ) is
  begin
    if not Is_Null(sols) then
      declare
        temp : Solution_List := sols;
        n : integer32 := Head_Of(sols).n;
        s : Solution(n);
        l : Link_To_Solution;
      begin
        while not Is_Null(temp) loop
          l := Head_Of(temp);
          s := l.all;
          Scale(basis,sccff,s);
          Clear(l);
          l := new Solution'(s);
          Set_Head(temp,l);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Scale;

end DoblDobl_Scaling;
