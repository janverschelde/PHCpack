with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;      use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Jaco_Matrices;     use Multprec_Complex_Jaco_Matrices;

procedure ts_jaco is

-- DESCRIPTION :
--   This routine provides basic testing routines for Jacobian matrices.

  procedure Test_Standard_Creation is

  -- DESCRIPTION :
  --   Tests the creation of a standard Jacobian matrix.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    n : integer32;

  begin
    new_line;
    put_line("Testing the creation of standard Jacobian matrices.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    n := lp'last;
    declare
      jm : constant Standard_Complex_Jaco_Matrices.Jaco_Mat(1..n,1..n)
         := Create(lp.all);
    begin
      put_line("The Jacobian matrix : ");
      for i in jm'range(1) loop
        for j in jm'range(2) loop
          put(jm(i,j));
        end loop;
        new_line;
      end loop;
    end;
  end Test_Standard_Creation;

  function Expand2 ( jm : Standard_Complex_Jaco_Matrices.Jaco_Mat;
                     r1,r2,c1,c2 : integer32 )
                   return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :                 [ jm(r1,c1)  jm(r1,c2) ]
  --   Returns the determinant of  [                      ].
  --                               [ jm(r2,c1)  jm(r2,c2) ]

    res,acc : Standard_Complex_Polynomials.Poly;

  begin
    res := jm(r1,c1)*jm(r2,c2);
    acc := jm(r2,c1)*jm(r1,c2);
    Standard_Complex_Polynomials.Sub(res,acc);
    Standard_Complex_Polynomials.Clear(acc);
    return res;
  end Expand2;

  type Boolean_Array is array ( integer32 range <> ) of boolean;

  procedure Expandn ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                      k,n : in integer32; active : in out Boolean_Array;
                      dj : out Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Returns in dj the determinant of the Jacobian matrix starting at
  --   column k of jm.  If active(i), then the row has not been deleted.

    acc : Poly;
    minsgn : boolean := false;
    r1,r2 : integer32;

  begin
    if k = n-1 then
      r1 := 0;
      for i in active'range loop
        if active(i) then
          if r1 = 0
           then r1 := i;
           else r2 := i;
          end if;
        end if;
      end loop;
      dj := Expand2(jm,r1,r2,n-1,n);
    else
      dj := Null_Poly;
      for i in active'range loop
        if active(i) then
          if jm(i,k) /= Null_Poly then
            active(i) := false;
            Expandn(jm,k+1,n,active,acc);
            Mul(acc,jm(i,k));
            if minsgn
             then Sub(dj,acc);
             else Add(dj,acc);
            end if;
            Clear(acc);
            active(i) := true;
          end if;
          minsgn := not minsgn;
        end if;
      end loop;
    end if;
  end Expandn;

  function Expand ( jm : Standard_Complex_Jaco_Matrices.Jaco_Mat )
                  return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;
    n : constant integer32 := jm'length(1);

  begin
    if n = 2
     then res := Expand2(jm,jm'first(1),jm'first(1)+1,
                            jm'first(2),jm'first(2)+1);
     else declare
            actrows : Boolean_Array(1..n) := (1..n => true);
          begin
            Expandn(jm,1,n,actrows,res);
          end;
    end if;
    return res;
  end Expand;

  procedure Test_Expand_by_Eval 
               ( n : in integer32;
                 jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                 dj : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Generates a random vector and computes the determinant of the
  --   Jacobian once via the determinant of the evaluated Jacobi matrix,
  --   and once via the evaluation of the expanded determinant.
  --   The results of this evaluations should agree.

    x : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    m : Standard_Complex_Matrices.Matrix(1..n,1..n) := Eval(jm,x);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    y : Standard_Complex_Numbers.Complex_Number;

  begin
    put_line("Evaluation at random vector : ");
    lufac(m,n,ipvt,info);
    if info /= 0 then
      y := Create(0.0);
    else
      y := Create(1.0);
      for i in m'range(1) loop
        if ipvt(i) > i
         then y := -y;
        end if;
        y := y*m(i,i);
      end loop;
    end if;
    put("as det of matrix : "); put(y); new_line;
    y := Eval(dj,x);
    put("in expanded form : "); put(y); new_line;
  end Test_Expand_by_Eval;

  procedure Test_Standard_Expansion is

  -- DESCRIPTION :
  --   Tests the creation of a standard Jacobian matrix.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    dj : Standard_Complex_Polynomials.Poly;
    n : integer32;

  begin
    new_line;
    put_line("Testing the creation of standard Jacobian matrices.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    n := lp'last;
    declare
      jm : constant Standard_Complex_Jaco_Matrices.Jaco_Mat(1..n,1..n)
         := Create(lp.all);
    begin
      put_line("The Jacobian matrix : ");
      for i in jm'range(1) loop
        for j in jm'range(2) loop
          put(jm(i,j));
        end loop;
        new_line;
      end loop;
      dj := Expand(jm);
      put_line("The determinant of the Jacobian matrix : "); put(dj);
      new_line;
      Test_Expand_by_Eval(n,jm,dj);
    end;
  end Test_Standard_Expansion;

  procedure Test_Multprec_Creation is

  -- DESCRIPTION :
  --   Tests the creation of a multi-precision Jacobian matrix.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    n : integer32;

  begin
    new_line;
    put_line("Testing the creation of multi-precision Jacobian matrices.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    n := lp'last;
    declare
      mp : constant Multprec_Complex_Poly_Systems.Poly_Sys(1..n)
         := Convert(lp.all);
      jm : constant Multprec_Complex_Jaco_Matrices.Jaco_Mat(1..n,1..n)
         := Create(mp);
    begin
      put_line("The Jacobian matrix : ");
      for i in jm'range(1) loop
        for j in jm'range(2) loop
          put(jm(i,j));
        end loop;
        new_line;
      end loop;
    end;
  end Test_Multprec_Creation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of the operations on complex polynomials.");
    loop
      new_line;
      put_line("Choose one of the following :                              ");
      put_line("  0. Exit this program.                                    ");
      put_line("  1. Test creation of standard Jacobian matrices.          ");
      put_line("  2. Test expansion of standard Jacobian matrices.         ");
      put_line("  3. Test creation of multi-precision Jacobian matrices.   ");
      put("Type 0,1,2 or 3 to select : "); get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Test_Standard_Creation;
        when '2' => Test_Standard_Expansion;
        when '3' => Test_Multprec_Creation;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_jaco;
