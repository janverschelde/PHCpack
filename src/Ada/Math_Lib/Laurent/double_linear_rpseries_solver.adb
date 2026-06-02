with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Standard_Floating_Matrices_IO;     use Standard_Floating_Matrices_IO;
with Standard_Complex_Singular_Values;
with Double_Puiseux_Operations;
with Double_Real_Power_Series;
with Double_Real_Power_Series_IO;

package body Double_Linear_rpSeries_Solver is

  function Right_Hand_Side
             ( A : Double_rpSeries_Matrices.Matrix;
               x : in Double_rpSeries_Vectors.Vector ) 
             return Double_rpSeries_Vectors.Vector is

    res : Double_rpSeries_Vectors.Vector(A'range);

    use Double_Real_Power_Series;

  begin
    for i in A'range(1) loop
      res(i) := A(i,A'first(2))*x(x'first);
      for j in A'first(2)+1..A'last(2) loop
        declare
          prd : Link_to_Series := A(i,j)*x(j);
          sum : Link_to_Series := res(i) + prd;
        begin
          Copy(sum,res(i)); Clear(prd); clear(sum);
        end;
      end loop;
    end loop;
    return res;
  end Right_Hand_Side;

  function Extract_Constants
             ( A : Double_rpSeries_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).cff(0);
      end loop;
    end loop;
    return res;
  end Extract_Constants;

  function Extract_Leading_Powers
             ( A : Double_rpSeries_Matrices.Matrix )
             return Standard_Floating_Matrices.Matrix is

    res : Standard_Floating_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).pwt(1);
      end loop;
    end loop;
    return res;
  end Extract_Leading_Powers;

  function Extract_Constants 
             ( v : Double_rpSeries_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := v(i).cff(0);
    end loop;
    return res;
  end Extract_Constants;

  function Inverse ( A : Standard_Complex_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix is

    wrk : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2)) := A;
    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Standard_Complex_Vectors.Vector(1..mm);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Standard_Complex_Matrices.Matrix(1..n,1..n);
    v : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Standard_Complex_Singular_Values.SVD(wrk,n,p,s,e,u,v,job,info);
    res := Standard_Complex_Singular_Values.Inverse(u,v,s);
    return res;
  end Inverse;

  function Matrix_Multiply
             ( A : Standard_Complex_Matrices.Matrix;
               x : Double_rpSeries_Vectors.Vector ) 
             return Double_rpSeries_Vectors.Vector is

    res : Double_rpSeries_Vectors.Vector(A'range);

    use Double_Real_Power_Series;

  begin
    for i in A'range(1) loop
      res(i) := A(i,A'first(2))*x(x'first);
      for j in A'first(2)+1..A'last(2) loop
        declare
          prd : Link_to_Series := A(i,j)*x(j);
          sum : Link_to_Series := res(i) + prd;
        begin
          Copy(sum,res(i)); Clear(prd); clear(sum);
        end;
      end loop;
    end loop;
    return res;
  end Matrix_Multiply;

  function Is_In ( A : Standard_Floating_Matrices.Matrix;
                   nbr : double_float; tol : double_float := 1.0E-12 )
                 return boolean is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if abs(A(i,j) - nbr) < tol
         then return true;
        end if;
      end loop;
    end loop;
    return false;
  end Is_In;

  procedure Leading_Term
              ( invAb : in Double_rpSeries_Vectors.Vector;
                rA : in Standard_Floating_Matrices.Matrix;
                leadidx : out integer32; 
                leadpow : out double_float; leadcff : out complex_number;
                tol : in double_float := 1.0E-12 ) is

    sb : Double_Real_Power_Series.Link_to_Series;
    foundfirst : boolean;

  begin
    for i in invAb'range loop
      sb := invAb(i);
      foundfirst := false;
      for j in sb.pwt'range loop
        if not Is_In(rA,sb.pwt(j),tol) then
          foundfirst := true;
          if i = invAb'first then
            leadidx := 1;
            leadpow := sb.pwt(j);
            leadcff := sb.cff(j);
          elsif sb.pwt(j) < leadpow then
            leadidx := i;
            leadpow := sb.pwt(j);
            leadcff := sb.cff(j);
          end if;
        end if;
        exit when foundfirst;
      end loop;
    end loop;
  end Leading_Term;

  procedure Next_Term
              ( invAb : in Double_rpSeries_Vectors.Vector;
                rA : in Standard_Floating_Matrices.Matrix;
                powers : in Standard_Floating_Vectors.Vector;
                leadidx : out integer32; 
                leadpow : out double_float; leadcff : out complex_number;
                tol : in double_float := 1.0E-12 ) is

    sb : Double_Real_Power_Series.Link_to_Series;
    found,foundfirst : boolean;

  begin
    leadidx := -1;
    for i in invAb'range loop
      if powers(i) < 0.0 then -- not yet computed i-th index
        sb := invAb(i);
        foundfirst := false;
        for j in sb.pwt'range loop
          if not Is_In(rA,sb.pwt(j),tol) then
            found := false;
            for k in powers'range loop
              if powers(k) > 0.0
               then found := Is_In(rA,sb.pwt(j)-powers(k),tol);
              end if;
              exit when found;
            end loop;
            if not found then
              foundfirst := true;
              if leadidx = -1 then
                leadidx := i;
                leadpow := sb.pwt(j);
                leadcff := sb.cff(j);
              elsif sb.pwt(j) < leadpow then
                leadidx := i;
                leadpow := sb.pwt(j);
                leadcff := sb.cff(j);
              end if;
            end if;
          end if;
          exit when foundfirst;
        end loop;
      end if;
    end loop;
  end Next_Term;

  procedure Real_Power_Series_Solver 
              ( A : in Double_rpSeries_Matrices.Matrix;
                b : in Double_rpSeries_Vectors.Vector;
                z0,z1c : out Standard_Complex_Vectors.Vector;
                z1p : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    A0 : constant Standard_Complex_Matrices.Matrix(A'range(1),A'range(2))
       := Extract_Constants(A);
    rA : constant Standard_Floating_Matrices.Matrix(A'range(1),A'range(2))
       := Extract_Leading_Powers(A);
    b0 : constant Standard_Complex_Vectors.Vector(b'range)
       := Extract_Constants(b);
    rcond : double_float;
    invA0 : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));
    invAb : Double_rpSeries_Vectors.Vector(b'range);
    leadidx : integer32;
    leadcff : complex_number;
    leadpow : double_float;

  begin
    Double_Puiseux_Operations.Solve_Constant_Linear_System
      (A0,b0,z0,rcond,vrblvl-1);
    if vrblvl > 0 then
      put_line("the computed constant coefficients :");
      put_line(z0);
      put("rcond :"); put(rcond,3); new_line;
    end if;
    invA0 := Inverse(A0);
    invAb := Matrix_Multiply(invA0,b);
    if vrblvl > 0 then
      put_line("multiplied right hand side vector :"); 
      Double_Real_Power_Series_IO.Write(invAb);
      put_line("leading powers of the coefficient matrix :"); put(rA);
    end if;
    Leading_Term(invAb,rA,leadidx,leadpow,leadcff);
    z1p := (z1p'range => -1.0);
    z1p(leadidx) := leadpow;
    z1c(leadidx) := leadcff;
    if vrblvl > 0 then
      put("leading index : "); put(leadidx,1); new_line;
      put(leadcff); put(" * t**"); put(leadpow,1,14,3); new_line;
    end if;
    for k in 2..z0'last loop
      if vrblvl > 0
       then put("computing term "); put(k,1); put_line(" ...");
      end if;
      Next_Term(invAb,rA,z1p,leadidx,leadpow,leadcff);
      if vrblvl > 0 then
        put("next index : "); put(leadidx,1); new_line;
        put(leadcff); put(" * t**"); put(leadpow,1,14,3); new_line;
        z1p(leadidx) := leadpow;
        z1c(leadidx) := leadcff;
      end if;
    end loop;
  end Real_Power_Series_Solver;

end Double_Linear_rpSeries_Solver;
