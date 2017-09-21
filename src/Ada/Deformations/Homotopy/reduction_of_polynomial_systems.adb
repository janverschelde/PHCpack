with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with Reduction_of_Polynomials;           use Reduction_of_Polynomials;
with Standard_Linear_Reduction;
with DoblDobl_Linear_Reduction;
with QuadDobl_Linear_Reduction;

package body Reduction_of_Polynomial_Systems is

  procedure Reduce ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean ) is

    use Standard_Complex_Matrices;
    use Standard_Complex_Polynomials;
    use Standard_Linear_Reduction;

    n : constant integer32 := p'length;
    max_terms : constant integer32 := integer32(Sum_Number_of_Terms(p));
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural32 := 0;
    mat : Matrix(p'range,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns,diagonal);
    if diagonal then
      inconsistent := false;
      infinite := false;
    else
      declare
        coeffmat : Matrix(p'range,1..integer32(numb_columns));
        tol : constant double_float := 10.0**(-8);
      begin
        for i in coeffmat'range(1) loop
          for j in coeffmat'range(2) loop
            coeffmat(i,j) := mat(i,j);
          end loop;
        end loop;
        Triangulate(coeffmat,tol,n,integer32(numb_columns));
        Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                               inconsistent,infinite);
        for i in 1..integer32(numb_Columns) loop
          Standard_Natural_Vectors.Clear
            (Standard_Natural_Vectors.Link_to_Vector(columns(i)));
        end loop;
      end;
    end if;
  end Reduce;

  procedure Reduce ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean ) is

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Linear_Reduction;

    n : constant integer32 := p'length;
    max_terms : constant integer32 := integer32(Sum_Number_of_Terms(p));
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural32 := 0;
    mat : Matrix(p'range,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns,diagonal);
    if diagonal then
      inconsistent := false;
      infinite := false;
    else
      declare
        coeffmat : Matrix(p'range,1..integer32(numb_columns));
        tol : constant double_float := 10.0**(-8);
      begin
        for i in coeffmat'range(1) loop
          for j in coeffmat'range(2) loop
            coeffmat(i,j) := mat(i,j);
          end loop;
        end loop;
        Triangulate(coeffmat,tol,n,integer32(numb_columns));
        Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                               inconsistent,infinite);
        for i in 1..integer32(numb_Columns) loop
          Standard_Natural_Vectors.Clear
            (Standard_Natural_Vectors.Link_to_Vector(columns(i)));
        end loop;
      end;
    end if;
  end Reduce;

  procedure Reduce ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean ) is

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Linear_Reduction;

    n : constant integer32 := p'length;
    max_terms : constant integer32 := integer32(Sum_Number_of_Terms(p));
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural32 := 0;
    mat : Matrix(p'range,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns,diagonal);
    if diagonal then
      inconsistent := false;
      infinite := false;
    else
      declare
        coeffmat : Matrix(p'range,1..integer32(numb_columns));
        tol : constant double_float := 10.0**(-8);
      begin
        for i in coeffmat'range(1) loop
          for j in coeffmat'range(2) loop
            coeffmat(i,j) := mat(i,j);
          end loop;
        end loop;
        Triangulate(coeffmat,tol,n,integer32(numb_columns));
        Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                               inconsistent,infinite);
        for i in 1..integer32(numb_Columns) loop
          Standard_Natural_Vectors.Clear
            (Standard_Natural_Vectors.Link_to_Vector(columns(i)));
        end loop;
      end;
    end if;
  end Reduce;

  procedure Sparse_Reduce ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                            inconsistent,infinite : in out boolean ) is

    use Standard_Complex_Matrices;
    use Standard_Complex_Polynomials;
    use Standard_Linear_Reduction;

    n : constant integer32 := p'length;
    max_terms : constant integer32 := integer32(Sum_Number_of_Terms(p));
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural32 := 0;
    mat : Matrix(1..n,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns);
    declare
      coeffmat : Matrix(p'range,1..integer32(numb_columns));
    begin
      for i in coeffmat'range(1) loop
        for j in coeffmat'range(2) loop
          coeffmat(i,j) := mat(i,j);
        end loop;
      end loop;
      Diagonalize(coeffmat,n,integer32(numb_Columns));
      Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                             inconsistent,infinite);
      for i in 1..integer32(numb_columns) loop
        Standard_Natural_Vectors.Clear
          (Standard_Natural_Vectors.Link_to_Vector(columns(i)));
      end loop;
    end;
  end Sparse_Reduce;

  procedure Sparse_Reduce ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                            inconsistent,infinite : in out boolean ) is

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Linear_Reduction;

    n : constant integer32 := p'length;
    max_terms : constant integer32 := integer32(Sum_Number_of_Terms(p));
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural32 := 0;
    mat : Matrix(1..n,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns);
    declare
      coeffmat : Matrix(p'range,1..integer32(numb_columns));
    begin
      for i in coeffmat'range(1) loop
        for j in coeffmat'range(2) loop
          coeffmat(i,j) := mat(i,j);
        end loop;
      end loop;
      Diagonalize(coeffmat,n,integer32(numb_Columns));
      Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                             inconsistent,infinite);
      for i in 1..integer32(numb_columns) loop
        Standard_Natural_Vectors.Clear
          (Standard_Natural_Vectors.Link_to_Vector(columns(i)));
      end loop;
    end;
  end Sparse_Reduce;

  procedure Sparse_Reduce ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                            inconsistent,infinite : in out boolean ) is

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Linear_Reduction;

    n : constant integer32 := p'length;
    max_terms : constant integer32 := integer32(Sum_Number_of_Terms(p));
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural32 := 0;
    mat : Matrix(1..n,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns);
    declare
      coeffmat : Matrix(p'range,1..integer32(numb_columns));
    begin
      for i in coeffmat'range(1) loop
        for j in coeffmat'range(2) loop
          coeffmat(i,j) := mat(i,j);
        end loop;
      end loop;
      Diagonalize(coeffmat,n,integer32(numb_Columns));
      Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                             inconsistent,infinite);
      for i in 1..integer32(numb_columns) loop
        Standard_Natural_Vectors.Clear
          (Standard_Natural_Vectors.Link_to_Vector(columns(i)));
      end loop;
    end;
  end Sparse_Reduce;

-- NONLINEAR REDUCTION :

  function LEQ ( d1,d2 : Standard_Complex_Polynomials.Degrees )
               return boolean is

  -- DESCRIPTION :
  --   Returns true if all degrees of d1 are lower than
  --   or equal to the degrees of d2

  begin
    for i in d1'range loop
      if d1(i) > d2(i)
       then return false;
      end if;
    end loop;
    return true;
  end LEQ;

  function Leading_Term
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Returns the leading term of the polynomial p.

    use Standard_Complex_Polynomials;

    tf : Term;

    procedure First_Term (t : in Term; continue : out boolean) is
    begin
      Copy(t,tf);
      continue := false;
    end First_Term;
    procedure Get_First_Term is new Visiting_Iterator (First_Term);
  begin
    Get_First_Term(p);
    return tf;
  end Leading_Term;

  function Can_Be_Eliminated
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; j : integer32 )
             return boolean is

  -- DESCRIPTION :
  --   returns true if the degree of the j-th unknown in each equation 
  --   is zero.

    use Standard_Complex_Polynomials;

  begin
    for i in p'range loop
      if Degree(p(i),j) > 0
       then return false;
      end if;
    end loop;
    return true;
  end Can_Be_Eliminated;

  procedure Shift_Null_Polynomial
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   The null polynomial in the system p will be shifted down
  --   towards the end.

    use Standard_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) = Null_Poly then
        for j in i..(p'last-1) loop
          Copy(p(j+1),p(j));
          Clear(p(j+1));
        end loop;
      end if;
    end loop;
  end Shift_Null_Polynomial;

  procedure Eliminate
              ( p : in out Standard_Complex_Polynomials.Poly;
                j : in integer32 ) is

  -- DESCRIPTION :
  --   The j-th unknown will be eliminated out of the polynomial p

    use Standard_Complex_Polynomials;

    n : constant integer32 := integer32(Number_Of_Unknowns(p));

    procedure Eliminate_Term (t : in out Term; continue : out boolean) is

      d : constant Standard_Complex_Polynomials.Degrees
        := new Standard_Natural_Vectors.Vector(1..(n-1));

    begin
      for i in 1..(j-1) loop
        d(i) := t.dg(i);
      end loop;
      for i in j..(n-1) loop
        d(i) := t.dg(i+1);
      end loop;
      Clear(t);
      t.dg := d;
      continue := true;
    end Eliminate_Term;
    procedure Eliminate_Terms is new Changing_Iterator(Eliminate_Term);

  begin
    Eliminate_Terms(p);
  end Eliminate;
        
  procedure Eliminate
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                j : in integer32 ) is

  -- DESCRIPTION :
  --   The j-th unknown will be eliminated out of each equation.

  begin
    for i in p'range loop
      Eliminate(p(i),j);
    end loop;
  end Eliminate;
    
  procedure Replace
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                pp : in Standard_Complex_Polynomials.Poly;
                i : in integer32 ) is

  -- DESCRIPTION :
  --   This procedure replaces the i-th polynomial in the system p
  --   by the polynomial pp.  If pp is a null polynomial then the procedure
  --   tries to eliminate an unknown, in order to have as much equations
  --   as there are unknowns.

    use Standard_Complex_Polynomials;

    tmp : natural32;

  begin
    if (pp = Null_Poly) or else (Number_Of_Unknowns(pp) = 0) then
      -- try to eliminate an unknown
      tmp := Number_Of_Unknowns(p(1));
      Clear(p(i)); p(i) := Null_Poly;
      for j in reverse 1..integer32(Number_Of_Unknowns(p(1))) loop
        if Can_Be_Eliminated(p,j)
         then Eliminate(p,j);
        end if;
      end loop;
      Shift_Null_Polynomial(p);
    else
      Clear(p(i)); Copy(pp,p(i));
    end if;
  end Replace;

  function red ( p,b1,b2 : Standard_Complex_Polynomials.Poly )
               return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    Rpb1 : Poly := Rpoly(p,b1);

  begin
    if Number_Of_Unknowns(Rpb1) = 0 then
      return Null_Poly;
    else
      declare
        Rpb2 : Poly := Rpoly(Rpb1,b2);
      begin
        Clear(Rpb1);
        return Rpb2;
      end;
    end if;
  end red;

  function Reduce ( p,b1,b2 : Standard_Complex_Polynomials.Poly )
                  return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   returns p mod < b1,b2 >

    use Standard_Complex_Polynomials;

    temp : Poly := red(p,b1,b2);

  begin
    if Number_Of_Unknowns(temp) = 0 then
      return Null_Poly;
    else
      Clear(temp);
      return red(p,b2,b1); 
    end if;
  end Reduce;

  function Simple_Criterium
             ( p1,p2 : Standard_Complex_Polynomials.Poly ) return boolean is

  -- DESCRIPTION :
  --   returns true if lcm(in(p1),in(p2)) = in(p1), if in(p2) | in(p1).

    use Standard_Complex_Polynomials;

    lt1,lt2 : Term;
    res : boolean;

  begin
    lt1 := Leading_Term(p1);
    lt2 := Leading_Term(p2);
    res := LEQ(lt2.dg,lt1.dg);
    Clear(lt1); Clear(lt2);
    return res;
  end Simple_Criterium;

  procedure Rpoly_Criterium
              ( p,b1,b2 : in Standard_Complex_Polynomials.Poly;
                cnt : in out natural32; res : out boolean ) is

  -- DESCRIPTION :
  --   Applies the R-polynomial criterium and counts the number of 
  --   R-polynomials computed.

    use Standard_Complex_Polynomials;

    Rpb1 : Poly := Rpoly(p,b1);
    Rpb2 : Poly;

  begin
    cnt := cnt + 1;
    if Number_Of_Unknowns(Rpb1) = 0 then
      res := true;
    else
      Rpb2 := Rpoly(Rpb1,b2);
      cnt := cnt + 1;
      Clear(Rpb1);
      if Number_of_Unknowns(Rpb2) = 0 then
        res := true;
      else
        Clear(Rpb2);
        Rpb2 := Rpoly(p,b2);
        cnt := cnt + 1;
        if Number_of_Unknowns(Rpb2) = 0 then
          res := true;
        else
          Rpb1 := Rpoly(Rpb2,b1);
          cnt := cnt + 1;
          Clear(Rpb2);
          if Number_of_Unknowns(Rpb1) = 0
           then res := true;
           else res := false;
                Clear(Rpb1);
          end if;
        end if;
      end if;
    end if;
  end Rpoly_Criterium;

  function Criterium
             ( p,q,s : Standard_Complex_Polynomials.Poly ) return boolean is

  -- DESCRIPTION :
  --   returns true if p may be replaced by s.

    use Standard_Complex_Polynomials;

  begin
    if Simple_Criterium(p,q) then
      return true;
    else
      declare
        temp : Poly := Reduce(p,q,s);
        res : constant boolean := (Number_Of_Unknowns(temp) = 0);
      begin
        Clear(temp);
        return res;
      end;
    end if;
  end Criterium;

  procedure Criterium
              ( p,q,s : in Standard_Complex_Polynomials.Poly;
                cnt : in out natural32; res : out boolean ) is

  -- DESCRIPTION :
  --   returns true if p may be replaced by s.

  begin
    if Simple_Criterium(p,q)
     then res := true;
     else Rpoly_Criterium(p,q,s,cnt,res);
    end if;
  end Criterium;

  procedure Reduce ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     res : in out Standard_Complex_Poly_Systems.Poly_Sys;
                     cnt_eq : in out natural32; max_eq : in natural32;
                     cnt_sp : in out natural32; max_sp : in natural32;
                     cnt_rp : in out natural32; max_rp : in natural32 ) is

    use Standard_Complex_Polynomials;
 
    S : Poly;
    n : constant integer32 := p'last - p'first + 1;
    dS,dpi,dpj : integer32;
    ok : boolean;

    procedure try ( i,dpi : in integer32 ) is

    -- DESCRIPTION : try to replace p_i by S

      p_red : Standard_Complex_Poly_Systems.Poly_Sys(1..n);

    begin
      if cnt_eq > max_eq then return; end if;
      if cnt_sp > max_sp then return; end if;
      Standard_Complex_Poly_Systems.Clear(p_red);
      Standard_Complex_Poly_Systems.Copy(p,p_red);
      Replace(p_red,S,i);
     -- put("replaced polynomial p("); put(i,1); put_line(").");
      if dS = 0 then
        return;
      elsif Total_Degree(p_red) < Total_Degree(res) then
        Standard_Complex_Poly_Systems.Copy(p_red,res); 
        Reduce(p_red,res,cnt_eq,max_eq,cnt_sp,max_sp,cnt_rp,max_rp);
      elsif cnt_eq <= max_eq then
        cnt_eq := cnt_eq + 1;
        Reduce(p_red,res,cnt_eq,max_eq,cnt_sp,max_sp,cnt_rp,max_rp);
      end if;
      Standard_Complex_Poly_Systems.Clear(p_red);
    end try;

  begin
    if cnt_eq > max_eq then return; end if;
    if cnt_sp > max_sp then return; end if;
    if cnt_rp > max_rp then return; end if;
    for i in 1..n loop
      for j in (i+1)..n loop
        if (p(i) /= Null_Poly) and (p(j) /= Null_Poly) then
          Clear(S); S := Spoly(p(i),p(j));
          cnt_sp := cnt_sp + 1;
          dS  := Degree(S); dpi := Degree(p(i)); dpj := Degree(p(j));
          -- put("S-polynomial of p("); put(i,1); put(") and p(");
          --   put(j,1); put_line(") computed.");
          if dS <= dpi and then dpi > dpj
             and then Criterium(p(i),p(j),S) then
            try(i,dpi);
          elsif dS <= dpj and then dpi < dpj
                 and then Criterium(p(j),p(i),S) then
            try(j,dpj);
          else -- dpi = dpj
            if dS <= dpi then
              Criterium(p(i),p(j),S,cnt_rp,ok);
              if ok then try(i,dpi); end if;
            end if;
            if dS <= dpj
             then Criterium(p(j),p(i),S,cnt_rp,ok);
                  if ok then try(j,dpj); end if;
            end if;
          end if;
          Clear(S);
        end if;
        exit when (dS = 0);
      end loop;
    end loop;
  end Reduce;

  procedure Sparse_Reduce
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                res : in out Standard_Complex_Poly_Systems.Poly_Sys;
                cnt_eq : in out natural32; max_eq : in natural32 ) is

    use Standard_Complex_Polynomials;

    S : Poly;
    n : constant integer32 := p'last - p'first + 1;
    dS,dpi,dpj : integer32;

    procedure try ( i,dpi : in integer32 ) is

    -- DESCRIPTION : try to replace p_i by S

      p_red : Standard_Complex_Poly_Systems.Poly_Sys(1..n);
      inconsistent,infinite : boolean := false;

    begin
      if cnt_eq > max_eq then return; end if;
      Standard_Complex_Poly_Systems.Clear(p_red);
      Standard_Complex_Poly_Systems.Copy(p,p_red);
      Replace(p_red,S,i);
      if dS /= 0
       then Sparse_Reduce(p_red,inconsistent,infinite);
      end if;
      if dS = 0 or inconsistent then
        cnt_eq := max_eq + 1; return;
      elsif Total_Degree(p_red) < Total_Degree(res) then
        Standard_Complex_Poly_Systems.Copy(p_red,res); 
        Sparse_Reduce(p_red,res,cnt_eq,max_eq);
      else
        cnt_eq := cnt_eq + 1;
        Sparse_Reduce(p_red,res,cnt_eq,max_eq);
      end if;
      Standard_Complex_Poly_Systems.Clear(p_red);
    end try;

  begin
    if cnt_eq > max_eq then return; end if;
    for i in 1..n loop
      for j in (i+1)..n loop
        if (p(i) /= Null_Poly) and (p(j) /= Null_Poly) then
          Clear(S); S := Spoly(p(i),p(j));
          dS  := Degree(S);
          dpi := Degree(p(i));
          dpj := Degree(p(j));
          if dS <= dpi and then dpi > dpj 
                       and then Criterium(p(i),p(j),S) then
            try(i,dpi);
          elsif dS <= dpj and then dpi < dpj 
                          and then Criterium(p(j),p(i),S) then
            try(j,dpj);
          else -- dpi = dpj
            if dS <= dpi
             and then Criterium(p(i),p(j),S) then try(i,dpi); end if;
            if dS <= dpj
              and then Criterium(p(j),p(i),S) then try(j,dpj); end if;
          end if;
          Clear(S);
        end if;
      end loop;
    end loop;
  end Sparse_Reduce;

end Reduction_of_Polynomial_Systems;
