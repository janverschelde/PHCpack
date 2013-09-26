with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Affine_Binomial_Iterator;

package body Standard_Affine_Binomials is

  procedure Extract_Two_Terms
              ( p : in Standard_Complex_Laurentials.Poly;
                t1,t2 : out Standard_Complex_Laurentials.Term;
                fail : out boolean ) is

    use Standard_Complex_Laurentials;
    first_time : boolean := true;

    procedure Extract_Term ( t : in Term; continue : out boolean ) is
    begin
      if first_time then
        t1 := t;
        first_time := false;
        continue := true;
      else
        t2 := t;
        continue := false;
      end if;
    end Extract_Term;
    procedure Extract_Terms is new Visiting_Iterator(Extract_Term);

  begin
    fail := (Number_of_Terms(p) /= 2);
    if not fail
     then Extract_Terms(p);
    end if;
  end Extract_Two_Terms;

  procedure Incidence_Matrix
               ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 A : out Standard_Integer_Matrices.Matrix;
                 fail : out boolean ) is

    use Standard_Complex_Laurentials;
    n : constant integer32 := A'last(2);
    t1,t2 : Term;
    ind : integer32 := 0;

  begin
    for i in p'range loop
      Extract_Two_Terms(p(i),t1,t2,fail);
      exit when fail;
      for k in 1..n loop
        if t1.dg(k) > 0 
         then A(ind+1,k) := 1;
         else A(ind+1,k) := 0;
        end if;
        if t2.dg(k) > 0 
         then A(ind+2,k) := 1;
         else A(ind+2,k) := 0;
        end if;
      end loop;
      ind := ind + 2;
     -- Clear(t1); Clear(t2); --> sharing !!!
    end loop;
  end Incidence_Matrix;

  procedure Nonzero_Binomials
             ( A : in Standard_Integer_Matrices.Matrix;
               s : in Standard_Integer_Vectors.Vector;
               e : out Standard_Integer_Vectors.Vector;
               cnt : out integer32; valid : out boolean ) is
 
    nq : constant integer32 := A'last(1)/2;
    m : integer32;

  begin
    cnt := 0; valid := true;
    for i in 1..nq loop
      m := 2*(i-1) + 1;  -- index of first monomial in i-th equation
      if Affine_Binomial_Iterator.Set_to_Zero(A,m,s) then
        if Affine_Binomial_Iterator.Set_to_Zero(A,m+1,s)
         then e(i) := 0;
         else e(i) := -1; cnt := cnt + 1; valid := false;
        end if;
      else
        if Affine_Binomial_Iterator.Set_to_Zero(A,m+1,s)
         then e(i) := -1; valid := false;
         else e(i) := 1; cnt := cnt + 1;
        end if;
      end if;
    end loop;
  end Nonzero_Binomials;

  function Free_Variables
             ( A : Standard_Integer_Matrices.Matrix;
               s : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(s'range);

  begin
    for i in res'range loop
      if s(i) = 1 then
        res(i) := 0;        -- not free if set to zero
      else
        res(i) := 1;                 -- assume is free
        for k in A'range(1) loop     -- check monomial k
          if A(k,i) /= 0 then        -- occurs in monomial
            if not Affine_Binomial_Iterator.Set_to_Zero(A,k,s)
             then res(i) := 0;       -- not free if monomial remains
            end if;
          end if;
          exit when (res(i) = 0);
        end loop;
      end if;
    end loop;
    return res;
  end Free_Variables;

  function Subsystem
             ( p : Standard_Complex_Laur_Systems.Laur_Sys;
               nq : integer32;
               eq : Standard_Integer_Vectors.Vector ) 
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(1..nq);
    ind : integer32 := 0;

  begin
    for i in eq'range loop
      if eq(i) = 1
       then ind := ind+1; res(ind) := p(i);
      end if;
    end loop;
    return res;
  end Subsystem;

  function Eliminate_Variables
             ( x,s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(x'first..x'last-s_cnt);
    ind : integer32 := x'first-1;

  begin
    for i in x'range loop
      if s(i) = 0 then
        ind := ind + 1;
        res(ind) := x(i);
      end if;
    end loop;
    return res;
  end Eliminate_Variables;

  function Eliminate_Variables
             ( p : Standard_Complex_Laurentials.Poly;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Standard_Complex_Laurentials.Poly := Null_Poly;

    procedure Eliminate ( t : in Term; continue : out boolean ) is

      nt : Term;
      ind : integer32 := 0;

    begin
      nt.cf := t.cf;
      nt.dg := new Standard_Integer_Vectors.Vector(1..t.dg'last-s_cnt);
      for i in t.dg'range loop
        if s(i) = 0 then
          ind := ind + 1;
          nt.dg(ind) := t.dg(i);
        end if;
      end loop;
      Add(res,nt);
      continue := true;
    end Eliminate;
    procedure Eliminate_Terms is new Visiting_Iterator(Eliminate);

  begin
    Eliminate_Terms(p);
    return res;
  end Eliminate_Variables;

  function Eliminate_Variables
             ( p : Standard_Complex_Laur_Systems.Laur_Sys;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Eliminate_Variables(p(i),s,s_cnt);
    end loop;
    return res;
  end Eliminate_Variables;

  function Insert_Zero_Values
             ( c : Standard_Complex_Vectors.Vector;
               s : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(s'range);
    ind : integer32 := c'first-1;

  begin
    for i in s'range loop
      if s(i) = 1 then -- selected i-th coordinate to be zero
        res(i) := Standard_Complex_Numbers.Create(0.0);
      else
        ind := ind + 1;         -- take next component from c
        res(i) := c(ind);
      end if;
    end loop;
    return res;
  end Insert_Zero_Values;

  function Insert_Zero_Values
             ( c : Standard_Complex_Solutions.Solution;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(c.n+s_cnt);

  begin
    put("  c.n = "); put(c.n,1); 
    put("  s_cnt = "); put(s_cnt,1); 
    put("  res.n = "); put(res.n,1); new_line; 
    res.t := c.t;
    res.m := c.m;
    res.err := c.err;
    res.rco := c.err;
    res.res := c.res;
    res.v := Insert_Zero_Values(c.v,s);
    return res;
  end Insert_Zero_Values;

  function Insert_Zero_Values
             ( c : Standard_Complex_Solutions.Solution_List;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    tmp : Solution_List := c;
    res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Insert_Zero_Values(ls.all,s,s_cnt));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert_Zero_Values;

  function Insert_Unit_Vectors
             ( M : Standard_Integer_Matrices.Matrix;
               s : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(s'range,s'range);
    row,col : integer32;

  begin
    row := M'first(1) - 1;
    for i in s'range loop
      if s(i) = 1 then
        for j in res'range(2) loop
          res(i,j) := 0;
        end loop;
        res(i,i) := 1;
      else
        row := row + 1;
        col := M'first(2) - 1;
        for j in s'range loop
          if s(j) = 0 then
            col := col + 1;
            res(i,j) := M(row,col);
          else
            if j = i
             then res(i,j) := 1;
             else res(i,j) := 0;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    return res;
  end Insert_Unit_Vectors;

end Standard_Affine_Binomials;
