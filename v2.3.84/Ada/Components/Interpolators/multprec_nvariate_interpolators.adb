with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;

package body Multprec_Nvariate_Interpolators is

-- AUXILIARIES TO CONSTRUCTORS :

  procedure Adjust_Divided_Differences 
               ( divisor : in Complex_Number; v : in NesVec;
                 w : in out NesVec ) is 

  -- DESCRIPTION :
  --   Computes divided differences across the nested vectors.

  begin
    if v.n = 1 then
      for i in v.v'range loop
        Min(w.v(i));
        Add(w.v(i),v.v(i));
        Div(w.v(i),divisor);
      end loop;
    else
      for i in v.w'range loop
        Adjust_Divided_Differences(divisor,v.w(i).all,w.w(i).all);
      end loop;
    end if;
  end Adjust_Divided_Differences;

-- CONSTRUCTORS :

  function Create ( n,d : natural32; x : VecVec; y : NesVec ) return NesVec is

    res : NesVec(n,0,integer32(d));

    procedure Rec_Create ( k,b : in integer32; v : in out NesVec ) is

    -- DESCRIPTION :
    --   Going through the levels of the nested vector v, the generalized
    --   divided differences are computed and stored in v.
    --   No extra higher order divided differences are computed.

    -- ON ENTRY :
    --   k         current level in the nested vector;
    --   b         last entry in the triangular grid;
    --   v         nested vector, originally with function values.

    -- ON RETURN :
    --   v         updated with divided differences.

      divisor : Complex_Number;

    begin
      if v.n = 1 then
        for i in 1..b loop         -- for square grid : b <-> v.b
          for j in 0..(i-1) loop
            divisor := x(k)(j) - x(k)(i);
            Min(v.v(i));
            Add(v.v(i),v.v(j));
            Div(v.v(i),divisor);
            Clear(divisor);
          end loop;
        end loop;
      else 
        for i in v.a..b loop       -- for square grid : b <-> v.b
          Rec_Create(k+1,b-i,v.w(i).all);
        end loop;
        for ii in 1..b loop
          for jj in 0..(ii-1) loop
            divisor := x(k)(jj) - x(k)(ii);
            Adjust_Divided_Differences(divisor,v.w(jj).all,v.w(ii).all);
            Clear(divisor);
          end loop;
        end loop;
      end if;
    end Rec_Create;

  begin
    Copy(y,res);
    Rec_Create(1,y.b,res);
    return res;
  end Create;

  function Create_on_Square
             ( n,d : natural32; x : VecVec; y : NesVec ) return NesVec is

    res : NesVec(n,0,integer32(d));

    procedure Rec_Create ( k,b : in integer32; v : in out NesVec ) is

    -- DESCRIPTION :
    --   Going through the levels of the nested vector v, the generalized
    --   divided differences are computed and stored in v.
    --   For testing we also compute higher order differences.

    -- ON ENTRY :
    --   k         current level in the nested vector;
    --   b         last entry in the triangular grid;
    --   v         nested vector, originally with function values.

    -- ON RETURN :
    --   v         updated with divided differences.

      divisor : Complex_Number;

    begin
      if v.n = 1 then
        for i in 1..v.b loop       -- for triangular grid : v.b <-> b
          for j in 0..(i-1) loop
            divisor := x(k)(j) - x(k)(i);
            Min(v.v(i));
            Add(v.v(i),v.v(j));
            Div(v.v(i),divisor);
            Clear(divisor);
          end loop;
        end loop;
      else 
        for i in v.a..v.b loop     -- for triangular grid : v.b <-> b
          Rec_Create(k+1,v.b-i,v.w(i).all);
        end loop;
       for ii in 1..v.b loop
          for jj in 0..(ii-1) loop
            divisor := x(k)(jj) - x(k)(ii);
            Adjust_Divided_Differences(divisor,v.w(jj).all,v.w(ii).all);
            Clear(divisor);
          end loop;
        end loop;
      end if;
    end Rec_Create;

  begin
    Copy(y,res);
    Rec_Create(1,y.b,res);
    return res;
  end Create_on_Square;

  function Expand ( f : NesVec; x : VecVec ) return Poly is
 
    res : Poly := Null_Poly;
    n : constant integer32 := x'length;
    acc : Standard_Integer_Vectors.Vector(x'range) := (x'range => 0);
    fac : Term;
 
    procedure Rec_Expand (  v : in Nesvec; k,b : in integer32 ) is
 
    -- DESCRIPTION : recursive elaboration of the evaluation order.
 
    -- ON ENTRY :
    --   v         nested vector with generalized divided differences;
    --   k         current entry in the accumulating index acc;
    --   b         last entry in the triangular grid.
 
      trm,sum : Poly;
 
    begin
      if v.n = 1 then
        for i in reverse v.a..b loop
          acc(k) := i;
          Copy(v.v(i),fac.cf);
          trm := Create(fac);
          for j1 in acc'range loop
            if acc(j1) > 0 then
              for j2 in 0..acc(j1)-1 loop
                Clear(fac.cf);
                fac.cf := Create(integer(1));
                fac.dg(j1) := 1;
                sum := Create(fac);
                fac.dg(j1) := 0;
                Clear(fac.cf);
                fac.cf := -x(j1)(j2);
                Add(sum,fac);
                Mul(trm,sum);  -- trm := trm*(a(j1)-x(j1)(j2))
                Clear(sum);
              end loop;
            end if;
          end loop;
          Add(res,trm);
        end loop;
      else
        for i in reverse v.a..b loop
          acc(k) := i;
          Rec_Expand(v.w(i).all,k+1,b-i);
        end loop;
      end if;
    end Rec_Expand;
 
  begin
    fac.dg := new Standard_natural_Vectors.Vector'(1..n => 0);
    Rec_Expand(f,x'first,f.b);
    Clear(fac);
    return res;
  end Expand;

-- EVALUATORS WITHOUT HORNER :

  function Eval0 ( f : NesVec; x : VecVec; a : Vector )
                 return Complex_Number is

    res : Complex_Number := Create(integer(0));
    acc : Standard_Integer_Vectors.Vector(x'range) := (x'range => 0);

    procedure Rec_Eval ( v : in Nesvec; k,b : in integer32 ) is

    -- DESCRIPTION : recursive elaboration of the evaluation order.

    -- ON ENTRY :
    --   v         nested vector with generalized divided differences;
    --   k         current entry in the accumulating index acc;
    --   b         last entry in the triangular grid.

      term,dif : Complex_Number;

    begin
      if v.n = 1 then
        for i in reverse v.a..b loop
          acc(k) := i;
          Copy(v.v(i),term);
          for j1 in acc'range loop
            if acc(j1) > 0 then
              for j2 in 0..acc(j1)-1 loop
                dif := a(j1)-x(j1)(j2);
                Mul(term,dif);
                Clear(dif);
              end loop;
            end if;
          end loop;
          Add(res,term);
          Clear(term);
        end loop;
      else
        for i in reverse v.a..b loop
          acc(k) := i;
          Rec_Eval(v.w(i).all,k+1,b-i);
        end loop;
      end if;
    end Rec_Eval;

  begin
    Rec_Eval(f,x'first,f.b);
    return res;
  end Eval0;

  function Eval0 ( file : file_type; f : NesVec; x : VecVec; a : Vector )
                 return Complex_Number is

    res : Complex_Number := Create(integer(0));
    acc : Standard_Integer_Vectors.Vector(x'range) := (x'range => 0);

    procedure Rec_Eval ( v : in Nesvec; k,b : in integer32 ) is

    -- DESCRIPTION : recursive elaboration of the evaluation order.

    -- ON ENTRY :
    --   v         nested vector with generalized divided differences;
    --   k         current entry in the accumulating index acc;
    --   b         last entry in the triangular grid.

      term,dif : Complex_Number;

    begin
      if v.n = 1 then
        for i in reverse v.a..b loop
          acc(k) := i;
          put(file,"+"); put(file,acc,1);
          Copy(v.v(i),term);
          for j1 in acc'range loop
            if acc(j1) > 0 then
              for j2 in 0..acc(j1)-1 loop
                put(file,"*(pt("); put(file,j1,1);
                put(file,")-x(");  put(file,j1,1);
                put(file,")(");    put(file,j2,1);
                put(file,"))");
                dif := a(j1)-x(j1)(j2);
                Mul(term,dif);
              end loop;
            end if;
          end loop;
          new_line(file);
          Add(res,term);
          Clear(term);
        end loop;
      else
        for i in reverse v.a..b loop
          acc(k) := i;
          Rec_Eval(v.w(i).all,k+1,b-i);
        end loop;
      end if;
    end Rec_Eval;

  begin
    Rec_Eval(f,x'first,f.b);
    return res;
  end Eval0;

-- EVALUATORS WITH HORNER :

  function Eval ( f : NesVec; x : VecVec; a : Vector ) return Complex_Number is

    res,term : Complex_Number := Create(integer(0));
    acc : Standard_Integer_Vectors.Vector(x'range) := (x'range => 0);
    open : integer32 := 0;
    linear : boolean := true;

    procedure Rec_Eval ( v : in Nesvec; k,b,lj,rj : in integer32 ) is

    -- DESCRIPTION : recursive elaboration of the evaluation order.

    -- ON ENTRY :
    --   v         nested vector with generalized divided differences;
    --   k         current entry in the accumulating index acc;
    --   b         last entry in the triangular grid;
    --   lj        leftmost nonzero entry in the index acc is needed for
    --             the Horner nesting, which occurs when lj = rj;
    --   rj        rightmost nonzero entry in the index acc is needed for
    --             multiplication with (a(rj) - x(rj)(acc(rj)-1)).

      nlj,nrj : integer32;
      close : boolean := false;

    begin
      if v.n = 1 then
        for i in reverse v.a..b loop
          acc(k) := i;
          if i > 0 then            -- determine new leftmost and rightmost
            nrj := k;
            if lj = 0
             then nlj := k;
             else nlj := lj;
            end if;
          else
            nrj := rj;
            nlj := lj;
          end if;
          if nlj > 0 and nlj = nrj and open = 0 then       -- open bracket
            open := acc(nlj)-1;
            linear := false;
          else
            if nlj > 0 and nlj = nrj and open > 0 then
              close := true;                   -- close bracket
              open := open-1;
            else
              close := false;
            end if;
          end if;
          if nrj > 0 then
            if close then
              term := (term + v.v(i))*(a(nrj)-x(nrj)(acc(nrj)-1));
              res := res + term;
              term := Create(integer(0));
            else
              term := term + v.v(i)*(a(nrj)-x(nrj)(acc(nrj)-1));
            end if;
          elsif linear then
            res := term + v.v(i);
          else
            res := res + v.v(i);
          end if;
        end loop;
      else
        for i in reverse v.a..b loop
          acc(k) := i;
          if i > 0 then
            if lj = 0
             then Rec_Eval(v.w(i).all,k+1,b-i,k,k);
             else Rec_Eval(v.w(i).all,k+1,b-i,lj,k);
            end if;
          else
            Rec_Eval(v.w(i).all,k+1,b-i,lj,rj);
          end if;
        end loop;
      end if;
    end Rec_Eval;

  begin
   -- Rec_Eval(f,x'first,f.b,0,0);
    res := Eval0(f,x,a);
    return res;
  end Eval;

  function Eval ( file : file_type;
                  f : NesVec; x : VecVec; a : Vector ) return Complex_Number is

    res,term : Complex_Number := Create(integer(0));
    acc : Standard_Integer_Vectors.Vector(x'range) := (x'range => 0);
    open : integer32 := 0;
    linear : boolean := true;

    procedure Rec_Eval ( v : in Nesvec; k,b,lj,rj : in integer32 ) is

    -- DESCRIPTION : recursive elaboration of the evaluation order.

    -- ON ENTRY :
    --   v         nested vector with generalized divided differences;
    --   k         current entry in the accumulating index acc;
    --   b         last entry in the triangular grid;
    --   lj        leftmost nonzero entry in the index acc is needed for
    --             the Horner nesting, which occurs when lj = rj;
    --   rj        rightmost nonzero entry in the index acc is needed for
    --             multiplication with (a(rj) - x(rj)(acc(rj)-1)).

      nlj,nrj : integer32;
      close : boolean := false;

    begin
     -- put(file,v.n,1);
     -- put(file," "); put(file,v.a,1); 
     -- put(file," "); put(file,b,1); new_line(file);
      if v.n = 1 then
        for i in reverse v.a..b loop
          acc(k) := i;
          if i > 0 then             -- determine new leftmost and rightmost
            nrj := k;
            if lj = 0
             then nlj := k;
             else nlj := lj;
            end if;
          else
            nrj := rj;
            nlj := lj;
          end if;
          if nlj > 0 and nlj = nrj and open = 0 then       -- open bracket
            put(file,"+");
            open := acc(nlj)-1;
            linear := false;
            for ii in 1..open loop
              put(file,"(");
            end loop;
            put(file,acc,1);
          else
            put(file," +");
            put(file,acc,1);
            if nlj > 0 and nlj = nrj and open > 0 then
              put(file,")");                  -- close bracket
              close := true;
              open := open-1;
            else
              put(file," ");
              close := false;
            end if;
          end if;
          if nrj > 0 then
            put(file,"*(pt("); put(file,nrj,1); put(file,")-");
            put(file,"x("); put(file,nrj,1); put(file,")(");
            put(file,acc(nrj)-1,1); put_line(file,"))");
            if close then
              term := (term + v.v(i))*(a(nrj)-x(nrj)(acc(nrj)-1));
              res := res + term;
              term := Create(integer(0));
            else
              term := term + v.v(i)*(a(nrj)-x(nrj)(acc(nrj)-1));
            end if;
          else
            new_line(file);
            if linear
             then res := term + v.v(i);
             else res := res + v.v(i);
            end if;
          end if;
          -- put(file,"lj : "); put(file,nlj,1);
          -- put(file,"  rj : "); put(file,nrj,1); new_line(file);
        end loop;
      else
        for i in reverse v.a..b loop
          acc(k) := i;
          if i > 0 then
            if lj = 0
             then Rec_Eval(v.w(i).all,k+1,b-i,k,k);
             else Rec_Eval(v.w(i).all,k+1,b-i,lj,k);
            end if;
          else
            Rec_Eval(v.w(i).all,k+1,b-i,lj,rj);
          end if;
        end loop;
      end if;
    end Rec_Eval;

  begin
   -- Rec_Eval(f,x'first,f.b,0,0);
    res := Eval0(f,x,a);
    return res;
  end Eval;

-- DIAGNOSTICS :

  procedure Eval_Grid ( file : in file_type;
                        x : in VecVec; y,f : in NesVec;
                        maxerr : out Floating_Number ) is

    point : Vector(x'range);
    eva : Complex_Number;

    procedure Rec_Eval_Grid ( k : in integer32; v : in NesVec ) is

      dif : Complex_Number;
      err : Floating_Number;

    begin
      put(file,v.n,1);
      put(file," "); put(file,v.a,1);
      put(file," "); put(file,v.b,1); new_line(file);
      if v.n = 1 then
        for i in v.v'range loop
          put(file,"Sample : "); put(file,v.v(i)); new_line(file);
          point(k) := x(k)(i);
          eva := Eval(f,x,point);
          put(file,"Valued : "); put(file,eva); new_line(file);
          dif := v.v(i)-eva;
          err := AbsVal(dif);
          put(file,"-> Difference between Sample and Valued : ");
          put(file,err); new_line(file);
          if err > maxerr
           then Copy(err,maxerr);
          end if;
          Clear(err); Clear(eva); Clear(dif);
        end loop;
      else
        for i in v.w'range loop
          point(k) := x(k)(i);
          Rec_Eval_Grid(k+1,v.w(i).all);
        end loop;
      end if;
    end Rec_Eval_Grid;

  begin
    maxerr := Create(integer(0));
    Rec_Eval_Grid(x'first,y);
    put(file,"The maximal error on the grid : ");
    put(file,maxerr); new_line(file);
  end Eval_Grid;

  function Maximal_Error ( x : VecVec; y,f : NesVec ) return Floating_Number is

    res : Floating_Number := Create(integer(0));
    point : Vector(x'range);
    eva : Complex_Number;

    procedure Rec_Eval_Grid ( k : in integer32; v : in NesVec ) is

      err : Floating_Number;

    begin
      if v.n = 1 then
        for i in v.v'range loop
          point(k) := x(k)(i);
          eva := Eval(f,x,point);
          err := AbsVal(v.v(i)-eva);
          if err > res
           then Copy(err,res);
          end if;
          Clear(eva); Clear(err);
        end loop;
      else
        for i in v.w'range loop
          point(k) := x(k)(i);
          Rec_Eval_Grid(k+1,v.w(i).all);
        end loop;
      end if;
    end Rec_Eval_Grid;

  begin
    Rec_Eval_Grid(x'first,y);
    return res;
  end Maximal_Error;

end Multprec_Nvariate_Interpolators;
