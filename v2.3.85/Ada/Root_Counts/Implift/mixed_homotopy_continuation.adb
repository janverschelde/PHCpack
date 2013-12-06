with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Integer_Support_Functions;          use Integer_Support_Functions;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Integer32_Vectors_Utilities;        use Integer32_Vectors_Utilities;
with Lists_of_Vectors32_Utilities;       use Lists_of_Vectors32_Utilities;
with Arrays_of_Lists_Utilities;          use Arrays_of_Lists_Utilities;
with Standard_Integer32_Transformations;
 use Standard_Integer32_Transformations;
with Transforming_Solutions;             use Transforming_Solutions;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Volumes;
with Standard_Durand_Kerner;             use Standard_Durand_Kerner;

package body Mixed_Homotopy_Continuation is

-- INVARIANT CONDITION :
--   The procedures and functions in this package `mirror' corresponding
--   routines in the package Volumes.

-- AUXILIARIES :

  function Interchange2 ( p : Laur_Sys; index : integer32 ) return Laur_Sys is

  -- DESCRIPTION :
  --   Returns a polynomial system where the first equation is interchanged
  --   with the equation given by the index.

    res : Laur_Sys(p'range);

  begin
    if index = p'first then
      res := p;
    else
      res(res'first) := p(index);
      res(index) := p(p'first);
      res(res'first+1..index-1) := p(p'first+1..index-1);
      res(index+1..res'last) := p(index+1..p'last);
    end if;
    return res;
  end Interchange2;

  procedure Interchange2 ( adl : in out Array_of_Lists; index : integer32 ) is

  -- DESCRIPTION :
  --   Interchanges the first list with the list given by the index.

    tmp : List;

  begin
    if index /= adl'first then
      tmp := adl(adl'first);
      adl(adl'first) := adl(index);
      adl(index) := tmp;
    end if;
  end Interchange2;

  function Permute ( perm : Standard_Integer_Vectors.Vector; p : Laur_Sys )
                   return laur_Sys is

  -- DESCRIPTION :
  --   Returns a permuted polynomial system.

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := p(perm(i));
    end loop;
    return res;
  end Permute;

  function Initial_Degrees ( p : Poly ) return Degrees is

   -- DESCRIPTION :
   --   Returns the initial degrees of the polynomial p.

    init : Degrees;

    procedure Init_Term ( t : in Term; cont : out boolean ) is
    begin
      init := new Standard_Integer_Vectors.Vector'(t.dg.all);
      cont := false;
    end Init_Term;
    procedure Initial_Term is new Visiting_Iterator (Init_Term);

  begin
    Initial_Term(p);
    return init;
  end Initial_Degrees;

  procedure Binomial ( p : in Poly;
                       d : out Standard_Integer_Vectors.Link_to_Vector;
                       k : in out integer32; c : in out Complex_Number ) is

   -- DESCRIPTION :
   --   p consists of two terms, this procedure gets the degrees d and
   --   the constant c of the binomial equation.

    first : boolean := true;
    dd : Degrees;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is
    begin
      if first then
        dd := new Standard_Integer_Vectors.Vector'(t.dg.all);
        c := t.cf;
        first := false;
      else
        k := dd'first - 1;
        for i in dd'range loop
          dd(i) := dd(i) - t.dg(i);
	  if k < dd'first and then dd(i) /= 0
  	   then k := i;
          end if;
         end loop;
	 c := -t.cf/c;
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator (Scan_Term);

  begin
    Scan_Terms(p);
    d := new Standard_Integer_Vectors.Vector'(dd.all);
    Standard_Integer_Vectors.Clear(Standard_Integer_Vectors.Link_to_Vector(dd));
  end Binomial;

  procedure Normalize ( p : in Laur_Sys; dl : in out List; 
			wp : in out Laur_Sys; shifted : out boolean ) is

  -- DESCRIPTION :
  --   Makes sure that the first element of dl contains all zeroes.

  -- REQUIRED :
  --   The list dl is not empty.

  -- ON ENTRY :
  --   p           a Laurent polynomial system;
  --   dl          the power list of p(p'first).

  -- ON RETURN :
  --   dl          the first element of contains all zeroes,
  --                therefore dl has been shifted;
  --   wp          a Laurent polynomial system where
  --                dl is the power list of wp(wp'first);
  --   shifted     is true if dl has been shifted.

    use Standard_Integer_Vectors;

    first : constant Link_to_Vector := Head_Of(dl);
    nullvec : constant Vector(first'range) := (first'range => 0); 
    shiftvec : Link_to_Vector;

  begin
    if not Is_In(dl,nullvec) then
      shiftvec := Graded_Max(dl);
      Shift(dl,shiftvec);
      Clear(shiftvec);
      Copy(p,wp);
      for i in p'range loop
        Shift(wp(i));
      end loop;
    else
      wp := p;
      shifted := false;
    end if;
    Move_to_Front(dl,nullvec);
  end Normalize;

  function Re_Arrange ( p : poly ) return Poly is

  -- DESCRIPTION :
  --   Returns a polynomial whose terms are sorted
  --   in graded lexicographical ordening.
   
    res : Poly := Null_Poly;

    procedure Re_Arrange_Term ( t : in Term; cont : out boolean ) is
    begin
      Add(res,t);
      cont := true;
    end Re_Arrange_Term;
    procedure Re_Arrange_Terms is new Visiting_Iterator (Re_Arrange_Term);

  begin
    Re_Arrange_Terms(p);
    return res;
  end Re_Arrange;

  function Substitute ( p : Poly; v : Standard_Complex_Vectors.Vector;
                        k : integer32 ) return Poly is

   -- DESCRIPTION :
   --   Substitutes the values in v into the polynomial p,
   --   starting at the last unknowns of p.

    init : Degrees := Initial_Degrees(p);
    index : integer32 := init'last;
    res,tmp : Poly;

  begin
    if index = k then
      index := index - 1;
      res := Eval(p,v(v'last),index);
    else 
      res := Eval(p,v(v'last),index);
    end if;
    for i in reverse v'first..(v'last-1) loop
      index := index - 1;
      if index = k
       then index := index - 1;
      end if;
      tmp := Eval(res,v(i),index);
      Clear(res); Copy(tmp,res); Clear(tmp);
    end loop;
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(init));
    return res;
  end Substitute;

--  procedure Refine_and_Concat 
--             ( p : in Laur_Sys;
--               newsols,sols,sols_last : in out Solution_List ) is
--
--  -- DESCRIPTION :
--  --   This procedure refines the solutions of a given
--  --   polynomial system and adds them to the solution list.
--  --   This can be very useful for eliminating rounding errors
--  --   after transformating the solutions.
--
--    pp : Poly_Sys(p'range) := Laurent_to_Polynomial_System(p);
--    numb : natural := 0;
--
--  begin
--    Silent_Root_Refiner
--      (pp,newsols,10.0**(-12),10.0**(-12),10.0**(-8),numb,5,false);
--    Concat(sols,sols_last,newsols);
--    Clear(pp); Shallow_Clear(newsols);
--  end Refine_and_Concat;

  procedure Refine_and_Concat 
	       ( file : in file_type; p : in Laur_Sys;
	         newsols,sols,sols_last  : in out Solution_List ) is

  -- DESCRIPTION :
  --   This procedure refines the solutions of a given
  --   polynomial system and adds them to the solution list.
  --   This can be very useful for eliminating rounding errors
  --   after transformating the solutions.

    pp : Poly_Sys(p'range) := Laurent_to_Polynomial_System(p);
    numb : natural32 := 0;
    deflate : boolean := false;

  begin
    Reporting_Root_Refiner
      (file,pp,newsols,10.0**(-12),10.0**(-12),10.0**(-8),numb,5,deflate,false);
    Concat(sols,sols_last,newsols);
    Clear(pp); Shallow_Clear(newsols);
  end Refine_and_Concat;

  procedure Write_Direction
                  ( file : in file_type;
                    v : in Standard_Integer_Vectors.Link_to_Vector ) is
  begin
    put(file,"++++  considering direction "); put(file,v);
    put_line(file,"  ++++");
  end Write_Direction;

-- INTERMEDIATE LAYER :

  procedure Mixed_Continue
             ( file : in file_type; p : in Laur_Sys;
               k : in integer32; m : in Standard_Integer_Vectors.Vector;
               sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   This continuation routine computes a part of the solution list 
  --   of a Laurent polynomial system.

  -- ON ENTRY :
  --   file      a file where the intermediate results are written;
  --   p         the transformed Laurent polynomial system to be solved;
  --   k         the index;
  --   m         m(k) = p1(v), m(k/=l) = Maximal_Degree(p(l),k);
  --   sols      the start solutions.

  -- ON RETURN :
  --   sols      the computed solutions.

    h : Laur_Sys(p'range);
    hp : Poly_Sys(h'range);
    hpe : Eval_Poly_Sys(hp'range);
    j : Jaco_Mat(p'range,p'first..p'last+1);
    je : Eval_Jaco_Mat(j'range(1),j'range(2));

    function Construct_Homotopy
               ( p : Laur_Sys; k : integer32;
                 m : Standard_Integer_Vectors.Vector ) return Laur_Sys is

      res : Laur_Sys(p'range);
      ran : Complex_Number;
      re : boolean;
      zeroes : Degrees := new Standard_Integer_Vectors.Vector'(p'range => 0);

      function Construct_First_Polynomial
                 ( pp : Poly; max : integer32 ) return Poly is

        r : Poly := Null_Poly;

        procedure Construct_Term ( t : in Term; cont : out boolean ) is

          rt : Term;

        begin
          rt.cf := t.cf;
          rt.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last+1);
          rt.dg(t.dg'range) := t.dg.all;
          rt.dg(k) := -t.dg(k) + max;
          if Equal(t.dg,zeroes)
           then rt.dg(rt.dg'last) := 0;
                re := ( IMAG_PART(rt.cf) + 1.0 = 1.0 );
           else rt.dg(rt.dg'last) := rt.dg(k);
          end if;
          Add(r,rt);
          Clear(rt);
          cont := true;
        end Construct_Term;
        procedure Construct_Terms is new Visiting_Iterator(Construct_Term);

      begin
        Construct_Terms(pp);
        Standard_Integer_Vectors.Clear
          (Standard_Integer_Vectors.Link_to_Vector(zeroes));
        return r;
      end Construct_First_Polynomial;

      function Construct_Polynomial
                 ( pp : Poly; max : integer32) return Poly is

        r : Poly := Null_Poly;

        procedure Construct_Term ( t : in Term; cont : out boolean ) is

          rt : Term;

        begin
          rt.cf := t.cf;
          rt.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last+1);
          rt.dg(t.dg'range) := t.dg.all;
          rt.dg(k) := -t.dg(k) + max;
          rt.dg(rt.dg'last) := rt.dg(k);
          Add(r,rt);
          Clear(rt);
          cont := true;
        end Construct_Term;
        procedure Construct_Terms is new Visiting_Iterator(Construct_Term);

      begin
        Construct_Terms(pp);
        return r;
      end Construct_Polynomial;

    begin
      res(res'first) := Construct_First_Polynomial(p(p'first),m(m'first));
      for i in p'first+1..p'last loop
        res(i) := Construct_Polynomial(p(i),m(i));
      end loop;
      if re then
        for i in res'range loop
          ran := Random;
          Mul(res(i),ran);
        end loop;
      end if;
      return res;
    end Construct_Homotopy;

    function To_Be_Removed ( flag : in integer32 ) return boolean is
    begin
      return ( flag /= 1 );
    end To_Be_Removed;
    procedure Extract_Regular_Solutions is
      new Standard_Complex_Solutions.Delete(To_Be_Removed);

  begin
    h := Construct_Homotopy(p,k,m);               -- CONSTRUCTION OF HOMOTOPY
    hp := Laurent_to_Polynomial_System(h);
    hpe := Create(hp);
    j := Create(hp);
    je := Create(j);
    declare                                                   -- CONTINUATION 

      use Standard_Complex_Vectors;

      function Eval ( x : Vector; t : Complex_Number ) return Vector is

        xt : Vector(x'first..x'last+1);

      begin
        xt(x'range) := x;
        xt(xt'last) := t;
        return Eval(hpe,xt);
      end Eval;

      function dHt ( x : Vector; t : Complex_Number ) return Vector is

        xt : Vector(x'first..x'last+1);
        res : Vector(p'range);

      begin
        xt(x'range) := x;
        xt(xt'last) := t;
        for i in res'range loop
          res(i) := Eval(je(i,xt'last),xt);
        end loop;
        return res;
      end dHt;

      function dHx ( x : Vector; t : Complex_Number ) return Matrix is

        xt : Vector(x'first..x'last+1);
        m : Matrix(x'range,x'range);

      begin
        xt(x'range) := x;
        xt(xt'last) := t;
        for i in m'range(1) loop
          for j in m'range(1) loop
            m(i,j) := Eval(je(i,j),xt);
          end loop;
        end loop;
        return m;
      end dHx;

      procedure Invert ( k : in integer32; sols : in out Solution_List ) is

      -- DESCRIPTION :
      --   For all solutions s in sols : s.v(k) := 1/s.v(k).
      
      tmp : Solution_List := sols;

      begin
        while not Is_Null(tmp) loop
          declare
	    l : constant Link_to_Solution := Head_Of(tmp);
          begin
            l.v(k) := Create(1.0)/l.v(k);
            l.t := Create(0.0);
          end;
          tmp := Tail_Of(tmp);
        end loop;
      end Invert;

      procedure Cont is new Reporting_Continue(Max_Norm,Eval,dHt,dHx);

    begin
      Invert(k,sols);
      Cont(file,sols,false);
      Invert(k,sols);
      Extract_Regular_Solutions(sols);
    end;
    Clear(h); Clear(hp); Clear(hpe); Clear(j); Clear(je);
  end Mixed_Continue;

-- THE SOLVERS : 

  function One_Unknown_Solve
             ( p : Poly ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the solution vector of p, a polynomial in one unknown.
  --   p will be solved by using the method of Durand-Kerner.

    p1 : Poly := Re_Arrange(p);
    init : Degrees := Initial_Degrees(p1);
    index : constant integer32 := init'first;
    min : constant integer32 := -Minimal_Degree(p1,index);
    pv : Standard_Complex_Vectors.Vector(0..Degree(p1)+min);
    z,res : Standard_Complex_Vectors.Vector(1..pv'last);
    maxsteps : constant natural32 := 10;
    eps : constant double_float := 10.0**(-10);
    nb : integer32 := pv'last + 1;
    fail : boolean;

    procedure Store_Coeff ( t : in Term; cont : out boolean ) is
    begin
      nb := nb - 1;
      if t.dg(index) = (nb - min) then
        pv(nb) := t.cf;
      else
        for i in reverse nb..(t.dg(index)+min+1) loop
          pv(i) := Create(0.0);
        end loop;
        nb := t.dg(index) + min;
        pv(nb) := t.cf;
      end if;
      cont := true;
    end Store_Coeff;
    procedure Polynomial_To_Vector is new Visiting_Iterator (Store_Coeff);

  begin
    for i in pv'range loop               -- initialize coefficient vector
      pv(i) := Create(0.0);
    end loop;
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(init));
    Polynomial_To_Vector(p1);
    Clear(p1);
    for i in z'range loop
      z(i) := Random;
      res(i) := z(i);
    end loop;
   -- put_line("Applying Durand Kerner to "); put_line(pv);
    Silent_Durand_Kerner(pv,z,res,maxsteps,eps,natural32(nb),fail);
    return z;
  end One_Unknown_Solve;

  procedure One_Unknown_Solve ( p : in Poly; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   If p is a polynomial in one unknown,
  --   p can be solved efficiently by the application of Durand-Kerner.

    init : Degrees := Initial_Degrees(p);

    function Make_Solutions ( x : in Standard_Complex_Vectors.Vector )
                            return Solution_List is

      res,res_last : Solution_List;
      s : Solution(1);

    begin
      s.m := 1;
      s.t := Create(0.0);
      for i in x'range loop
	s.v(init'first) := x(i);
	Append(res,res_last,s);
      end loop;
      return res;
    end Make_Solutions;

  begin
    sols := Make_Solutions(One_Unknown_Solve(p));
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(init));
  end One_Unknown_Solve;

  procedure Two_Terms_Solve
               ( file : in file_type; p : in Laur_Sys;
                 tv : in Tree_of_Vectors; bkk : out natural32;
                 sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   The first polynomial of p consists of two terms.
  --   A binomial system can be solved efficiently by 
  --   transforming and using de Moivre's rule.

    d : Standard_Integer_Vectors.Link_to_Vector;
    c : Complex_Number := Create(0.0);
    k : integer32 := 0;
    sols_last : Solution_List;

  begin
   -- put_line(file,"Applying Two_Terms_Solve on "); put_line(file,p);
    Binomial(p(p'first),d,k,c);
   -- put(file,"After binomial : k = "); put(file,k,1); new_line(file);
   -- put(file,"the degrees :"); put(file,d.all); new_line(file);
   -- put(file,"c = "); put(file,c); new_line(file);
    if k < d'first then
      bkk := 0;
      Standard_Integer_Vectors.Clear(d); return;
    elsif ( c + Create(1.0) = Create(1.0) ) then
      bkk := 0;
      Standard_Integer_Vectors.Clear(d); return;
    elsif d(k) < 0 then
      Standard_Integer_Vectors.Min(d);
      c := Create(1.0)/c;
    end if;
    declare
      t : Transfo := Rotate(d,k);
      tmp_bkk : natural32 := 0;
    begin
      Apply(t,d);
      declare
        v : Standard_Complex_Vectors.Vector(1..d(k));
        tmp : Poly;
        rtp : Laur_Sys(p'first..(p'last-1));
        rtp_bkk : natural32;
        rtp_sols : Solution_List;
      begin
        for i in v'range loop
          v(i) := Root(c,natural32(d(k)),natural32(i));
          for j in rtp'range loop
            tmp := Transform(t,p(j+1));
            rtp(j) := Eval(tmp,v(i),k);
            Clear(tmp);
          end loop;
          Solve(file,rtp,tv,rtp_bkk,rtp_sols);
          Clear(rtp);
          tmp_bkk := tmp_bkk + rtp_bkk;
          Insert(v(i),k,rtp_sols);
          Transform(t,rtp_sols);
         --Concat(sols,sols_last,rtp_sols);
          Refine_and_Concat(file,p,rtp_sols,sols,sols_last);
        end loop;
      end;
      Clear(t);
      bkk := tmp_bkk;
    end;
    Standard_Integer_Vectors.Clear(d);
   -- put_line(file,"The solutions found : ");  put(file,sols);
  end Two_Terms_Solve;

  function Project_on_First_and_Solve
                  ( p : Poly; k : integer32; sols : Solution_List )
                  return Solution_List is

  -- ON ENTRY :
  --   p          a Laurent polynomial in n unknowns x1,..,xk,..,xn;
  --   sols       contains values for x1,..,xn, except xk.

   -- ON RETURN :
   --   a solution list for p, obtained after substition of the values
   --   for x1,..,xn into the polynomial p.

    tmp,res,res_last : Solution_List;
    init : Degrees := Initial_Degrees(p);

  begin
   -- put_line(file,"Calling Project_on_First_and_Solve");
   -- put_line(file," with polynomial : ");
   -- put(file,Laurent_Polynomial_to_Polynomial(p)); new_line(file);
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
	p1 : Poly := Substitute(p,Head_Of(tmp).v,k);
	sols1 : Solution_List;
      begin
       -- put(file,"k : "); put(file,k,1); new_line(file);
       -- put(file,"v : "); put(file,Head_Of(tmp).v,3,3,3); new_line(file);
       -- put_line(file,"After substitution : "); Write(file,p1);
	if Number_of_Terms(p1) < 2 then
	  null;
        elsif Number_of_Terms(p1) = 2 then
          declare
            d : Standard_Integer_Vectors.Link_to_Vector;
	    L : integer32 := 0;
	    c : Complex_Number := Create(0.0);
          begin
            Binomial(p1,d,L,c);
	    if L < init'first then
	      null;
            elsif ( c + Create(1.0) = Create(1.0) ) then
	      null;
            else
              if d(L) < 0
	       then d(L) := -d(L); c := Create(1.0)/c;
              end if;
              declare
	        v : Standard_Complex_Vectors.Vector(1..d(L));
              begin
	        for i in v'range loop
	          v(i) := Root(c,natural32(d(L)),natural32(i));
                end loop;
	        sols1 := Insert(v,k,Head_Of(tmp).all);
               -- put_line(file,"Sols1 after Binomial :");
               -- put(file,sols1);
                Concat(res,res_last,sols1);
                Shallow_Clear(sols1);
              end;
            end if;
            Standard_Integer_Vectors.Clear(d);
          end;
        else
          sols1 := Insert(One_Unknown_Solve(p1),k,Head_Of(tmp).all);
          Concat(res,res_last,sols1);
         -- put_line(file,"Sols1 :"); put(file,sols1);
          Shallow_Clear(sols1);
        end if;
        Clear(p1);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(init));
    return res;
  end Project_on_First_and_Solve;

  procedure Project_and_Solve
                ( file : in file_type; p : in Laur_Sys; k : in integer32;
                  m : in out Standard_Integer_Vectors.Vector;
                  nd : in node; bkk : out natural32;
                  sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves the projected start system of p along a direction v.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a Laurent polynomial system;
  --   k          entry in the degrees of p;
  --   m          m(m'first) equals Maximal_Support(p(p'first),v) > 0;
  --   nd         a node in the tree of useful directions.

  -- ON RETURN :
  --   m          m(m'first+1..m'last) contains the maximal degrees
  --               of the last n-1 equations of p in xk;
  --   bkk        the BKK bound of the projected system;
  --   sols       the solutions of the projected system.

    g_v : Laur_Sys(p'first..(p'last-1));
    bkk_g_v : natural32;
    sols_g_v : Solution_List;

  begin
   -- put_line(file,"Applying Project_and_Solve on"); Write(file,p);
    for i in g_v'range loop
      m(i+1) := Maximal_Degree(p(i+1),k);
      g_v(i) := Face(k,m(i+1),p(i+1));
      Reduce(k,g_v(i));
    end loop;
    if (nd.ltv = null) or else Is_Null(nd.ltv.all)
     then Solve(file,g_v,bkk_g_v,sols_g_v);
     else Solve(file,g_v,nd.ltv.all,bkk_g_v,sols_g_v);
    end if;
   -- put(file,"After Solve (without tv) bkk_g_v = "); put(file,bkk_g_v,1);
   -- new_line(file);
    declare
      p0 : Poly := Re_Arrange(p(p'first));
      p1 : Poly := Face(k,m(m'first),p0);
      cnst : Term;
    begin
      cnst.dg := new Standard_Integer_Vectors.Vector'(p'range => 0);
      if Coeff(p1,cnst.dg) = Create(0.0) then
        cnst.cf := Coeff(p0,cnst.dg);
        Add(p1,cnst);
      end if;
      Standard_Integer_Vectors.Clear
        (Standard_Integer_Vectors.Link_to_Vector(cnst.dg));
      Clear(p0);
      sols := Project_on_First_and_Solve(p1,k,sols_g_v);
     -- sols := Project_on_First_and_Solve(file,p1,k,sols_g_v);
      Set_Continuation_Parameter(sols,Create(0.0));
      Clear(p1);
    end;
    bkk := natural32(m(m'first))*bkk_g_v;
    Clear(sols_g_v);
    Clear(g_v);
   -- put_line(file,"The solutions found : "); put(file,sols);
  end Project_and_Solve;

  procedure Unmixed_Solve
                ( file : in file_type; p : in Laur_Sys; dl : in List;
                  tv : in Tree_of_Vectors;
                  bkk : out natural32; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves a Laurent polynomial system where all polytopes are the same.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a Laurent polynomial system;
  --   dl         the list of powers for p;
  --   tv         the tree of degrees containing the useful directions.

  -- ON RETURN :
  --   bkk        the bkk bound;
  --   sols       the list of solutions.

    sols_last : Solution_List;
    tmp_bkk : natural32 := 0;
    tmp : Tree_of_Vectors := tv;

  begin
    tmp_bkk := 0;
    tmp := tv;
    while not Is_Null(tmp) loop
      declare
        nd : constant node := Head_Of(tmp);
        v : constant Standard_Integer_Vectors.Link_to_Vector := nd.d;
        i : constant integer32 := Pivot(v);
        pv : constant integer32 := Maximal_Support(dl,v.all);
        t : Transfo := Build_Transfo(v,i);
        tp : Laur_Sys(p'range) := Transform(t,p);
        bkk_tp : natural32;
        sols_tp : Solution_List;
        max : Standard_Integer_Vectors.Vector(p'range);
      begin
        Write_Direction(file,v);
        max(max'first) := pv;
       -- if (nd.ltv = null) or else Is_Null(nd.ltv.all)
       --  then Projected_Solve(file,tp,i,max,bkk_tp,sols_tp);
       --  else Projected_Solve
       --           (file,tp,i,max,nd.ltv.all,bkk_tp,sols_tp);
       -- end if;
        Project_and_Solve(file,tp,i,max,nd,bkk_tp,sols_tp);
        Mixed_Continue(file,tp,i,max,sols_tp);
        tmp_bkk := tmp_bkk + bkk_tp;
        Transform(t,sols_tp);
        --Concat(sols,sols_last,sols_tp);
        Refine_and_Concat(file,p,sols_tp,sols,sols_last);
        Clear(t); Clear(tp);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    bkk := tmp_bkk;
  end Unmixed_Solve;

  procedure Unmixed_Solve
                ( file : in file_type; p : in Laur_Sys; dl : in List;
                  bkk : out natural32; sols : in out Solution_List ) is

    tv : Tree_of_Vectors;

  begin
    Volumes.Volume(natural32(p'last),dl,tv,bkk);
    Unmixed_Solve(file,p,dl,tv,bkk,sols);
    Clear(tv);
  end Unmixed_Solve;

  procedure Mixed_Solve 
               ( file : in file_type; p : in Laur_Sys;
                 adl : in out Array_of_Lists; tv : in Tree_of_Vectors;
                 bkk : out natural32; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Computes the solutions of the Laurent polynomial system p,
  --   where p has more than one equation.

  -- NOTE :
  --   This procedure mirrors the procedure Minkowski_Sum in the body
  --   of the package Volumes.

    tmp_bkk,len : natural32;
    tmp : Tree_of_Vectors;
    index : constant integer32 := Index2(adl);
    wp : Laur_Sys(p'range);
    sols_last : Solution_List;
    shifted : boolean;
    perm,mix : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Interchange2(adl,index);
    len := Length_Of(adl(adl'first));
   -- put_line(file,"Applying Mixed_Solve on"); Write(file,p);
    if len = 2 then
      wp := Interchange2(p,index);
      Two_Terms_Solve(file,wp,tv,bkk,sols);
    elsif len > 2 then  -- INITIALIZATION :
      Mixture(adl,perm,mix);
      wp := Permute(perm.all,p);
      declare
        zeroes : Degrees
               := new Standard_Integer_Vectors.Vector'(p'range => 0);
        tmpwpi : Poly;
      begin
        if Coeff(wp(wp'first),zeroes) = Create(0.0) then
          shifted := true;
         -- wp(wp'first) := Shift(p(p'first));
          Copy(p(index),tmpwpi); wp(wp'first) := tmpwpi;
          Shift(wp(wp'first));
        else 
          shifted := false;
        end if;
        Standard_Integer_Vectors.Clear
          (Standard_Integer_Vectors.Link_to_Vector(zeroes));
      end; -- MIXED HOMOTOPY CONTINUATION :
      tmp_bkk := 0;
      tmp := tv;
      while not Is_Null(tmp) loop
        declare
          nd : constant node := Head_Of(tmp);
          v : constant Standard_Integer_Vectors.Link_to_Vector := nd.d;
          k : constant integer32 := Pivot(v);
          pv : constant integer32 := Maximal_Support(wp(wp'first),v);
          t : Transfo := Build_Transfo(v,k);
          twp : Laur_Sys(wp'range) := Transform(t,wp);
          bkk_twp : natural32;
          sols_twp : Solution_List;
          m : Standard_Integer_Vectors.Vector(wp'range);
        begin
          Write_Direction(file,v);
          m(m'first) := pv;
         -- if (nd.ltv = null) or else Is_Null(nd.ltv.all)
         --  then Projected_Solve(file,twp,k,m,bkk_twp,sols_twp);
         --  else Projected_Solve
         --           (file,twp,k,m,nd.ltv.all,bkk_twp,sols_twp);
         -- end if;
          Project_and_Solve(file,twp,k,m,nd,bkk_twp,sols_twp);
          Mixed_Continue(file,twp,k,m,sols_twp);
          tmp_bkk := tmp_bkk + bkk_twp;
          Transform(t,sols_twp);
         --Concat(sols,sols_last,sols_twp);
          Refine_and_Concat(file,wp,sols_twp,sols,sols_last);
          Clear(t); Clear(twp);
        end;
        tmp := Tail_Of(tmp);
      end loop;
      bkk := tmp_bkk;
      Standard_Integer_Vectors.Clear(perm);
      Standard_Integer_Vectors.Clear(mix);
      if shifted
       then Clear(wp(wp'first));
      end if;
    else -- len < 2
      bkk := 0;
    end if;
  end Mixed_Solve;

  procedure Mixed_Solve
               ( file : in file_type; p : in Laur_Sys;
                 adl : in out Array_of_Lists; bkk : out natural32;
                 sols : in out Solution_List ) is

    tv : Tree_of_Vectors;

  begin
    Volumes.Mixed_Volume(natural32(adl'last),adl,tv,bkk);
    Mixed_Solve(file,p,adl,tv,bkk,sols);
    Clear(tv);
  end Mixed_Solve;

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    bkk : out natural32; sols : in out Solution_List ) is

    al : Array_of_Lists(p'range) := Create(p);
    tv : Tree_of_Vectors;

  begin
    Volumes.Mixed_Volume(natural32(p'last),al,tv,bkk);
    Solve(file,p,tv,bkk,sols);
    Deep_Clear(al); Clear(tv);
  end Solve;

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    tv : in Tree_of_Vectors; bkk : out natural32;
                    sols : in out Solution_List ) is

  -- NOTE :
  --   This procedure mirrors the procedure Volumes.Mixed_Volume,
  --   with a tree of useful directions on entry.

  begin
    if p'first = p'last
     then One_Unknown_Solve(p(p'first),sols);
          bkk := Length_Of(sols);
     else --if Is_Fewnomial_System(p)
          -- then
          --   declare
          --     fail : boolean;
          --   begin
          --     Fewnomials.Solve(p,sols,fail);
          --     if fail
          --      then bkk := 0;  Clear(sols);
          --      else bkk := Length_Of(sols);
          --     end if;
          --   end;
          -- else
      declare
        adl : Array_of_Lists(p'range) := Create(p);
      begin
        if All_Equal(adl) then
          for i in (adl'first+1)..adl'last loop
            Deep_Clear(adl(i));
          end loop;
          declare
            wp : Laur_Sys(p'range);
            shifted : boolean;
          begin
            Normalize(p,adl(adl'first),wp,shifted);
            if Is_Null(tv)
             then Unmixed_Solve(file,wp,adl(adl'first),bkk,sols);
             else Unmixed_Solve(file,wp,adl(adl'first),tv,bkk,sols);
            end if;
            if shifted
             then Clear(wp);
            end if;
          end;
        elsif Is_Null(tv) then
          Mixed_Solve(file,p,adl,bkk,sols);
        else
          Mixed_Solve(file,p,adl,tv,bkk,sols);
        end if;
      end;
    end if;
  end Solve;

end Mixed_Homotopy_Continuation;
