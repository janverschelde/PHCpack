with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
--with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
--with Standard_Complex_Substitutors;      use Standard_Complex_Substitutors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Simpomial_Solvers;
with Polyhedral_Coefficient_Homotopies;  use Polyhedral_Coefficient_Homotopies;
with BKK_Bound_Computations;             use BKK_Bound_Computations;

package body Integer_Polyhedral_Continuation is

--  procedure Write ( file : in file_type;
--                    c : in Standard_Complex_VecVecs.VecVec ) is
--  begin
--    for i in c'range loop
--      put(file,i,1); put(file," : ");
--      for j in c(i)'range loop
--        if REAL_PART(c(i)(j)) = 0.0
--         then put(file," 0");
--         else put(file,c(i)(j),2,3,0);
--        end if;
--      end loop;
--      new_line(file);
--    end loop;
--  end Write;

 -- pragma inline(Eval);

-- HOMOTOPY CONSTRUCTOR :

  function Select_Subsystem ( p : Laur_Sys; mix : Vector; mic : Mixed_Cell )
                            return Laur_Sys is

    res : Laur_Sys(p'range);
    cnt : integer32 := 0;

  begin
    for k in mix'range loop
      for l in 1..mix(k) loop
        cnt := cnt + 1;
        res(cnt) := Select_Terms(p(cnt),mic.pts(k));
      end loop;
    end loop;
    return res;
  end Select_Subsystem;

  function Construct_Homotopy ( p : Laur_Sys; normal : Vector )
                              return Laur_Sys is

  -- DESCRIPTION :
  --   Given a Laurent polynomial system of dimension n*(n+1) and a
  --   normal, a homotopy will be constructed, with t = x(n+1)
  --   and so that the support of the start system corresponds with
  --   all points which give the minimal product with the normal.

    res : Laur_Sys(p'range);
    n : constant integer32 := p'length;
    use Standard_Complex_Laurentials;
  
    function Construct_Polynomial ( p : Poly; v : Vector ) return Poly is

      res : Poly := Null_Poly;

      procedure Construct_Term ( t : in Term; cont : out boolean ) is

        rt : term;

      begin
        rt.cf := t.cf;
        rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
        rt.dg(n+1) := t.dg.all*v;
        Add(res,rt);
        Clear(rt);
        cont := true;
      end Construct_Term;
      procedure Construct_Terms is 
        new Standard_Complex_Laurentials.Visiting_Iterator(Construct_Term);

    begin
      Construct_Terms(p);
      return res;
    end Construct_Polynomial;

  begin
   -- SUBSTITUTIONS :
    for k in p'range loop
      res(k) := Construct_Polynomial(p(k),normal);
    end loop;
   -- SHIFT :
    for k in res'range loop
      declare
        d : constant integer32 := Minimal_Degree(res(k),n+1);
        t : Term;
      begin
        t.cf := Create(1.0);
        t.dg := new Standard_Integer_Vectors.Vector'(1..n+1 => 0);
        t.dg(n+1) := -d;
        Mul(res(k),t);
        Clear(t);
      end;
    end loop;
    return res;
  end Construct_Homotopy;

--  function Determine_Power ( n : natural; h : Laur_Sys ) return positive is
--
--  -- DESCRIPTION :
--  --   Returns the smallest power of the last unknown,
--  --   over all polynomials in h.
--
--    res : positive := 1;
--    first : boolean := true;
--    d : integer;
--    use Standard_Complex_Laurentials;
--
--    procedure Scan ( t : in Term; cont : out boolean ) is
--    begin
--      if (t.dg(n+1) > 0) and then (t.dg(n+1) < d)
--       then d := t.dg(n+1);
--      end if;
--      cont := true;
--    end Scan;
--    procedure Search_Positive_Minimum is new Visiting_Iterator(Scan);
--
--  begin
--    for i in h'range loop
--      d := Maximal_Degree(h(i),n+1);
--      if d > 0
--       then Search_Positive_Minimum(h(i));
--            if first
--             then res := d;
--                  first := false;
--             elsif d < res
--                 then res := d;
--            end if;
--      end if;
--      exit when (d=1);
--    end loop;
--    if res = 1
--     then return res;
--     else return 2;
--    end if;
--  end Determine_Power;

  procedure Extract_Regular ( sols : in out Solution_List ) is

    function To_Be_Removed ( flag : in integer32 ) return boolean is
    begin
      return ( flag /= 1 );
    end To_Be_Removed;
    procedure Extract_Regular_Solutions is
      new Standard_Complex_Solutions.Delete(To_Be_Removed);

  begin
    Extract_Regular_Solutions(sols);
  end Extract_Regular;

--  procedure Refine ( file : in file_type; p : in Laur_Sys;
--                     sols : in out Solution_List ) is
--
--  -- DESCRIPTION :
--  --   Given a polyhedral homotopy p and a list of solution for t=1,
--  --   this list of solutions will be refined.
--
--    pp : Poly_Sys(p'range) := Laurent_to_Polynomial_System(p);
--    n : constant integer32 := p'length;
--   -- eps : constant double_float := 10.0**(-12);
--   -- tolsing : constant double_float := 10.0**(-8);
--   -- max : constant integer32 := 3;
--   -- numb : natural := 0;
--
--  begin
--    pp := Laurent_to_Polynomial_System(p);
--    Substitute(n+1,Create(1.0),pp);
--   -- Reporting_Root_Refiner(file,pp,sols,eps,eps,tolsing,numb,max,false);
--    Clear(pp); Extract_Regular(sols);
--  end Refine;

-- FIRST LAYER OF CONTINUATION ROUTINES :

  procedure Mixed_Continuation
                 ( p : in Laur_Sys; normal : in Vector;
                   sols : in out Solution_List ) is

    h   : Laur_Sys(p'range) := Construct_Homotopy(p,normal);
    hpe : Eval_Laur_Sys(h'range) := Create(h);
    j  : Jaco_Mat(h'range,h'first..h'last+1) := Create(h);
    je : Eval_Jaco_Mat(j'range(1),j'range(2)) := Create(j);

    function Eval ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      return Eval(hpe,xt);
    end Eval;

    function dHt ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(p'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range loop
        res(i) := Eval(je(i,xt'last),xt);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return matrix is

      m : Matrix(x'range,x'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

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

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Silent_Continue(Max_Norm,Eval,dHt,dHx);

  begin
   -- Continuation_Parameters.power_of_t := Determine_Power(h'length,h);
    Laur_Cont(sols,false);
    Clear(h); Clear(hpe); Clear(j); Clear(je);
    Extract_Regular(sols);
  end Mixed_Continuation;

  procedure Mixed_Continuation
                 ( file : in file_type; p : in Laur_Sys;
                   normal : in Vector; sols : in out Solution_List ) is

    h   : Laur_Sys(p'range) := Construct_Homotopy(p,normal);
    hpe : Eval_Laur_Sys(h'range) := Create(h);
    j  : Jaco_Mat(h'range,h'first..h'last+1) := Create(h);
    je : Eval_Jaco_Mat(j'range(1),j'range(2)) := Create(j);

    function Eval ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      return Eval(hpe,xt);
    end Eval;

    function dHt ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(p'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range loop
        res(i) := Eval(je(i,xt'last),xt);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return Matrix is

      m : Matrix(x'range,x'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

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

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Reporting_Continue(Max_Norm,Eval,dHt,dHx);

  begin
   -- Continuation_Parameters.power_of_t := Determine_Power(h'length,h);
    Laur_Cont(file,sols,false);
    Clear(h); Clear(hpe); Clear(j); Clear(je);
    Extract_Regular(sols);
  end Mixed_Continuation;

  procedure Mixed_Continuation
                ( mix : in Standard_Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  normal : in Vector; sols : in out Solution_List ) is

    pow : Standard_Integer_VecVecs.VecVec(c'range)
        := Power_Transform(e,lifted,mix,normal);
   -- scapow : Standard_Floating_VecVecs.VecVec(c'range) := Scale(pow);
    ctm : Standard_Complex_VecVecs.VecVec(c'range);

    function Eval ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is
    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      return Eval(h,ctm,x);
    end Eval;

    function dHt ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(h'range);
      xtl : constant integer32 := x'last+1;

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for i in res'range loop
        res(i) := Eval(j(i,xtl),m(i,xtl).all,ctm(i).all,x);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return matrix is

      mt : Matrix(x'range,x'range);

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for k in mt'range(1) loop
        for l in mt'range(2) loop
          mt(k,l) := Eval(j(k,l),m(k,l).all,ctm(k).all,x);
        end loop;
      end loop;
      return mt;
    end dHx;

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Silent_Continue(Max_Norm,Eval,dHt,dHx);

  begin
    for i in c'range loop
      ctm(i) := new Standard_Complex_Vectors.Vector(c(i)'range);
    end loop;
    Laur_Cont(sols,false);
    Standard_Integer_VecVecs.Clear(pow);
   -- Standard_Floating_VecVecs.Clear(scapow);
    Standard_Complex_VecVecs.Clear(ctm);
    Extract_Regular(sols);
  end Mixed_Continuation;

  procedure Mixed_Continuation
                ( file : in file_type; mix : in Standard_Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  normal : in Vector; sols : in out Solution_List ) is

    pow : Standard_Integer_VecVecs.VecVec(c'range) 
        := Power_Transform(e,lifted,mix,normal);
   -- scapow : Standard_Floating_VecVecs.VecVec(c'range) := Scale(pow);
    ctm : Standard_Complex_VecVecs.VecVec(c'range);

    function Eval ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is
    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      return Eval(h,ctm,x);
    end Eval;

    function dHt ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(h'range);
      xtl : constant integer32 := x'last+1;

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for i in res'range loop
        res(i) := Eval(j(i,xtl),m(i,xtl).all,ctm(i).all,x);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                 return Matrix is

      mt : Matrix(x'range,x'range);

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for k in m'range(1) loop
        for l in m'range(1) loop
          mt(k,l) := Eval(j(k,l),m(k,l).all,ctm(k).all,x);
        end loop;
      end loop;
      return mt;
    end dHx;

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Reporting_Continue(Max_Norm,Eval,dHt,dHx);

  begin
   -- put(file,"The normal : "); put(file,normal); new_line(file);
   -- put_line(file,"The exponent vector : ");
   -- for i in pow'range loop
   --   put(file,pow(i)); new_line(file);
   -- end loop;
   -- put_line(file,"The coefficient vector : "); Write(file,c);
    for i in c'range loop
      ctm(i) := new Standard_Complex_Vectors.Vector(c(i)'range);
    end loop;
    Laur_Cont(file,sols,false);
    Standard_Integer_VecVecs.Clear(pow);
   -- Standard_Floating_VecVecs.Clear(scapow);
    Standard_Complex_VecVecs.Clear(ctm);
    Extract_Regular(sols);
  end Mixed_Continuation;

-- UTILITIES FOR SECOND LAYER :

  function Sub_Lifting ( q : Laur_Sys; mix : Vector; mic : Mixed_Cell )
                       return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the lifting used to subdivide the cell.

    res : Array_of_Lists(mix'range);
    sup : Array_of_Lists(q'range);
    n : constant integer32 := q'last;
    cnt : integer32 := sup'first;

  begin
    for i in mic.pts'range loop
      sup(cnt) := Reduce(mic.pts(i),q'last+1);
      for j in 1..(mix(i)-1) loop
        Copy(sup(cnt),sup(cnt+j));
      end loop;
      cnt := cnt + mix(i);
    end loop;
    res := Induced_Lifting(n,mix,sup,mic.sub.all);
    Deep_Clear(sup);
    return res;
  end Sub_Lifting;

  procedure Refined_Mixed_Solve
               ( q : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                 qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves a subsystem q using the refinement of the cell mic.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lifq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);

  begin
    Mixed_Solve(lifq,mix,mic.sub.all,qsols);
    Deep_Clear(lif); Clear(lifq);
  end Refined_Mixed_Solve;

  procedure Refined_Mixed_Solve
               ( file : in file_type;
                 q : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                 qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves a subsystem q using the refinement of the cell mic.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lifq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);

  begin
    Mixed_Solve(file,lifq,mix,mic.sub.all,qsols);
    Deep_Clear(lif); Clear(lifq);
  end Refined_Mixed_Solve;

  function Sub_Polyhedral_Homotopy
               ( l : List; e : Standard_Integer_VecVecs.VecVec;
                 c : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   For every vector in e that does not belong to l, the corresponding
  --   index in c will be set to zero, otherwise it is copied to the result.

    res : Standard_Complex_Vectors.Vector(c'range);
    found : boolean;
    lif : integer32;

  begin
    for i in e'range loop
      Search_Lifting(l,e(i).all,found,lif);
      if not found
       then res(i) := Create(0.0);
       else res(i) := c(i);
      end if;
    end loop;
    return res;
  end Sub_Polyhedral_Homotopy;

  function Sub_Polyhedral_Homotopy
               ( mix : Vector; mic : Mixed_Cell;
                 e : Exponent_Vectors_Array;
                 c : Standard_Complex_VecVecs.VecVec )
               return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Given a subsystem q of p and the coefficient vector of p, the
  --   vector on return will have only nonzero entries for coefficients
  --   that belong to q.

    res : Standard_Complex_VecVecs.VecVec(c'range);
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      ind := ind+1;
      declare
        cri : constant Standard_Complex_Vectors.Vector
            := Sub_Polyhedral_Homotopy(mic.pts(i),e(ind).all,c(ind).all);
      begin
        res(ind) := new Standard_Complex_Vectors.Vector'(cri);
        for j in 1..(mix(i)-1) loop
          declare
            crj : constant Standard_Complex_Vectors.Vector
                := Sub_Polyhedral_Homotopy(mic.pts(i),e(i+j).all,c(i+j).all);
          begin
            ind := ind+1;
            res(ind) := new Standard_Complex_Vectors.Vector'(crj);
          end;
        end loop;
      end;
    end loop;
    return res;
  end Sub_Polyhedral_Homotopy;

  procedure Refined_Mixed_Solve
                ( q : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Polyhedral coeffient-homotopy for subsystem q.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);
    cq : Standard_Complex_VecVecs.VecVec(c'range)
       := Sub_Polyhedral_Homotopy(mix,mic,e,c);

  begin
    Mixed_Solve(lq,lif,h,cq,e,j,m,mix,mic.sub.all,qsols);
    Standard_Complex_VecVecs.Clear(cq); Deep_Clear(lif); Clear(lq);
  end Refined_Mixed_Solve;

  procedure Refined_Mixed_Solve
                ( file : in file_type; q : in Laur_Sys; mix : in Vector;
                  mic : in Mixed_Cell; h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Polyhedral coeffient-homotopy for subsystem q.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);
    cq : Standard_Complex_VecVecs.VecVec(c'range)
       := Sub_Polyhedral_Homotopy(mix,mic,e,c);

  begin
    Mixed_Solve(file,lq,lif,h,cq,e,j,m,mix,mic.sub.all,qsols);
    Standard_Complex_VecVecs.Clear(cq); Deep_Clear(lif); Clear(lq);
  end Refined_Mixed_Solve;

-- TARGET ROUTINES FOR SECOND LAYER :

  procedure Mixed_Solve 
               ( p : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                 sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    fail,zero_y : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    if mic.sub = null then 
      Standard_Simpomial_Solvers.Solve(sq,tol_zero,qsols,fail,zero_y);
    else
      fail := true;
    end if;
    if fail then
      if mic.sub = null
       then pq := Laurent_to_Polynomial_System(sq);
            qsols := Solve_by_Static_Lifting(pq); Clear(pq);
       else Refined_Mixed_Solve(q,mix,mic,qsols);
      end if;
      Set_Continuation_Parameter(qsols,Create(0.0));
    end if;
    len := Length_Of(qsols);
    if len > 0
     then Mixed_Continuation(p,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

  procedure Mixed_Solve 
               ( file : in file_type; p : in Laur_Sys;
                 mix : in Vector; mic : in Mixed_Cell;
                 sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    fail,zero_y : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    if mic.sub = null then
      Standard_Simpomial_Solvers.Solve(sq,tol_zero,qsols,fail,zero_y);
    else 
      fail := true;
    end if;
    if not fail then
      put_line(file,"It is a simplex system.");
    else
      put_line(file,"No simplex system.");
      if mic.sub = null
       then put_line(file,"Calling the black box solver.");
            pq := Laurent_to_Polynomial_System(sq);
            qsols := Solve_by_Static_Lifting(file,pq); Clear(pq);
       else put_line(file,"Using the refinement of the cell.");
            Refined_Mixed_Solve(file,q,mix,mic,qsols);
      end if;
      Set_Continuation_Parameter(qsols,Create(0.0));
    end if;
    len := Length_Of(qsols);
    put(file,len,1); put_line(file," solutions found.");
    if len > 0
     then Mixed_Continuation(file,p,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    fail,zero_y : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    if mic.sub = null then
      Standard_Simpomial_Solvers.Solve(sq,tol_zero,qsols,fail,zero_y);
    else
      fail := true;
    end if;
    if fail then
      if mic.sub = null
       then pq := Laurent_to_Polynomial_System(sq);
            qsols := Solve_by_Static_Lifting(pq); Clear(pq);
       else Refined_Mixed_Solve(q,mix,mic,h,c,e,j,m,qsols);
      end if;
      Set_Continuation_Parameter(qsols,Create(0.0));
    end if;
    len := Length_Of(qsols);
    if len > 0
     then Mixed_Continuation(mix,lifted,h,c,e,j,m,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    fail,zero_y : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    if mic.sub = null then
      Standard_Simpomial_Solvers.Solve(sq,tol_zero,qsols,fail,zero_y);
    else
      fail := true;
    end if;
    if not fail then
      put_line(file,"It is a simplex system.");
    else
      put_line(file,"No simplex system.");
      if mic.sub = null
       then put_line(file,"Calling the black box solver.");
            pq := Laurent_to_Polynomial_System(sq);
            qsols := Solve_by_Static_Lifting(file,pq); Clear(pq);
       else put_line(file,"Using the refinement of the cell.");
            Refined_Mixed_Solve(file,q,mix,mic,h,c,e,j,m,qsols);
      end if;
      Set_Continuation_Parameter(qsols,Create(0.0));
    end if;
    len := Length_Of(qsols);
    put(file,len,1); put_line(file," solutions found.");
    if len > 0
     then Mixed_Continuation(file,mix,lifted,h,c,e,j,m,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

-- THIRD LAYER :

  procedure Mixed_Solve
               ( p : in Laur_Sys;
                 mix : in Vector; mixsub : in Mixed_Subdivision;
                 sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Mixed_Solve(p,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

  procedure Mixed_Solve 
               ( file : in file_type; p : in Laur_Sys;
                 mix : in Vector; mixsub : in Mixed_Subdivision;
                 sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;
    cnt : integer32 := 0;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      cnt := cnt + 1;
      new_line(file);
      put(file,"*** PROCESSING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Mixed_Solve(file,p,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Mixed_Solve(p,lifted,h,c,e,j,m,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;
    cnt : integer32 := 0;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      cnt := cnt + 1;
      new_line(file);
      put(file,"*** PROCESSING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Mixed_Solve(file,p,lifted,h,c,e,j,m,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

end Integer_Polyhedral_Continuation;
