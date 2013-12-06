--with text_io,integer_io;                use text_io,integer_io;
--with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
--with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
--with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
--with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;

with Basis_Exchanges;                   use Basis_Exchanges;

package body Linear_Minimization is

-- AUXILIARIES FOR PIVOTING :

  function Is_In ( n : integer32;
                   v : Standard_Integer_Vectors.Vector; k : integer32 )
                 return boolean is
  
  -- DESCRIPTION :
  --   Returns true if there is an index in 1..n such that v(index) = k.

  begin
    for i in 1..n loop
      if v(i) = k
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Complement ( n,m : integer32;
                        active : Standard_Integer_Vectors.Vector )
                      return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of length m-n which contains complementary
  --   indices than the ones that are in active.

    res : Standard_Integer_Vectors.Vector(1..m-n);
    ind : integer32 := 0;

  begin
    for i in 1..m loop
      if not Is_In(n,active,i)
       then ind := ind+1; res(ind) := i;
      end if;
    end loop;
    return res;
  end Complement;

-- PIVOTING ROUTINES :

  function Leave_Variable
             ( cost : Standard_Floating_Vectors.Vector;
               binv : Standard_Floating_Matrices.Matrix;
               tol : double_float ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of the variable that will become passive.
  --   This is an index of the vector in the columns of binv that
  --   makes a most negative inner product with the cost vector.
  --   This ensures to have the largest cost reduction,
  --   because we consider the update of solution x1 = x0 + s*d(i);
  --   with d the i-th column in binv.
  --   If c*d(i) >= 0 for all i, then the solution is already optimal,
  --   and the integer on return will be 0.

  -- REQUIRED : cost'range = binv'range(1).

    ind : integer32 := binv'first(1);               -- index of min column
    min : double_float := 0.0;                      -- min inner product
    wrk : double_float;                             -- working variable

  begin
    for i in cost'range loop                        -- initialization
      min := min + binv(ind,i)*cost(i);             -- first column of binv
    end loop;
    for i in binv'first(1)+1..binv'last(1) loop     -- iterate columns of binv
      wrk := 0.0;
      for j in cost'range loop                      -- compute inner product
        wrk := wrk + binv(i,j)*cost(j);
      end loop;
      if wrk < min                                  -- ensure invariance
       then min := wrk; ind := i;
      end if;
    end loop;
    if min < -tol
     then return ind;
     else return 0;                                 -- optimality reached
    end if;
  end Leave_Variable;

  procedure Enter_Variable
              ( cff,binv : in Standard_Floating_Matrices.Matrix;
                rhs,x : in Standard_Floating_Vectors.Vector;
                passive : in Standard_Integer_Vectors.Vector;
                leave_var : in integer32; tol : in double_float;
                enter_var : out integer32; s : out double_float;
                degenerate : out boolean ) is

  -- DESCRIPTION :
  --   Returns the index of the contraint that will become active.
  --   This index i belongs to the range of passive for which
  --   (cff(j)*x - rhs(j))/cff(j)*binv(leave_var) is minimal and positive,
  --   with j = passive(i).  This ensures feasibility.
  --   If all denominators in the quotient above are zero,
  --   then there is no feasibility contraint and the solution is unbounded.
  --   In that case, the index enter_var on return is passive'last+1.
  --   If there is one denominator in the quotient above that is positive,
  --   then s < 0 and there is no improvement possible.
  --   In that case, the index enter_var on return is passive'first-1;
  --   In case the entering variable corresponds to a degenerate vertex,
  --   i.e.: a vertex satisfying more than n equalities, then the out
  --   variable degenerate is true on return.

  -- REQUIRED : cff'range(1) = x'range

    ind : integer32 := passive'last+1;     -- initialize as unbounded
    pai : integer32;
    cycind : integer32 := 0;
    min,val,wrk : double_float := 0.0;

  begin
    for i in passive'range loop
      pai := passive(i);
      wrk := rhs(pai);
      val := 0.0;
      for j in x'range loop
        wrk := wrk - x(j)*cff(j,pai);
        val := val + binv(leave_var,j)*cff(j,pai);
      end loop;
      if cycind = 0 then
        if abs(val) < tol and abs(wrk) < tol
         then cycind := i;
        end if;
      end if;
      if val < -tol then                 -- otherwise no feasibility contraint
        wrk := wrk/val;
        if ((ind > passive'last) or else (wrk < min))
         then min := wrk; ind := i;
        end if;
      end if;
      if min < -tol                      -- this should not happen :
       then ind := passive'first-1;      -- pivoting violates feasibility
      end if;
      exit when (ind = passive'first-1);
    end loop;
    if ((ind >= passive'first) and then (ind <= passive'last)) then
      s := min;
      degenerate := false;
    else 
      s := 0.0;                      -- unbounded or optimality reached 
      if cycind /= 0 then
        ind := cycind;               -- allow pivoting without decrease
        degenerate := true;
      else
        degenerate := false;
      end if;
    end if;
    enter_var := ind;
  end Enter_Variable;

  procedure Pivoting
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  cost,rhs : in Standard_Floating_Vectors.Vector;
                  binv : in out Standard_Floating_Matrices.Matrix;
                  x : in out Standard_Floating_Vectors.Vector;
                  active,passive : in out Standard_Integer_Vectors.Vector;
                  tol : in double_float; stop,infty,fail : out boolean ) is

  -- DESCRIPTION :
  --   Does one pivoting step with the current solution.

  -- ON ENTRY :
  --   n          number of variables <= number of rows in cff;
  --   m          number of constraints <= number of columns in cff;
  --   cff        coefficient matrix of the constraints;
  --   cost       cost vector to minimize;
  --   rhs        right-hand side vector;
  --   binv       inverse of the current basis vector;
  --   x          current feasible solution;
  --   active     indices to active constraints;
  --   passive    complement of active, constraints that are not active;
  --   tol        tolerance on the absolute value of a nonzero number.

  -- ON RETURN :
  --   binv       updated inverse of the basis;
  --   x          updated feasible solution;
  --   active     updated active constraints;
  --   passive    updated passive constraints;
  --   stop       pivoting stopped because either optimal or unbounded;
  --   infty      minimization problem is unbounded;
  --   fail       update of the inverse basis matrix failed.

    leave_ind,enter_ind,ind : integer32;
    ds : double_float;
    dg : boolean;

  begin
    stop := false;
    infty := false;
    leave_ind := Leave_Variable(cost,binv,tol);
    if (leave_ind < binv'first(1)) then
      stop := true;                    -- solution is already optimal
    else
      Enter_Variable(cff,binv,rhs,x,passive,leave_ind,tol,enter_ind,ds,dg);
      if enter_ind > passive'last
       then infty := true;
       else infty := false;
      end if;
      if ((enter_ind < passive'first) or (enter_ind > passive'last)) then
        stop := true;
      else
        ind := passive(enter_ind);        -- ind for update of binv
        passive(enter_ind) := active(leave_ind);
        active(leave_ind) := ind;
        for i in x'range loop
          x(i) := x(i) + ds*binv(leave_ind,i);
        end loop;
        Update(n,m,binv,cff,active,leave_ind,ind,tol,fail);
      end if;
    end if;
  end Pivoting;

-- EVALUATORS :

  function Eval ( n,k : integer32; 
                  cff : Standard_Floating_Matrices.Matrix;
                  x : Standard_Floating_Vectors.Vector )
                return double_float is

    res : double_float := 0.0;

  begin
    for i in 1..n loop
      res := res + x(i)*cff(i,k);
    end loop;
    return res;
  end Eval;

  function Eval ( n : integer32; cost,x : Standard_Floating_Vectors.Vector )
                return double_float is

    res : double_float := 0.0;

  begin
    for i in 1..n loop
      res := res + cost(i)*x(i);
    end loop;
    return res;
  end Eval;

  procedure Eval ( n,m : in integer32;
                   cff : in Standard_Floating_Matrices.Matrix;
                   rhs,cost,x : in Standard_Floating_Vectors.Vector;
                   res : out Standard_Floating_Vectors.Vector;
                   val : out double_float ) is
  begin
    for j in 1..m loop
      res(j) := Eval(n,j,cff,x) - rhs(j);
    end loop;
    val := Eval(n,cost,x);
  end Eval;
             
-- SET UP ROUTINES :

  procedure Feasibility_Model
              ( n,m : in integer32;
                cff : in Standard_Floating_Matrices.Matrix;
                rhs : in Standard_Floating_Vectors.Vector;
                cols : in out Standard_Integer_Vectors.Vector;
                init : in out Standard_Floating_Vectors.Vector;
                extcff : out Standard_Floating_Matrices.Matrix;
                extrhs,cost : out Standard_Floating_Vectors.Vector ) is

    prd : double_float;

  begin
    for j in 1..m loop                  -- find a value for the slack var s
      extrhs(j) := rhs(j);
      for i in 1..n loop
        extcff(i,j) := cff(i,j);
      end loop;
      if Is_In(n,cols,j) then
        extcff(n+1,j) := 0.0;       -- active row, no slack var needed
      else
        extcff(n+1,j) := 1.0;       -- passive row, check residual
        prd := 0.0;
        for i in 1..n loop
          prd := prd + init(i)*cff(i,j);
        end loop;
        prd := rhs(j) - prd;
        if prd > init(n+1)
         then init(n+1) := prd; cols(n+1) := j;
        end if;
      end if;
    end loop;
    for i in 1..n loop                  -- add s >= 0 and define cost vector
      extcff(i,m+1) := 0.0;
      cost(i) := 0.0;
    end loop;
    extrhs(m+1) := 0.0;
    extcff(n+1,m+1) := 1.0;
    cost(n+1) := 1.0; 
  end Feasibility_Model;

  procedure Initialize 
              ( n,m : in integer32;
                cff : in Standard_Floating_Matrices.Matrix;
                rhs,x : in Standard_Floating_Vectors.Vector;
                active : out Standard_Integer_Vectors.Vector;
                extcff,basinv : out Standard_Floating_Matrices.Matrix;
                extrhs : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in 1..n loop
      active(i) := i;
      extrhs(i) := -x(i);
      for j in 1..n loop
        if i = j
         then basinv(i,j) := -1.0; extcff(i,j) := -1.0;
         else basinv(i,j) := 0.0; extcff(i,j) := 0.0;
        end if;
      end loop;
    end loop;
    for j in 1..m loop
      extrhs(n+j) := rhs(j);
      for i in 1..n loop
        extcff(i,n+j) := cff(i,j);
      end loop;
    end loop;
  end Initialize;

-- FEASIBILITY TEST :

  procedure Test_Basis ( bas,cff : in Standard_Floating_Matrices.Matrix;
                         cols : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Makes the product of the basis with the columns in
  --   the coefficient matrix.

    prd : double_float;

  begin
    for i in bas'range(1) loop
      for j in bas'range(2) loop
        prd := 0.0;
        for k in cff'range(1) loop
          prd := prd + bas(i,k)*cff(k,cols(j));
        end loop;
       -- put(prd,3);
      end loop;
     -- new_line;
    end loop;
  end Test_Basis;

  procedure Feasible ( n,m : in integer32; 
                       cff : in Standard_Floating_Matrices.Matrix;
                       rhs : in Standard_Floating_Vectors.Vector;
                       tol : in double_float;
                       sol : out Standard_Floating_Vectors.Vector;
                       binv : out Standard_Floating_Matrices.Matrix;
                       active : out Standard_Integer_Vectors.Vector;
                       fail : out boolean ) is

    wrkbinv : Standard_Floating_Matrices.Matrix(1..n,1..n);
    cols : Standard_Integer_Vectors.Vector(1..n);
    fail_rank : boolean;                -- true if rank(cff) < n
    extcff : Standard_Floating_Matrices.Matrix(1..n+1,1..m+1);
    extrhs : Standard_Floating_Vectors.Vector(1..m+1);
    extbas : Standard_Floating_Matrices.Matrix(1..n+1,1..n+1);
    extcols : Standard_Integer_Vectors.Vector(1..n+1);
    cost,init : Standard_Floating_Vectors.Vector(1..n+1);
    infty : boolean;
    shift : integer32;

  begin
    Initial_Basis(n,m,cff,tol,wrkbinv,cols,fail_rank);
    if fail_rank
     then
       fail := true;  -- not necessarily so !!!
      -- put_line("Coefficient matrix is not of full rank.");
     else
       init(1..n) := Solve(n,wrkbinv,rhs,cols);
       init(n+1) := 0.0;
       extcols(1..n) := cols;
       extcols(n+1) := m+1;
       Feasibility_Model(n,m,cff,rhs,extcols,init,extcff,extrhs,cost);
       for i in 1..n loop                    -- extend the basis
         for j in 1..n loop
           extbas(i,j) := wrkbinv(i,j);
         end loop;
         extbas(i,n+1) := 0.0;
         for j in 1..n loop
           extbas(i,n+1) := extbas(i,n+1)
                          - wrkbinv(i,j)*extcff(j,extcols(n+1));
         end loop;
         extbas(n+1,i) := 0.0;
       end loop;
       extbas(n+1,n+1) := 1.0;
       if init(n+1) > tol
        then Minimize(n+1,m+1,extcff,cost,extrhs,init,extbas,extcols,tol,infty);
       end if;
       if init(n+1) < tol then
         fail := false;
         shift := 0;
         for i in 1..n loop
           if extcols(i) = m+1
            then shift := 1;
           end if;
           active(i) := extcols(i+shift);
           for j in 1..n loop
             binv(i,j) := extbas(i+shift,j);
           end loop;
         end loop;
         sol(1..n) := init(1..n);
       else
         fail := true;
       end if;
    end if;
  end Feasible;

  procedure Feasible ( n,m : in integer32; 
                       cff : in Standard_Floating_Matrices.Matrix;
                       rhs : in Standard_Floating_Vectors.Vector;
                       tol : in double_float;
                       sol : out Standard_Floating_Vectors.Vector;
                       fail : out boolean ) is

    binv : Standard_Floating_Matrices.Matrix(1..n,1..n);
    cols : Standard_Integer_Vectors.Vector(1..n);
    fail_rank : boolean;                -- true if rank(cff) < n
    extcff : Standard_Floating_Matrices.Matrix(1..n+1,1..m+1);
    extrhs : Standard_Floating_Vectors.Vector(1..m+1);
    extbas : Standard_Floating_Matrices.Matrix(1..n+1,1..n+1);
    extcols : Standard_Integer_Vectors.Vector(1..n+1);
    cost,init : Standard_Floating_Vectors.Vector(1..n+1);
    infty : boolean;

  begin
    Initial_Basis(n,m,cff,tol,binv,cols,fail_rank);
    if fail_rank
     then
       fail := true;  -- not necessarily so !!!
     else
       init(1..n) := Solve(n,binv,rhs,cols);
       init(n+1) := 0.0;
       extcols(1..n) := cols;
       extcols(n+1) := m+1;
       Feasibility_Model(n,m,cff,rhs,extcols,init,extcff,extrhs,cost);
       for i in 1..n loop                    -- extend the basis
         for j in 1..n loop
           extbas(i,j) := binv(i,j);
         end loop;
         extbas(i,n+1) := 0.0;
         for j in 1..n loop
           extbas(i,n+1) := extbas(i,n+1) - binv(i,j)*extcff(j,extcols(n+1));
         end loop;
         extbas(n+1,i) := 0.0;
       end loop;
       extbas(n+1,n+1) := 1.0;
       if init(n+1) > tol then
         Minimize(n+1,m+1,extcff,cost,extrhs,init,extbas,extcols,tol,infty);
       end if;
       if init(n+1) < tol
        then fail := false; sol(1..n) := init(1..n);
        else fail := true;
       end if;
    end if;
  end Feasible;

-- MINIMIZERS :

  procedure Minimize ( n,m : in integer32;
                       cff : in Standard_Floating_Matrices.Matrix;
                       cost,rhs : in Standard_Floating_Vectors.Vector;
                       x : in out Standard_Floating_Vectors.Vector;
                       tol : in double_float; infty : out boolean ) is

    active : Standard_Integer_Vectors.Vector(1..n);
    extcff : Standard_Floating_Matrices.Matrix(1..n,1..m+n);
    extrhs : Standard_Floating_Vectors.Vector(1..m+n);
    binv : Standard_Floating_Matrices.Matrix(1..n,1..n);

  begin
    Initialize(n,m,cff,rhs,x,active,extcff,binv,extrhs);
    Minimize(n,m+n,extcff,cost,extrhs,x,binv,active,tol,infty);
  end Minimize;
             
  procedure Minimize ( n,m : in integer32;
                       cff : in Standard_Floating_Matrices.Matrix;
                       cost,rhs : in Standard_Floating_Vectors.Vector;
                       x : in out Standard_Floating_Vectors.Vector;
                       binv : in out Standard_Floating_Matrices.Matrix;
                       active : in out Standard_Integer_Vectors.Vector;
                       tol : in double_float; infty : out boolean ) is

    passive : Standard_Integer_Vectors.Vector(1..m-n)
            := Complement(n,m,active);
    stop : boolean := false;
    fail : boolean := false;

  begin
    infty := false;
    loop
      Pivoting(n,m,cff,cost,rhs,binv,x,active,passive,tol,stop,infty,fail);
      exit when stop;
    end loop;
  end Minimize;

  procedure Reporting_Minimize
               ( n,m : in integer32;
                 cff : in Standard_Floating_Matrices.Matrix;
                 cost,rhs : in Standard_Floating_Vectors.Vector;
                 x : in out Standard_Floating_Vectors.Vector;
                 tol : in double_float; infty : out boolean ) is

    active : Standard_Integer_Vectors.Vector(1..n);
    passive : Standard_Integer_Vectors.Vector(1..m);  -- after extension !
    extcff : Standard_Floating_Matrices.Matrix(1..n,1..m+n);
    extrhs : Standard_Floating_Vectors.Vector(1..m+n);
    binv : Standard_Floating_Matrices.Matrix(1..n,1..n);
    stop : boolean := false;
    fail : boolean := false;
    cont : boolean := true;

  begin
    Initialize(n,m,cff,rhs,x,active,extcff,binv,extrhs);
    passive := Complement(n,m+n,active);
    infty := false;
    loop
      Pivoting(n,m,extcff,cost,extrhs,binv,x,active,passive,
               tol,stop,infty,fail);
      Report(x,active,cont);
      exit when stop or not cont;
    end loop;
  end Reporting_Minimize;

-- ENUMERATOR :

  function Valid_Solution ( n : integer32;
                            cff : Standard_Floating_Matrices.Matrix;
                            rhs,x : Standard_Floating_Vectors.Vector;
                            tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the vector x satisfies all inequalities,
  --   returns false if one inequality is not satisfied.

  begin
    for i in rhs'range loop
      if Eval(n,i,cff,x) < rhs(i)-tol
       then return false;
      end if;
    end loop;
    return true;
  end Valid_Solution;

  function Minimum ( v : Standard_Integer_Vectors.Vector;
                     k,b : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if v(i) < b for all i in v'range, except i = k.
 
  begin
    for i in v'range loop
      if i /= k then
        if v(i) >= b
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Minimum;

  procedure Enumerate_Feasible_Vertices
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  cost,rhs : in Standard_Floating_Vectors.Vector;
                  tol : in double_float; fail,infty : out boolean ) is

    active,previous : Standard_Integer_Vectors.Vector(1..n);
    passive : Standard_Integer_Vectors.Vector(1..m-n);
    binv : Standard_Floating_Matrices.Matrix(1..n,1..n);
    sol : Standard_Floating_Vectors.Vector(1..n);
    cont : boolean := true;

    use Standard_Floating_Vectors;

    procedure Reverse_Pivoting
                ( bas : in Standard_Floating_Matrices.Matrix;
                  act,pas,prev : in Standard_Integer_Vectors.Vector ) is

    -- DESCRIPTION :
    --   This procedure enumerates all possible combinations of pivots
    --   for entering/leaving the active set of constraints.
    --   A valid choice does not return back to prev, leads to a nonsingular
    --   basis, and when yields act again when pivoting is applied.
    --   For every valid choice, Report is invoked, and if cont is true,
    --   the reverse search continues.

      leave_ind,enter_ind,ind,new_leave_var,new_enter_var : integer32;
      newbas : Standard_Floating_Matrices.Matrix(bas'range(1),bas'range(2));
      newact : Standard_Integer_Vectors.Vector(act'range);
      newpas : Standard_Integer_Vectors.Vector(pas'range);
      newsol : Standard_Floating_Vectors.Vector(act'range);
      s : double_float;
      sing,degen : boolean;

    begin
      for i in active'range loop               -- enumerate leaving variables
        leave_ind := i;
        for j in passive'range loop            -- enumerate entering variables
          enter_ind := j;
          newact := act; newpas := pas;
          ind := newpas(enter_ind);
          newpas(enter_ind) := newact(leave_ind);
          newact(leave_ind) := ind;
          if not Standard_Integer_Vectors.Equal(newact,prev) then
            newbas := bas;
            Update(n,m,newbas,cff,newact,leave_ind,ind,tol,sing);
            if not sing then
              new_leave_var := Leave_Variable(cost,newbas,tol);
              if new_leave_var > 0 then
                newsol := Solve(n,newbas,rhs,newact);
                if Valid_Solution(n,cff,rhs,newsol,tol) then
                  Enter_Variable
                    (cff,newbas,rhs,newsol,newpas,new_leave_var,tol,
                     new_enter_var,s,degen);
                  if new_enter_var > 0 and new_enter_var < newpas'last+1 then
                    if not (degen and then Minimum(act,leave_ind,ind)) then
                      if ((act(leave_ind) = newpas(new_enter_var))
                          and then (pas(enter_ind) = newact(new_leave_var)))
                       then Report(newbas,newact,newsol,cont);
                            if cont
                             then Reverse_Pivoting(newbas,newact,newpas,act);
                            end if;
                      end if;
                    end if;
                  end if;
                end if;
              end if;
            end if;
          end if;
        end loop;
      end loop;
    end Reverse_Pivoting;

  begin
    Feasible(n,m,cff,rhs,tol,sol,binv,active,fail);
    if not fail then
      Minimize(n,m,cff,cost,rhs,sol,binv,active,tol,infty);
      Report(binv,active,sol,cont);
      passive := Complement(n,m,active);
      previous := (1..n => 0);
      Reverse_Pivoting(binv,active,passive,previous);
    end if;
  end Enumerate_Feasible_Vertices;

end Linear_Minimization;
