with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Brackets;                          use Brackets;
with Numeric_Minor_Equations;           use Numeric_Minor_Equations; 
with Determinantal_Systems;             use Determinantal_Systems;
with Specialization_of_Planes;          use Specialization_of_Planes;
with Curves_into_Grassmannian;          use Curves_into_Grassmannian;

package body Pieri_Homotopies is

-- AUXILIARIES TO THE QUANTUM CASE :

  function First_Bottom_Pivot
            ( mat : Standard_Complex_Matrices.Matrix ) return integer32 is

  -- DESCRIPTION :
  --   Returns the row index to the bottommost nonzero entry of first
  --   column of the matrix m.  If the whole first column is zero, then
  --   the index to the first row minus one is returned.

    tol : constant double_float := 10.0**(-12);

  begin
    for i in reverse mat'range(1) loop
      if AbsVal(mat(i,1)) > tol
       then return i;
      end if;
    end loop;
    return mat'first(1)-1;
  end First_Bottom_Pivot;

  procedure Multiply ( p : in out Poly;
                       var : in integer32; deg : natural32 ) is

  -- DESCRIPTION :
  --   Multiplies p with x(var)**deg.

    procedure Multiply_Term ( t : in out Term; continue : out boolean ) is
    begin
      t.dg(var) := t.dg(var) + deg;
      continue := true;
    end Multiply_Term;
    procedure Multiply_Terms is new Changing_Iterator(Multiply_Term);

  begin
    Multiply_Terms(p);
  end Multiply;

  procedure Multiply ( xpm : in out Standard_Complex_Poly_Matrices.Matrix;
                       row,col,var : in integer32; deg : in natural32 ) is

  -- DESCRIPTION :
  --   Multiplies all polynomials in the nonzero rows starting with row
  --   in the column col of xpm with x(var)**deg.

  begin
    for i in row..xpm'last(1) loop
      if xpm(i,col) /= Null_Poly
       then Multiply(xpm(i,col),var,deg);
      end if;
    end loop;
  end Multiply;

  procedure Add_First ( xpm : in out Standard_Complex_Poly_Matrices.Matrix;
                        col : in integer32 ) is

  -- DESCRIPTION :
  --   Adds the first column to the indicated column.

  begin
    for i in xpm'range loop
      if xpm(i,1) /= Null_Poly then
        if xpm(i,col) = Null_Poly
         then Copy(xpm(i,1),xpm(i,col));
         else Add(xpm(i,col),xpm(i,col));
        end if;
      end if;
    end loop;
  end Add_First;

--  function Moving_Parameter0 ( n,xk,tk : natural;
--                               start,target : Complex_Number ) return Poly is
--
--  -- DESCRIPTION :
--  --   Returns the equation (x-start)*(1-t) + (x-target)*t = 0 that describes
--  --   the motion of x from start to target as t goes from 0 to 1.
--  --   Note that this equation equals x + (start-target)*t - start = 0.
--  --   This is the older version without using a random constant.
--
--  -- ON ENTRY :
--  --   n         total number of variables, continuation parameter t included;
--  --   xk        index of the moving variable x;
--  --   tk        index of the continuation parameter
--  --   start     starting value for x;
--  --   target    target value for x.
--
--    res : Poly;
--    t : Term;
--
--  begin
--    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
--    t.cf := -start;
--    res := Create(t);            -- res = -start
--    t.cf := Create(1.0);
--    t.dg(xk) := 1;
--    Add(res,t);                  -- res = x - start
--    t.dg(xk) := 0; 
--    t.dg(tk) := 1;
--    t.cf := start - target;      -- res = x + (start - target)*t - start
--    Add(res,t);
--    Clear(t);
--    return res;
--  end Moving_Parameter0;

  function Moving_Parameter ( n,xk,tk : integer32;
                              start,target : Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns the equation c*(x-start)*(1-t) + (x-target)*t = 0 for
  --   the motion of x from start to target as t goes from 0 to 1.
  --   This version uses a constant c, randomly generated from within.

  -- ON ENTRY :
  --   n         total number of variables, continuation parameter t included;
  --   xk        index of the moving variable x;
  --   tk        index of the continuation parameter t;
  --   start     starting value for x;
  --   target    target value for x.

    res : Poly;
    t : Term;
    c : constant Complex_Number := Random1;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    t.cf := -c*start;
    res := Create(t);          -- res = -c*start
    t.cf := -t.cf;
    t.dg(tk) := 1;
    Add(res,t);                -- res = -c*start + c*start*t = -c*start*(1-t)
    t.cf := c;
    t.dg(tk) := 0;
    t.dg(xk) := 1;
    Add(res,t);                -- res = -c*start*(1-t) + c*x
    t.cf := -t.cf;
    t.dg(tk) := 1;           
    Add(res,t);                -- res = -c*start*(1-t) + c*x - c*x*t
    t.cf := Create(1.0);       --     = -c*start*(1-t) + c*x*(1-t)
    Add(res,t);                --     = c*(1-t)*(x-start)
    t.cf := -target;           -- res = c*(1-t)*(x-start) + x*t
    t.dg(xk) := 0;
    Add(res,t);                -- res = c*(1-t)*(x-start) + x*t - target*t
    Clear(t);                  --     = c*(1-t)*(x-start) + (x-target)*t
    return res;
  end Moving_Parameter;

  function Constant_Parameter ( n,i : integer32; s : Complex_Number )
                              return Poly is

  -- DESCRIPTION :
  --   Returns the polynomial x_i - s; n equals the number of variables.

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    t.dg(i) := 1;
    t.cf := Create(1.0);
    res := Create(t);
    t.dg(i) := 0;
    t.cf := -s;
    Add(res,t);
    Clear(t);
    return res;
  end Constant_Parameter;

-- TARGET ROUTINES :

  function One_Hypersurface_Pieri_Homotopy
                ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                  xpm : Standard_Complex_Poly_Matrices.Matrix;
                  planes : VecMat ) return Poly_Sys is

    res : Poly_Sys(1..nd.level);
    p : constant integer32 := nd.p;
    m : constant integer32 := n-p;
    special : Standard_Complex_Matrices.Matrix(1..n,1..m);
    target,start : Poly;

  begin
    case nd.tp is
      when top    => special := Special_Plane(m,nd.top);
                    -- special := Special_Top_Plane(m,nd.top);
      when bottom => special := Special_Plane(m,nd.bottom);
                    -- special := Special_Bottom_Plane(m,nd.bottom);
      when others => null;  -- mixed case treated separately 
    end case;
    for i in 1..nd.level-1 loop
      res(i) := Expanded_Minors(planes(i).all,xpm,expbp);
      Embed(res(i));
    end loop;
    target := Expanded_Minors(planes(nd.level).all,xpm,expbp);
    start := Expanded_Minors(special,xpm,expbp);
    res(nd.level) := Linear_Homotopy(target,start);
    Clear(target); Clear(start);
    return res;
  end One_Hypersurface_Pieri_Homotopy;

  function Two_Hypersurface_Pieri_Homotopy
                ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                  xpm : Standard_Complex_Poly_Matrices.Matrix;
                  planes : VecMat ) return Poly_Sys is

    res : Poly_Sys(1..nd.level);
    p : constant integer32 := nd.p;
    m : constant integer32 := n-p;
    top_special : constant Standard_Complex_Matrices.Matrix
                := Special_Plane(m,nd.top);
               -- := Special_Top_Plane(m,nd.top);
    bot_special : constant Standard_Complex_Matrices.Matrix
                := Special_Plane(m,nd.bottom);
               -- := Special_Bottom_Plane(m,nd.bottom);
    top_start,bot_start,top_target,bot_target : Poly;

  begin
    for i in 1..nd.level-2 loop
      res(i) := Expanded_Minors(planes(i).all,xpm,expbp);
      Embed(res(i));
    end loop;
    top_target := Expanded_Minors(planes(nd.level).all,xpm,expbp);
    top_start := Expanded_Minors(top_special,xpm,expbp);
    res(nd.level) := Linear_Homotopy(top_target,top_start);
    Clear(top_target); Clear(top_start);
    bot_target := Expanded_Minors(planes(nd.level-1).all,xpm,expbp);
    bot_start := Expanded_Minors(bot_special,xpm,expbp);
    res(nd.level-1) := Linear_Homotopy(bot_target,bot_start);
    Clear(bot_target); Clear(bot_start);
    return res;
  end Two_Hypersurface_Pieri_Homotopy;

  function One_General_Pieri_Homotopy
                ( n,ind : integer32; nd : Node; bs : Bracket_System;
                  start,target : Standard_Complex_Matrices.Matrix;
                  xpm : Standard_Complex_Poly_Matrices.Matrix;
                  planes : VecMat ) return Link_to_Poly_Sys is

    res : Link_to_Poly_Sys;
    nva : constant integer32 := n*nd.p + 1;
    moving_plane : Standard_Complex_Poly_Matrices.Matrix(1..n,target'range(2))
                 := Moving_U_Matrix(nva,start,target);
    moving : Poly_Sys(1..bs'last);

  begin
    for i in 1..ind-1 loop
      Concat(res,Polynomial_Equations(planes(i).all,xpm));
    end loop;
    if res /= null then
      for i in res'range loop
        Embed(res(i));
      end loop;
    end if;
    moving := Lifted_Expanded_Minors(moving_plane,xpm,bs); 
    Concat(res,moving);
    Standard_Complex_Poly_Matrices.Clear(moving_plane);
    return res;
  end One_General_Pieri_Homotopy;

  function Two_General_Pieri_Homotopy
              ( n,ind : integer32; nd : Node; top_bs,bot_bs : Bracket_System;
                top_start,top_target,bot_start,bot_target
                  : Standard_Complex_Matrices.Matrix;
                xpm : Standard_Complex_Poly_Matrices.Matrix;
                planes : VecMat ) return Link_to_Poly_Sys is

  -- DESCRIPTION :
  --   Returns the homotopy for general linear subspace intersections,
  --   in case nd.tp = mixed.  The parameter ind indicates the plane
  --   planes(ind) towards the constructed homotopy works.

    res : Link_to_Poly_Sys;
    nva : constant integer32 := n*nd.p + 1;
    top_moving : Standard_Complex_Poly_Matrices.Matrix(1..n,top_target'range(2))
               := Moving_U_Matrix(nva,top_start,top_target);
    bot_moving : Standard_Complex_Poly_Matrices.Matrix(1..n,bot_target'range(2))
               := Moving_U_Matrix(nva,bot_start,bot_target);
    top_movsys : Poly_Sys(1..top_bs'last);
    bot_movsys : Poly_Sys(1..bot_bs'last);

  begin
    for i in 1..ind-1 loop
      Concat(res,Polynomial_Equations(planes(i).all,xpm));
    end loop;
    if res /= null
     then for i in res'range loop
            Embed(res(i));
          end loop;
    end if;
    top_movsys := Lifted_Expanded_Minors(top_moving,xpm,top_bs); 
    bot_movsys := Lifted_Expanded_Minors(bot_moving,xpm,bot_bs); 
    Concat(res,top_movsys);
    Concat(res,bot_movsys);
    Standard_Complex_Poly_Matrices.Clear(top_moving);
    Standard_Complex_Poly_Matrices.Clear(bot_moving);
    return res;
  end Two_General_Pieri_Homotopy;

  function One_Quantum_Pieri_Homotopy
                ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                  xpm : Standard_Complex_Poly_Matrices.Matrix;
                  planes : VecMat; s : Vector ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the Pieri homotopy that corresponds to the node.
  --   This homotopy is set up to work only in the hypersurface case,
  --   when the type of the node is either top or bottom.

    res : Poly_Sys(1..nd.level+1);
    p : constant integer32 := nd.p;
    m : constant integer32 := n-p;
    nvars : constant integer32 := nd.level + p + 2;
      -- nvars = level   #equations in the x_ij's
      --       + p       because not yet fixed the ones
      --       + 2       for s and t, note that t is continuation parameter
    special : Standard_Complex_Matrices.Matrix(1..n,1..m);
    target,start : Poly;
    map : Standard_Complex_Poly_Matrices.Matrix(1..n,1..p);

  begin
    case nd.tp is
      when top    => special := Special_Plane(m,Modulo(nd.top,natural32(m+p)));
                     Standard_Complex_Poly_Matrices.Copy(xpm,map);
                     Swap(map,nvars-1,nvars);
      when bottom => special := Special_Plane
                                  (m,Modulo(nd.bottom,natural32(m+p)));
                     map := xpm;
      when others => null;  -- mixed case treated separately
    end case;
    for i in 1..nd.level-1 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
            := Eval(xpm,s(i),Create(1.0));
      begin
        res(i) := Expanded_Minors(planes(i).all,eva,expbp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    target := Expanded_Minors(planes(nd.level).all,map,expbp);
    start := Expanded_Minors(special,map,expbp);
    res(nd.level) := Linear_Interpolation(target,start,nvars);
    if nd.tp = bottom
     then
       res(nd.level+1)
          := Moving_Parameter(nvars,nvars-1,nvars,Create(1.0),s(nd.level));
     else
       res(nd.level+1)
          := Moving_Parameter(nvars,nvars-1,nvars,Create(1.0),
                              (Create(1.0)/s(nd.level)));
       Divide_Common_Factor(res(nd.level),nvars);
    end if;
    if nd.tp = top
     then Standard_Complex_Poly_Matrices.Clear(map);
    end if;
    Clear(target); Clear(start);
    return res;
  end One_Quantum_Pieri_Homotopy;

  function Two_Quantum_Pieri_Homotopy
                ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                  xpm : Standard_Complex_Poly_Matrices.Matrix;
                  planes : VecMat; s : Vector ) return Poly_Sys is

    res : Poly_Sys(1..nd.level+2);
    p : constant integer32 := nd.p;
    m : constant integer32 := n-p;
    nvars : constant integer32 := nd.level + p + 2;
    top_special : constant Standard_Complex_Matrices.Matrix
                := Special_Plane(m,Modulo(nd.top,natural32(m+p)));
    bot_special : constant Standard_Complex_Matrices.Matrix
                := Special_Plane(m,Modulo(nd.bottom,natural32(m+p)));
    top_start,bot_start,top_target,bot_target : Poly;
    map : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
        := Insert(xpm,nvars);

  begin
   -- first do bottom pivots with s1 = nvars-1
    bot_target := Expanded_Minors(planes(nd.level-1).all,map,expbp);
    bot_start := Expanded_Minors(bot_special,map,expbp);
    res(nd.level-1) := Linear_Interpolation(bot_target,bot_start,nvars+1);
    Divide_Common_Factor(res(nd.level-1),nvars+1);
    res(nd.level+1)
      := Moving_Parameter(nvars+1,nvars-1,nvars+1,Create(1.0),s(nd.level-1));
    Clear(bot_target); Clear(bot_start);
   -- swap s1 with s2 to deal with the fixed equations
    Swap(map,nvars-1,nvars);
    for i in 1..nd.level-2 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(map'range(1),map'range(2))
            := Eval(map,s(i),Create(1.0));
      begin
        res(i) := Expanded_Minors(planes(i).all,eva,expbp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
   -- swap t with s2 for top pivots, s2 = nvars
    Swap(map,nvars,nvars+1);
    top_target := Expanded_Minors(planes(nd.level).all,map,expbp);
    top_start := Expanded_Minors(top_special,map,expbp);
    res(nd.level) := Linear_Interpolation(top_target,top_start,nvars+1);
    res(nd.level+2)
      := Moving_Parameter(nvars+1,nvars,nvars+1,Create(1.0),
                          (Create(1.0)/s(nd.level)));
    Divide_Common_Factor(res(nd.level),nvars+1);
    Clear(top_target); Clear(top_start);
    Standard_Complex_Poly_Matrices.Clear(map);
    return res;
  end Two_Quantum_Pieri_Homotopy;

  function One_General_Quantum_Pieri_Homotopy
                  ( n,ind : integer32; nd : Node; s_mode : natural32;
                    bs : Bracket_System;
                    start,target : Standard_Complex_Matrices.Matrix;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat; s : Vector ) return Link_to_Poly_Sys is

    res : Link_to_Poly_Sys;
    p : constant integer32 := nd.p;
    nvars : constant integer32 := nd.level + p + 2;
      -- nvars = level   #equations in the x_ij's
      --       + p       because not yet fixed the ones
      --       + 2       for s and t, note that t is continuation parameter
    eva,map : Standard_Complex_Poly_Matrices.Matrix(1..n,1..p);
    moving_plane
      : constant Standard_Complex_Poly_Matrices.Matrix(1..n,target'range(2))
      := Moving_U_Matrix(nvars,start,target);
    moving : Poly_Sys(1..bs'last);
    movpar : Poly_Sys(1..1);
    inddiv,row : integer32;

  begin
   -- PART I : fixed equations
    case nd.tp is
      when top    => Standard_Complex_Poly_Matrices.Copy(xpm,map);
                     Swap(map,nvars-1,nvars);
      when bottom => map := xpm;
      when others => null;  -- mixed case treated separately
    end case;
    if s_mode = 0
     then row := First_Bottom_Pivot(start);
          if row > 1                            -- multiply 2nd column with t
           then --Multiply(map,row+1,2,nvars,1);  -- from row+1 on
                Multiply(map,row,2,nvars-1,1);
                Add_First(map,2);               -- add first column to second
          end if;
    end if;
    for i in 1..ind-1 loop
      eva := Eval(map,s(i),Create(1.0));
      Concat(res,Polynomial_Equations(planes(i).all,eva));
      Standard_Complex_Poly_Matrices.Clear(eva);
    end loop;
    if res = null
     then inddiv := 0;
     else inddiv := res'last;
    end if;
   -- PART II : moving equation for the plane
    --if s_mode = 1 or s_mode = 2
    -- then
    moving := Expanded_Minors(moving_plane,map,bs);
    -- else eva := Eval(map,Create(1.0),Create(0.0));
         -- if s_mode = 0
         --  then row := First_Bottom_Pivot(start);
         --       if row > 1                              -- multiply 2nd column
         --        then Multiply(eva,row+1,2,nvars-1,1);  -- with s, from row+1 
         --       end if;
         -- end if;
         -- moving := Expanded_Minors(moving_plane,eva,bs);
         -- Standard_Complex_Poly_Matrices.Clear(eva);
   -- end if;
    Concat(res,moving);
    if nd.tp = top
     then Standard_Complex_Poly_Matrices.Clear(map);
          for i in inddiv+1..res'last loop
            Divide_Common_Factor(res(i),nvars);
          end loop;
    end if;
   -- PART III : moving equation for the interpolation points
    case s_mode is
      when 0 =>                                    -- move s from 0 to 1
        case nd.tp is
          when bottom =>
            movpar(1)
              := Moving_Parameter(nvars,nvars-1,nvars,
                                  Create(0.0),Create(1.0));
          when others => null;
        end case;
      when 1 =>                                    -- leave s constant at 1
        movpar(1) := Constant_Parameter(nvars,nvars-1,Create(1.0));
      when 2 =>                                    -- move s from 1 to target
        case nd.tp is
          when top =>
            movpar(1)     -- move s from 1 to 1/s_ind (swap of s and t)
              := Moving_Parameter(nvars,nvars-1,nvars,Create(1.0),
                                  (Create(1.0)/s(ind)));
          when bottom =>
            movpar(1)     -- move s from 1 to s_ind
              := Moving_Parameter(nvars,nvars-1,nvars,Create(1.0),s(ind));
          when others => null;
        end case;
      when others => null;
    end case;
    Concat(res,movpar);
    return res;
  end One_General_Quantum_Pieri_Homotopy;

 -- function Two_General_Quantum_Pieri_Homotopy
 --              ( n,ind : natural; nd : Node; top_bs,bot_bs : Bracket_System;
 --                top_start,top_target,bot_start,bot_target
 --                  : Standard_Complex_Matrices.Matrix;
 --                xpm : Standard_Complex_Poly_Matrices.Matrix;
 --                planes : VecMat; s : Vector ) return Link_to_Poly_Sys is
 --
 --    res : Link_to_Poly_Sys;
 --
 --  begin
 --    return res;
 --  end Two_General_Quantum_Pieri_Homotopy;

end Pieri_Homotopies;
