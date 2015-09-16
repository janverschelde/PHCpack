with unchecked_deallocation;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors;
with QuadDobl_Random_Numbers;           use QuadDobl_Random_Numbers;
with QuadDobl_Complex_VecVecs;          use QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Sample_Points;            use QuadDobl_Sample_Points;
with QuadDobl_Rectangular_Sample_Grids; use QuadDobl_Rectangular_Sample_Grids;

package body QuadDobl_Divided_Differences is

-- DATA STRUCTURES :

  type Newton_Interpolator1_Rep ( n,d,d1 : integer32 ) is record
    rot : Vector(1..n);       -- normal to all parallel hyperplanes in grid
    pts : QuadDobl_Complex_Matrices.Matrix(0..d,0..d1);      -- data points
  end record;

 -- DATA REPRESENTATION CONVENTIONS :
 --   n = original #variables, d = degree, d1 = d+1;
 --   rot contains the normal to all parallel hyperplanes in the grid.
 --   pts(i,0) = x(i), pts(i,j) is y-value with x(i), for i in 0..d,
 --   pts(0,d+1) = c, pts(i,d+1) = p(x(i),c), for i in 1..d.
 --   So the interpolating polynomial p satisfies
 --       p(pts(i,0),pts(i,j)) = 0, for i in 0..d, j in 1..d;
 --       p(pts(0,0),pts(0,d+1)) = 1;
 --   and p(pts(i,0),pts(0,d+1)) = pts(i,d+1), for i in 1..d.

  type Newton_Form_Evaluator1_Rep ( n,d : integer32 ) is record
    rot : Vector(1..n);   -- normal to hyperplane, encodes rotation
    pts : Vector(0..d);   -- first coordinates of interpolation points
    dvd : VecVec(0..d);   -- coefficients of polynomials in 2nd variable
  end record;

  type Newton_Taylor_Form1_Rep ( n,d : integer32 ) is record
    rot : Vector(1..n);   -- normal to hyperplane, encodes rotation
    x,y : Vector(0..d);   -- x and y coordinates
                          -- dvd(i,j) = f[x0..xi;y0..yj]
    dvd : QuadDobl_Complex_Matrices.Matrix(0..d,0..d); 
  end record;

  type Newton_Taylor_Forms is
    array ( integer32 range <> ) of Newton_Taylor_Form;

  type Newton_Taylor_Form_Rep ( k,n,d : integer32 ) is record
    case k is   -- k = #slices; n = #variables; d = degree
      when 1      => q1 : Newton_Taylor_Form1_Rep(n,d);
      when others => rot : Vector(1..n);   -- k-th slice
                     pts : Vector(0..d);   -- evaluations at k-th slice
                     qp : Newton_Taylor_Forms(0..d);  -- deg(qp(i)) = i
    end case;
  end record;

-- NUMERICAL EVALUATION OF THE TABLE OF DIVIDED DIFFERENCES :

  procedure Numerical_Divided_Differences
               ( dvd : in out Vector; x2 : in Complex_Number;
                 pts : in QuadDobl_Complex_Matrices.Matrix;
                 i,j,d : in integer32 ) is

  -- DESCRIPTION :
  --   Updates the vector of divided differences for a curve.

  -- ON ENTRY :
  --   dvd       current last row in the table of divided differences;
  --   x2        second coordinate of the data point;
  --   pts       compact matrix form of the interpolation grid;
  --   i         index to current row in table of divided differences;
  --   j         index to current column in table and dvd;
  --   d         degree and size of the table.

  -- ON RETURN :
  --   dvd       dvd(j) will contained updated divided differences.

    acc : Complex_Number;

  begin
    if i = j then
      if i = 0
       then dvd(j) := Create(integer(1));
       else dvd(j) := pts(i,d+1);
      end if;
      for k in 1..d loop
        acc := (x2 - pts(i,k))/(pts(0,d+1) - pts(i,k));
        dvd(j) := dvd(j)*acc;
      end loop;
    else
      dvd(j) := (dvd(j) - dvd(j+1))/(pts(j,0) - pts(i,0));
    end if;
  end Numerical_Divided_Differences;

  function Eval ( d : integer32; pts : QuadDobl_Complex_Matrices.Matrix;
                  x : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   This routine returns the result of the Newton form at the points.

    res : Complex_Number;
    dvd : Vector(0..d);

  begin
    for i in 0..dvd'last loop
      for j in reverse 0..i loop
        Numerical_Divided_Differences(dvd,x(2),pts,i,j,d);
      end loop;
    end loop;
    res := dvd(0);
    for i in 1..d loop
      res := res*(x(1) - pts(i,0));
      res := res + dvd(i);
    end loop;
    return res;
  end Eval;

  function Eval ( q : Newton_Interpolator1_Rep; x : Vector )
                return Complex_Number is

    rpx : constant Vector(1..2) := Rotate_and_Project(q.rot,x);
    res : constant Complex_Number := Eval(q.d,q.pts,rpx);

  begin
    return res;
  end Eval;

  function Eval ( q : Newton_Taylor_Form1_Rep; x : Vector )
                return Complex_Number is

    res : Complex_Number := q.dvd(q.d,0);
    acc : Complex_Number;
    rpx : constant Vector(1..2) := Rotate_and_Project(q.rot,x);

  begin
    for i in reverse 0..(q.d-1) loop        -- reverse row wise in q.dvd
      acc := q.dvd(i,q.d-i);
      for j in reverse 0..(q.d-i-1) loop          -- reverse column wise
        acc := acc*(rpx(2)-q.y(j)) + q.dvd(i,j);
      end loop;
      res := res*(rpx(1)-q.x(i)) + acc;
    end loop;
    return res;
  end Eval;

  function Eval ( q : Newton_Taylor_Form_Rep; x : Vector )
                return Complex_Number is

    res,rpx : Complex_Number;

  begin
    if q.k = 1 then
      return Eval(q.q1,x);
    else
      rpx := q.rot(x'range)*x;
      res := Eval(q.qp(0),x);
      for i in 1..q.d loop
        res := res*(rpx - q.pts(q.d-i));
        res := res + Eval(q.qp(i),x);
      end loop;
      return res;
    end if;
  end Eval;

-- SYMBOLIC EVALUATION OF TABLE OF DIVIDED DIFFERENCES : 

  procedure Mul_Fac ( cf : in out Vector; y,c : in Complex_Number ) is

  -- DESCRIPTION :
  --   Multiplies the polynomial p(x) with coefficient in cf with
  --   the factor x - y and divides with the constant c - y.

  -- ON ENTRY :
  --   cf         coefficients of cf(0) + cf(1)*x + .. + cf(d)*x^d;
  --   y          particular y coordinate of interpolation point;
  --   c          random chosen complex number, different from y.

  -- ON RETURN :
  --   cf         updated coefficient vector of the polynomial.

    frac : constant Complex_Number := c - y;
    newc0,newc1 : Complex_Number;

  begin
    newc0 := -cf(0)*y;
    for i in 1..cf'last loop
      newc1 := cf(i-1) - cf(i)*y;
      cf(i-1) := newc0/frac;
      newc0 := newc1;
    end loop;
    cf(cf'last) := newc0/frac;
  end Mul_Fac;

  function Eval ( cf : Vector; x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the function value of x at the polynomial with
  --   coefficients in cf, evaluated with Horner scheme.

    res : Complex_Number := cf(cf'last);

  begin
    for i in reverse 0..cf'last-1 loop
      res := cf(i) + res*x;
    end loop;
    return res;
  end Eval;

  function Substitute ( cf : Vector; a1,a2 : Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns the substitution y := a2*x - a1*y in the polynomial with
  --   coefficient vector cf.

    res,s : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2 => 0);
    t.cf := a2;
    t.dg(1) := 1;
    s := Create(t);                         -- s = a2*x
    t.dg(1) := 0;
    t.dg(2) := 1;
    t.cf := -a1;                            -- t = -a1*y
    Add(s,t);                               -- s = a2*x - a1*y
    t.dg(2) := 0;
    t.cf := cf(cf'last);
    res := Create(t);
    for i in reverse 0..cf'last-1 loop      -- Horner scheme
      Mul(res,s);
      t.cf := cf(i);
      Add(res,t);
    end loop;
    Clear(t); Clear(s);
    return res;
  end Substitute;

  function Degree ( t : Term ) return integer32 is

  -- DESCRIPTION : returns the degree of the term.

    res : integer32 := 0;

  begin
    for i in t.dg'range loop
      res := res + integer32(t.dg(i));
    end loop;
    return res;
  end Degree;

  function Extend_and_Trunc ( p : Poly; n,d : integer32 ) return Poly is

  -- DESCRIPTION :
  --   Extends the range of the terms to n variables and truncates all
  --   terms of degree higher than d.

    res : Poly := Null_Poly;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      et : Term;

    begin
      if Degree(t) <= d then
        et.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
        et.dg(t.dg'range) := t.dg.all;
        et.cf := t.cf;
        Add(res,et);
        Clear(et);
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Extend_and_Trunc;

  function Extend_and_Trunc ( p : Poly_Sys; n : integer32 ) return Poly_Sys is

  -- DESCRIPTION :
  --   For the system res on return the following will hold: 
  --     Degree(res(i)) = i and Number_of_Unknowns(res(i)) = n.

    res : Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Extend_and_Trunc(p(i),n,i);
    end loop;
    return res;
  end Extend_and_Trunc;

  function Expand ( n : integer32; dvd : Poly_Sys; rot,pts : Vector )
                  return Poly is

  -- DESCRPTION :
  --   Returns the interpolator as a polynomial in n variables.

  -- ON ENTRY :
  --   n          number of variables in the polynomials in dvd;
  --   dvd        output of Extend_and_Trunc;
  --   rot        representation of the rotation;
  --   pts        first coordinates of the interpolation points.

    res,s,s1 : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    for i in 1..n loop
      t.dg(i) := 1;
      t.cf := rot(i);
      Add(s,t);
      t.dg(i) := 0;
    end loop;
    Copy(dvd(0),res);
    for i in 1..dvd'last loop
      t.cf := pts(i);
      s1 := s - t;
      Mul(res,s1);
      Add(res,dvd(i));
      Clear(s1);
    end loop;
    Clear(s); Clear(t);
    return res;
  end Expand;

  procedure Symbolic_Divided_Differences
               ( dvd : in out VecVec; pts : QuadDobl_Complex_Matrices.Matrix;
                 i,j,d : integer32 ) is

  -- DESCRIPTION :
  --   This procedure updates a vector of coefficient vectors of the
  --   polynomials in the second coordinate which occur in the symbolic
  --   table of divided differences.

  -- ON ENTRY :
  --   dvd       current last row in the table of divided differences;
  --   pts       compact matrix form of the interpolation grid;
  --   i         index to current row in table of divided differences;
  --   j         index to current column in table and dvd;
  --   d         degree and size of the table.

  -- ON RETURN :
  --   dvd       dvd(j) will contained updated divided differences.

    frac : Complex_Number;

  begin
    if i = j then
      dvd(j) := new Vector'(0..d => Create(integer(0)));
      if i = 0
       then dvd(j)(0) := Create(integer(1));
       else dvd(j)(0) := pts(i,d+1);
      end if;
      for k in 1..d loop
        Mul_Fac(dvd(j).all,pts(i,k),pts(0,d+1));
      end loop;
    else
      Sub(dvd(j).all,dvd(j+1).all);
      frac := Create(1.0)/(pts(j,0) - pts(i,0));
      Mul(dvd(j).all,frac);
    end if;
  end Symbolic_Divided_Differences;

-- INTERNAL TESTING ROUTINES :

  procedure Write_Divided_Differences
              ( file : in file_type; q : Newton_Taylor_Form ) is

  -- DESCRIPTION :
  --   This routine is auxiliary to test Truncate.

  begin
    if q.k = 1 then
      for i in 0..q.d loop
        for j in 0..q.d-i loop
          put(file,q.q1.dvd(i,j)); new_line(file);
        end loop;
      end loop;
    else
      for i in q.qp'range loop
        Write_Divided_Differences(file,q.qp(i));
      end loop;
    end if;
  end Write_Divided_Differences;

  procedure Test_Degrees ( file : in file_type; q : in Newton_Taylor_Form ) is
 
  -- DESCRIPTION :
  --   Writes down the actual degree information from the dvd info.

    tol : constant double_float := 10.0**(-8);
    deg : integer32;

  begin
    put(file,"Testing degrees at q.k = "); put(file,q.k,1);
    put_line(file," :");
    if q.k = 1 then
      deg := 0;
      for i in 1..q.q1.d loop
        for j in 0..i loop
          if AbsVal(q.q1.dvd(i-j,j)) > tol
           then deg := i; exit;
          end if;
          exit when (deg = (i-1));
        end loop;
      end loop;
      put(file,"Degree found : "); put(file,deg,1); new_line(file);
    else
      for i in 1..q.qp'last loop
        put(file,"Recursive call for i = "); put(file,i,1);
        put_line(file," :");
        Test_Degrees(file,q.qp(i));
      end loop;
    end if;
  end Test_Degrees;

  procedure Test_Evaluation1 
               ( file : in file_type; grid : in Stacked_Sample_Grid;
                 q : in Newton_Taylor_Form1_Rep ) is

  -- DESCRIPTION :
  --   Tests the evaluation of the 1-D Taylor form at the grid.
  --   Residuals are written to file.

  -- REQUIRED : grid.k = 1.

    tmp : QuadDobl_Sample_List;

  begin
    put_line(file,"Test evaluation of 1-D Taylor form at the grid.");
    for i in grid.g'range loop
      tmp := grid.g(i);
      while not Is_Null(tmp) loop
        put(file,Eval(q,Sample_Point(Head_Of(tmp)).v));
        new_line(file);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Test_Evaluation1;

  procedure Test_Evaluation
               ( file : in file_type; grid : in Stacked_Sample_Grid;
                 q : in Newton_Taylor_Form_Rep ) is

  -- DESCRIPTION :
  --   Tests the evaluation of the Newton-Taylor form at the grid.
  --   Residuals are written to file.

    tmp : QuadDobl_Sample_List;

  begin
    put(file,"Test evaluation of ");
    put(file,q.k,1); put(file,"-D Taylor form at ");
    put(file,grid.k,1); put_line(file,"-D grid :");
    if grid.k = 1 then
      for i in grid.g'range loop
        tmp := grid.g(i);
        while not Is_Null(tmp) loop
          put(file,Eval(q,Sample_Point(Head_Of(tmp)).v));
          new_line(file);
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    else
      for i in 1..grid.d loop
        Test_Evaluation(file,grid.a(i).all,q);
      end loop;
      put(file,"At the extra sample : ");
      put(file,Eval(q,Sample_Point(grid.spt).v));
      new_line(file);
    end if;
  end Test_Evaluation;

-- TARGET PROCEDURES :

  function Create ( grid : Array_of_QuadDobl_Sample_Lists )
                  return Newton_Interpolator1 is

    c : constant Complex_Number := Random1;

  begin
    return Create(grid,c);
  end Create;

  function Create ( grid : Array_of_QuadDobl_Sample_Lists;
                    c : Complex_Number ) return Newton_Interpolator1_Rep is

    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := grid'last;
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    res : Newton_Interpolator1_Rep(n,d,d+1);
    tmp : QuadDobl_Sample_List;
    spt : QuadDobl_Sample;
    acc : Complex_Number;
    rpx : Vector(1..2);

  begin
    res.rot := hyp(hyp'first)(1..n);
    for i in 0..d loop          -- process all samples on same slice #i
      tmp := grid(i);
      spt := Head_Of(tmp);
      rpx := Rotate_and_Project(res.rot,Sample_Point(spt).v);
      res.pts(i,0) := -Hyperplane_Sections(spt)(1)(0);   -- x coordinate
      res.pts(i,1) := rpx(2);                            -- y coordinate
      tmp := Tail_Of(tmp);
      for j in 2..d loop
        spt := Head_Of(tmp);
        res.pts(i,j) := Rotate_and_Project2(res.rot,Sample_Point(spt).v);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    res.pts(0,d+1) := c;
    for i in 1..d loop
      res.pts(i,d+1) := Create(integer(1));             -- will be p(x_i,c)
      for j in 1..d loop
        acc := (c - res.pts(i,j))/(c - res.pts(0,j));
        res.pts(i,d+1) := res.pts(i,d+1)*acc;
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( grid : Array_of_QuadDobl_Sample_Lists;
                    c : Complex_Number ) return Newton_Interpolator1 is

    res : Newton_Interpolator1;

  begin
    res := new Newton_Interpolator1_Rep'(Create(grid,c));
    return res;
  end Create;

  function Create ( q : Newton_Interpolator1 ) return Newton_Form_Evaluator1 is

    res : Newton_Form_Evaluator1;
    res_rep : Newton_Form_Evaluator1_Rep(q.n,q.d);

  begin
    res_rep.rot := q.rot;
    for i in 0..q.d loop
      res_rep.pts(i) := q.pts(i,0);
    end loop;
    for i in 0..q.d loop
      for j in reverse 0..i loop
        Symbolic_Divided_Differences(res_rep.dvd,q.pts,i,j,q.d);
      end loop;
    end loop;
    res := new Newton_Form_Evaluator1_Rep'(res_rep);
    return res;
  end Create;

  procedure Column_Divided_Differences
               ( dvd : in out QuadDobl_Complex_Matrices.Matrix;
                 x : in Vector; d,j : in integer32 ) is

  -- DESCRIPTION :
  --   Computes the j-th column in the table dvd of divided differences.

  -- ON ENTRY :
  --   dvd       j-th column contains f(x(i),y(j)), i=0,1,..,d-j;
  --   x         values for x at the data points;
  --   j         current column.

  -- ON RETURN :
  --   dvd       j-th column contains generalized divided differences:
  --               dvd(i,j) = f[x(0)..x(i);y(j)], for i=0,1,..,d-j.

  begin
    for k in 1..d-j loop
      for i in reverse k..d-j loop
        dvd(i,j) := (dvd(i,j) - dvd(i-1,j))/(x(i) - x(i-k));
      end loop;
    end loop;
  end Column_Divided_Differences;

  procedure Row_Divided_Differences
               ( dvd : in out QuadDobl_Complex_Matrices.Matrix;
                 y : in Vector; d,i : in integer32 ) is

  -- DESCRIPTION :
  --   Computes the i-th column in the table dvd of divided differences.

  -- ON ENTRY :
  --   dvd       i-th row contains f[x(0)..x(i),y(j)), j=0,1,..,d-i;
  --   y         values for y at the data points;
  --   d         degree of the interpolating polynomial;
  --   i         current row.

  -- ON RETURN :
  --   dvd       i-th row contains generalized divided differences:
  --               dvd(i,j) = f[x(0)..x(i);y(0)..y(j)], for j=0,1,..,d-i.

  begin
    for k in 1..(d-i) loop
      for j in reverse k..(d-i) loop
        dvd(i,j) := (dvd(i,j) - dvd(i,j-1))/(y(j) - y(j-k));
      end loop;
    end loop;
  end Row_Divided_Differences;

  function Create ( q : Newton_Interpolator1_Rep; y : Vector )
                  return Newton_Taylor_Form1_Rep is

    res : Newton_Taylor_Form1_Rep(q.n,q.d);
    xy : Vector(1..2);

  begin
    res.rot := q.rot;
    res.y := y;
    for i in 0..q.d loop
      res.x(i) := q.pts(i,0);
    end loop;
    for i in 0..q.d loop
      xy(1) := res.x(i);
      for j in 0..q.d-i loop
        xy(2) := y(j);
        res.dvd(i,j) := Eval(q,Inverse_Rotate(q.rot,xy));
      end loop;
    end loop;
    for j in 0..(q.d-1) loop
      Column_Divided_Differences(res.dvd,res.x,q.d,j);
    end loop;
    for i in 0..(q.d-1) loop
      Row_Divided_Differences(res.dvd,res.y,q.d,i);
    end loop;
    return res;
  end Create;

  function Create ( q : Newton_Interpolator1; y : Vector )
                  return Newton_Taylor_Form1 is

    res : Newton_Taylor_Form1;
    res_rep : Newton_Taylor_Form1_Rep(q.n,q.d);

  begin
    if q /= null
     then res_rep := Create(q.all,y);
          res := new Newton_Taylor_Form1_Rep'(res_rep);
    end if;
    return res;
  end Create;

  function Create_Constant
             ( n,deg : integer32; pts : Vector;
               spt : QuadDobl_Sample; qp : Newton_Taylor_Forms ) 
             return Newton_Taylor_Form is

  -- DESCRIPTION :
  --   Returns qp(0) in the interpolating polynomial for a k dimensional
  --   component of degree d in the ambient space of n dimensions.
  --   The points for the current variable are in pts.

  -- NOTE : since this is a constant, it will be stored as k = 1.

    res : Newton_Taylor_Form;
    res_rep : Newton_Taylor_Form_Rep(1,n,0);
    pt : constant Vector(1..n) := Sample_Point(spt).v;
    eva,dif,factor : Complex_Number;

  begin
    eva := Eval(qp(1),pt);
    factor := pts(deg) - pts(deg-1);
    for i in 2..deg loop
      dif := pts(deg) - pts(deg-i);
      factor := factor*dif;
      eva := eva*dif;
      eva := eva + Eval(qp(i),pt);
    end loop;
    res_rep.q1.dvd(0,0) := -eva/factor;
    res := new Newton_Taylor_Form_Rep'(res_rep);
    return res;
  end Create_Constant;

  procedure Adjust_Divided_Differences
               ( d : in integer32; divisor : in Complex_Number;
                 q2 : in Newton_Taylor_Form1_Rep;
                 q1 : in out Newton_Taylor_Form1_Rep ) is

  -- DESCRIPTION :
  --   Adjustment of the divided differences in q1 using q2.

  -- REQUIRED : k = 1, d > 0, and q2 is properly arranged,
  --            q1.dvd(d+1,0) /= 0.

  -- ON ENTRY :
  --   d         current degree;
  --   divisor   difference between the interpolation points;
  --   q2        interpolator of one degree up;
  --   q1        interpolator whose divided differences will be modified.

  -- ON RETURN :
  --   q1        has adjusted divided differences, up to degree <= d.

    fac : constant Complex_Number := q2.dvd(d+1,0)/q1.dvd(d+1,0);
    zero : constant quad_double := create(0.0);

  begin
    for i in 0..d+1 loop                   -- normalize table
      for j in 0..d+1-i loop
        q1.dvd(i,j) := q1.dvd(i,j)*fac;
      end loop;
    end loop;
    for i in 0..d loop                     -- compute divided differences
      for j in 0..d-i loop
        q1.dvd(i,j) := (q2.dvd(i,j) - q1.dvd(i,j))/divisor;
      end loop;
    end loop;
    for i in 0..d+1 loop                   -- set higher order ones to zero
      q1.dvd(i,d+1-i) := Create(zero);
    end loop;
  end Adjust_Divided_Differences;

  procedure Adjust_Divided_Differences
               ( file : in file_type; d : in integer32;
                 divisor : in Complex_Number;
                 q2 : in Newton_Taylor_Form1_Rep;
                 q1 : in out Newton_Taylor_Form1_Rep ) is

  -- DESCRIPTION :
  --   Adjustment of the divided differences in q1 using q2.
  --   This version prints the divided difference table and
  --   computes also those higher order divided differences
  --   that should evaluate to zero.

  -- REQUIRED : k = 1, d > 0, and q2 is properly arranged,
  --            q1.dvd(d+1,0) /= 0.

  -- ON ENTRY :
  --   file      for intermediate output;
  --   d         current degree;
  --   divisor   difference between the interpolation points;
  --   q2        interpolator of one degree up;
  --   q1        interpolator whose divided differences will be modified.

  -- ON RETURN :
  --   q1        has adjusted divided differences, up to degree <= d.

    fac : Complex_Number;

  begin
    put_line(file,"Tables of divided differences before normalization :");
    for i in 0..d+1 loop
      for j in 0..d+1-i loop
        put(file,"("); put(file,i,1); put(file,","); put(file,j,1);
        put(file,") at d : "); put(file,q1.dvd(i,j));
        new_line(file);
        put(file,"("); put(file,i,1); put(file,","); put(file,j,1);
        put(file,")  d+1 : "); put(file,q2.dvd(i,j));
        new_line(file);
      end loop;
    end loop;
    fac := q2.dvd(d+1,0)/q1.dvd(d+1,0);
    for i in 0..d+1 loop
      for j in 0..d+1-i loop
        q1.dvd(i,j) := q1.dvd(i,j)*fac;
      end loop;
    end loop;
    put_line(file,"Tables of divided differences after normalization :");
    for i in 0..d+1 loop
      for j in 0..d+1-i loop
        put(file,"("); put(file,i,1); put(file,","); put(file,j,1);
        put(file,") at d : "); put(file,q1.dvd(i,j));
        new_line(file);
        put(file,"("); put(file,i,1); put(file,","); put(file,j,1);
        put(file,")  d+1 : "); put(file,q2.dvd(i,j));
        new_line(file);
      end loop;
    end loop;
    for i in 0..d loop
      for j in 0..d-i loop
        put(file,"("); put(file,i,1); put(file,",");
        put(file,j,1); put(file,") : ");
        q1.dvd(i,j) := (q2.dvd(i,j) - q1.dvd(i,j))/divisor;
        put(file,q1.dvd(i,j)); new_line(file);
      end loop;
    end loop;
    put_line(file,"Test on higher order divided differences : ");
    for i in 0..d+1 loop
      put(file,"("); put(file,i,1); put(file,","); put(file,d+1-i,1);
      put(file,") : ");
      q1.dvd(i,d+1-i) := (q2.dvd(i,d+1-i) - q1.dvd(i,d+1-i))/divisor;
      put(file,q1.dvd(i,d+1-i)); new_line(file);
    end loop;
  end Adjust_Divided_Differences;

  procedure Access_to_Adjust
               ( d : in integer32; divisor : in Complex_Number;
                 upper : in Newton_Taylor_Form;
                 lower : in Newton_Taylor_Form ) is

  -- DESCRIPTION :
  --   With divided differences using the divisor, the tables in lower
  --   will be reduced in degree from combinations with upper.
  --   This procedure mediates the access to the difference tables
  --   at the lowest dimension.

    shift : integer32;

  begin
    if upper.k = 1 then
      Adjust_Divided_Differences(d,divisor,upper.q1,lower.q1);
    else
      shift := upper.qp'last-d;
      for i in shift..upper.qp'last loop
        Access_to_Adjust(i-shift,divisor,upper.qp(i),lower.qp(i));
      end loop;
    end if;
  end Access_to_Adjust;

  procedure Access_to_Adjust
               ( file : in file_type;
                 d : in integer32; divisor : in Complex_Number;
                 upper : in Newton_Taylor_Form;
                 lower : in Newton_Taylor_Form ) is

  -- DESCRIPTION :
  --   With divided differences using the divisor, the tables in lower
  --   will be reduced in degree from combinations with upper.
  --   This procedure mediates the access to the difference tables
  --   at the lowest dimension.  Intermediate output is written to file.

    shift : integer32;

  begin
    put(file,"Access to Adjust for k = "); put(file,upper.k,1);
    put(file," and d = "); put(file,d,1); new_line(file);
    if upper.k = 1 then
      Adjust_Divided_Differences(file,d,divisor,upper.q1,lower.q1);
    else
      shift := upper.qp'last-d;
      put(file,"upper.qp'last : "); put(file,upper.qp'last,1);
      put(file,"  d : "); put(file,d,1);
      put(file,"  shift : "); put(file,shift,1); new_line(file);
      for i in shift..upper.qp'last loop
        put(file,"Recursive call for i = "); put(file,i,1); new_line(file);
        Access_to_Adjust(file,i-shift,divisor,upper.qp(i),lower.qp(i));
      end loop;
    end if;
  end Access_to_Adjust;

  function Get_Constant ( q : Newton_Taylor_Form; d : integer32 )
                        return Complex_Number is

  -- DESCRIPTION :
  --   Returns the constant in the Newton form after adjustments
  --   have been done for degree d.

  begin
    if q.k = 2
     then return q.qp(d).q1.dvd(0,0);
     else return Get_Constant(q.qp(d),d);
    end if;
  end Get_Constant;

  procedure Shift_Truncate
              ( q : in out Newton_Taylor_Form; d : in integer32 ) is

  -- DESCRIPTION :
  --   Truncates the form to degree d.

    res : Newton_Taylor_Form_Rep(q.k,q.n,d);

  begin
    res.rot := q.rot;
    res.pts := q.pts(0..d);
    for i in 1..d loop
      res.qp(i) := q.qp(i+(q.d-d));
    end loop;
    if q.k > 2 then
      for i in 1..d loop
        Shift_Truncate(res.qp(i),i);
      end loop;
    end if;
    declare
      q1rep : Newton_Taylor_Form_Rep(1,q.n,0);
    begin
      q1rep.q1.dvd(0,0) := Get_Constant(q,q.d-d);
      res.qp(0) := new Newton_Taylor_Form_Rep'(q1rep);
    end;
   -- Clear(q);  -- too brutal, causes segmentation fault, sharing
    for i in 0..q.d-d loop
      Clear(q.qp(i));
    end loop;
    q := new Newton_Taylor_Form_Rep'(res);
  end Shift_Truncate;

  procedure Combine_Divided_Differences
              ( k,d : in integer32; pts : in Vector;
                qp : in out Newton_Taylor_Forms ) is

  -- DESCRIPTION :
  --   Combines the divided differences of the Newton forms in qp
  --   such that degree(qp(i)) = i.

    divisor : Complex_Number;

  begin
    for i in 1..d-1 loop                        -- reduce i degrees
      for j in 1..d-i loop
        divisor := pts(d-j-i) - pts(d-j);
        if k = 2 then
          Adjust_Divided_Differences(d-i,divisor,qp(j+1).q1,qp(j).q1);
        else
          for j2 in i..d loop
            Access_to_Adjust(j2-i,divisor,qp(j+1).qp(j2),qp(j).qp(j2));
          end loop;
        end if;
      end loop;
    end loop;
    if k >= 3 then
      for i in 1..d-1 loop
        Shift_Truncate(qp(i),i);
      end loop;
    end if;
  end Combine_Divided_Differences;

  procedure Combine_Divided_Differences
              ( file : in file_type; k,d : in integer32; pts : in Vector;
                qp : in out Newton_Taylor_Forms ) is

  -- DESCRIPTION :
  --   Combines the divided differences of the Newton forms in qp
  --   such that degree(qp(i)) = i.  Writes intermediate output to file.

    divisor : Complex_Number;

  begin
    for i in 1..d-1 loop                        -- reduce i degrees
      for j in 1..d-i loop
        put(file,"Combined Adjust for i = "); put(file,i,1);
        put(file,"  and j = "); put(file,j,1); new_line(file);
        divisor := pts(d-j-i) - pts(d-j);
        if k = 2 then
          put(file,"Adjust for k = "); put(file,k,1);
          put(file,"  and d = "); put(file,d,1); new_line(file);
          Adjust_Divided_Differences(file,d-i,divisor,qp(j+1).q1,qp(j).q1);
        else
          put(file,"Test degrees at j = "); put(file,j,1);
          put_line(file," :"); Test_Degrees(file,qp(j));
          for j2 in i..d loop
            put(file,"Adjust for k = "); put(file,k,1);
            put(file," and j2 = "); put(file,j2,1);
            put(file," and j2-i = "); put(file,j2-i,1);
            put(file," and qp(j).d = "); put(file,qp(j).d,1);
            new_line(file);
            Access_to_Adjust
              (file,j2-i,divisor,qp(j+1).qp(j2),qp(j).qp(j2));
          end loop;
          put(file,"Test degrees after adjust at j = "); put(file,j,1);
          put_line(file," :"); Test_Degrees(file,qp(j));
        end if;
      end loop;
    end loop;
    if k >= 3 then
      for i in 1..d-1 loop
        Shift_Truncate(qp(i),i);
      end loop;
    end if;
  end Combine_Divided_Differences;

  function Truncate ( d : integer32; q : Newton_Taylor_Form1_Rep )
                    return Newton_Taylor_Form1_Rep is

  -- DESCRIPTION :
  --   Truncates the representation of q to a form of degree d.

    res : Newton_Taylor_Form1_Rep(q.n,d);

  begin
    res.rot := q.rot;
    res.x := q.x(0..d);
    res.y := q.y(0..d);
    for i in 0..d loop
      for j in 0..d-i loop
        res.dvd(i,j) := q.dvd(i,j); 
      end loop;
    end loop;
    return res;
  end Truncate;

  function Truncate ( file : file_type;
                      d : integer32; q : Newton_Taylor_Form1_Rep )
                    return Newton_Taylor_Form1_Rep is

  -- DESCRIPTION :
  --   Truncates the representation of q to a form of degree d.
  --   Writes the truncated divided differences to file.

    res : Newton_Taylor_Form1_Rep(q.n,d);

  begin
    put(file,"Truncate on d = "); put(file,d,1);
    put(file," on k = 1 and q.d = "); put(file,q.d,1);
    new_line(file);
    res.rot := q.rot;
    res.x := q.x(0..d);
    res.y := q.y(0..d);
    for i in 0..d loop
      for j in 0..d-i loop
        res.dvd(i,j) := q.dvd(i,j); 
      end loop;
    end loop;
    if d < q.d then
      put_line(file,"Highest order divided differences truncated :");
      for i in d+1..q.d loop
        for j in 0..i loop
          put(file,q.dvd(j,i-j)); new_line(file);
        end loop;
      end loop;
    end if;
    return res;
  end Truncate;

  function Truncate ( d : integer32; q : Newton_Taylor_Form_Rep )
                    return Newton_Taylor_Form_Rep is

  -- DESCRIPTION :
  --   Truncates the representation of q to a form of degree d.

  -- NOTE : the constant res.qp(0) is treated separatedly.

    res : Newton_Taylor_Form_Rep(q.k,q.n,d);

  begin
    if q.k = 1 then
      res.q1 := Truncate(d,q.q1);
    else
      res.rot := q.rot;
      res.pts := q.pts(0..d);
      for i in 1..d loop
        declare
          qp_rep : constant Newton_Taylor_Form_Rep(q.k-1,q.n,i)
                 := Truncate(i,q.qp(i).all);
        begin
          res.qp(i) := new Newton_Taylor_Form_Rep'(qp_rep);
        end;
      end loop;
      res.qp(0) := new Newton_Taylor_Form_Rep'(q.qp(0).all);
    end if;
    return res;
  end Truncate;

  function Truncate ( file : file_type;
                      d : integer32; q : Newton_Taylor_Form_Rep )
                    return Newton_Taylor_Form_Rep is

  -- DESCRIPTION :
  --   Truncates the representation of q to a form of degree d.
  --   Writes truncated divided differences to file.

  -- NOTE : the constant res.qp(0) is treated separatedly.

    res : Newton_Taylor_Form_Rep(q.k,q.n,d);

  begin
    put(file,"Truncate for d = "); put(file,d,1);
    put(file," and q.k = "); put(file,q.k,1); new_line(file);
    if q.k = 1 then
      res.q1 := Truncate(file,d,q.q1);
    else
      res.rot := q.rot;
      res.pts := q.pts(0..d);
      for i in 1..d loop
        declare
          qp_rep : constant Newton_Taylor_Form_Rep(q.k-1,q.n,i)
                 := Truncate(file,i,q.qp(i).all);
        begin
          res.qp(i) := new Newton_Taylor_Form_Rep'(qp_rep);
        end;
      end loop;
      res.qp(0) := new Newton_Taylor_Form_Rep'(q.qp(0).all);
      put_line(file,"Truncated remainders :");
      for i in d+1..q.qp'last loop
        Write_Divided_Differences(file,q.qp(i));
      end loop;
    end if;
    return res;
  end Truncate;

  function Create ( grid : Stacked_Sample_Grid;
                    d : integer32; c : Complex_Number; y : Vector )
                  return Newton_Taylor_Form is

  -- DESCRIPTION :
  --   Recursive call, where d corresponds to the total degree of
  --   the polynomial under construction.

    res : Newton_Taylor_Form;
    res_rep : Newton_Taylor_Form_Rep(grid.k,integer32(grid.n),d);

  begin
    if grid.k = 1 then
      declare
        nq : constant Newton_Interpolator1_Rep
                        (integer32(grid.n),grid.d,grid.d+1)
           := Create(grid.g.all,c);
        tq : constant Newton_Taylor_Form1_Rep(integer32(grid.n),grid.d)
           := Create(nq,y);
      begin
        res_rep.q1.rot := tq.rot;
        res_rep.q1.x := tq.x(0..d);
        res_rep.q1.y := tq.y(0..d);
        for i in 0..d loop
          for j in 0..d loop
            res_rep.q1.dvd(i,j) := tq.dvd(i,j);
          end loop;
        end loop;
      end;
    else
      res_rep.rot := grid.hyp(grid.k)(1..integer32(grid.n));
      res_rep.pts := grid.pts;
      declare
        wrk : Newton_Taylor_Forms(0..grid.d);
      begin
        for i in 1..grid.d loop
          wrk(i) := Create(grid.a(i).all,grid.d,c,y);
        end loop;
        Combine_Divided_Differences(grid.k,grid.d,res_rep.pts,wrk);
        res_rep.qp := wrk;
        res_rep.qp(0)
          := Create_Constant(integer32(grid.n),grid.d,res_rep.pts,grid.spt,wrk);
      end; 
    end if;
    res := new Newton_Taylor_Form_Rep'(res_rep);
    return res;
  end Create;

  function Create ( file : file_type; grid : Stacked_Sample_Grid;
                    d : integer32; c : Complex_Number; y : Vector )
                  return Newton_Taylor_Form is

  -- DESCRIPTION :
  --   Recursive call, where d corresponds to the total degree of
  --   the polynomial under construction.

    res : Newton_Taylor_Form;
    res_rep : Newton_Taylor_Form_Rep(grid.k,integer32(grid.n),d);

  begin
    if grid.k = 1 then
      put(file,"Create at dimension 1, for degree ");
      put(file,grid.d,1); put_line(file,".");
      declare
        nq : constant Newton_Interpolator1_Rep
                        (integer32(grid.n),grid.d,grid.d+1)
           := Create(grid.g.all,c);
        tq : constant Newton_Taylor_Form1_Rep(integer32(grid.n),grid.d)
           := Create(nq,y);
      begin
        Test_Evaluation1(file,grid,tq);
        put(file,"Result is Newton-Taylor form of degree ");
        put(file,d,1); put_line(file,".");
        res_rep.q1.rot := tq.rot;
        res_rep.q1.x := tq.x(0..d);
        res_rep.q1.y := tq.y(0..d);
        for i in 0..d loop
          for j in 0..d loop
            res_rep.q1.dvd(i,j) := tq.dvd(i,j);
          end loop;
        end loop;
      end;
    else
      put(file,"Create at dimension "); put(file,grid.k,1);
      put(file,", for degrees 1 to ");
      put(file,grid.d,1); put_line(file,".");
      res_rep.rot := grid.hyp(grid.k)(1..integer32(grid.n));
      res_rep.pts := grid.pts;
      declare
        wrk : Newton_Taylor_Forms(0..grid.d);
      begin
        for i in 1..grid.d loop
          wrk(i) := Create(file,grid.a(i).all,grid.d,c,y);
        end loop;
        Combine_Divided_Differences(file,grid.k,grid.d,res_rep.pts,wrk);
        res_rep.qp := wrk;
        res_rep.qp(0)
          := Create_Constant(integer32(grid.n),grid.d,res_rep.pts,grid.spt,wrk);
        put_line(file,"Testing evaluation for res_rep.");
        Test_Evaluation(file,grid,res_rep);
      end; 
    end if;
    res := new Newton_Taylor_Form_Rep'(res_rep);
    return res;
  end Create;

  function Create ( grid : Stacked_Sample_Grid; c : Complex_Number;
                    y : Vector ) return Newton_Taylor_Form is

    res,wrk : Newton_Taylor_Form;
    res_rep : Newton_Taylor_Form_Rep(grid.k,integer32(grid.n),grid.d);

  begin
    wrk := Create(grid,grid.d,c,y);
    res_rep := Truncate(grid.d,wrk.all);
    res := new Newton_Taylor_Form_Rep'(res_rep);
    Clear(wrk);
    return res;
  end Create;

  function Create ( file : file_type; grid : Stacked_Sample_Grid;
                    c : Complex_Number; y : Vector )
                  return Newton_Taylor_Form is

    res,wrk : Newton_Taylor_Form;
    res_rep : Newton_Taylor_Form_Rep(grid.k,integer32(grid.n),grid.d);

  begin
    wrk := Create(file,grid,grid.d,c,y);
    put_line(file,"Checking the actual degrees of wrk : ");
    Test_Degrees(file,wrk);
    res_rep := Truncate(file,grid.d,wrk.all);
    res := new Newton_Taylor_Form_Rep'(res_rep);
    Clear(wrk);
    put_line(file,"Checking the actual degrees of res : ");
    Test_Degrees(file,res);
    put_line(file,"Testing evaluation after truncation : ");
    Test_Evaluation(file,grid,res_rep);
    return res;
  end Create;

  function Expand ( q : Newton_Form_Evaluator1 ) return Poly_Sys is

    res : Poly_Sys(0..q.d);

  begin
    for i in res'range loop
      res(i) := Substitute(q.dvd(i).all,q.rot(1),q.rot(2));
    end loop;
    return res;
  end Expand;

  function Expand ( q : Newton_Form_Evaluator1; dvd : Poly_Sys ) return Poly is

    res : Poly;
    extdvd : Poly_Sys(dvd'range) := Extend_and_Trunc(dvd,2);

  begin
    res := Expand(2,extdvd,q.rot,q.pts);
    Clear(extdvd);
    return res;
  end Expand;

  function Expand ( q : Newton_Form_Evaluator1 ) return Poly is

    res : Poly;
    dvd : Poly_Sys(0..q.d) := Expand(q);

  begin
    res := Expand(q,dvd);
    Clear(dvd);
    return res;
  end Expand;

-- SELECTORS :

  function Degree ( q : Newton_Interpolator1 ) return integer32 is
  begin
    return q.d;
  end Degree;

  function Degree ( q : Newton_Form_Evaluator1 ) return integer32 is
  begin
    return q.d;
  end Degree;

  function Dimension ( q : Newton_Interpolator1 ) return integer32 is
  begin
    return q.n;
  end Dimension;

  function Dimension ( q : Newton_Form_Evaluator1 ) return integer32 is
  begin
    return q.n;
  end Dimension;

-- EVALUATORS :

  function Eval ( q : Newton_Interpolator1; x : Vector )
                return Complex_Number is

    zero : constant quad_double := create(0.0);

  begin
    if q = null
     then return Create(zero);
     else return Eval(q.all,x);
    end if;
  end Eval;

  function Eval ( q : Newton_Form_Evaluator1; x : Vector )
                return Complex_Number is

    res : Complex_Number;
    rpx : Vector(1..2);
    zero : constant quad_double := create(0.0);

  begin
    if q = null then
      res := Create(zero);
    else
      rpx := Rotate_and_Project(q.rot,x);
      res := Eval(q.dvd(0).all,rpx(2));
      for i in 1..q.d loop
        res := res*(rpx(1) - q.pts(i));
        res := res + Eval(q.dvd(i).all,rpx(2));
      end loop;
    end if;
    return res;
  end Eval;

  function Eval ( q : Newton_Taylor_Form1; x : Vector )
                return Complex_Number is

    zero : constant quad_double := create(0.0);

  begin
    if q = null
     then return Create(zero);
     else return Eval(q.all,x);
    end if;
  end Eval;

  function Eval ( q : Newton_Taylor_Form; x : Vector )
                return Complex_Number is

    zero : constant quad_double := create(0.0);

  begin
    if q = null
     then return Create(zero);
     else return Eval(q.all,x);
    end if;
  end Eval;

-- DIAGNOSTICS :

  function Errors ( q : Newton_Interpolator1;
                    grid : Array_of_QuadDobl_Sample_Lists ) return Matrix is

    res : Matrix(grid'range,1..integer32(Length_Of(grid(grid'first))));
    tmp : QuadDobl_Sample_List;
    abseva : quad_double;

  begin
    for i in grid'range loop
      tmp := grid(i);
      for j in res'range(2) loop
        abseva := AbsVal(Eval(q,Sample_Point(Head_Of(tmp)).v));
        res(i,j) := hihi_part(abseva);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Errors;

  function Maximal_Error ( residuals : Matrix ) return double_float is

    max : double_float := residuals(residuals'first(1),residuals'first(2));

  begin
    for i in residuals'range(1) loop
      for j in residuals'range(2) loop
        if residuals(i,j) > max
         then max := residuals(i,j);
        end if;
      end loop;
    end loop;
    return max;   
  end Maximal_Error;

  function Maximal_Error 
             ( q : Newton_Interpolator1;
               grid : Array_of_QuadDobl_Sample_Lists ) return double_float is
  begin
    return Maximal_Error(Errors(q,grid));
  end Maximal_Error;

  function Maximal_Error 
             ( q : Newton_Taylor_Form; grid : Stacked_Sample_Grid )
             return double_float is

    res : double_float := -1.0;
    tmp : QuadDobl_Sample_List;
    eva : Complex_Number;
    qd_abseva : quad_double;
    sd_abseva : double_float;

  begin
    if grid.k = 1 then
      for i in grid.g'range loop
        tmp := grid.g(i);
        while not Is_Null(tmp) loop
          eva := Eval(q,Sample_Point(Head_Of(tmp)).v);
          qd_abseva := AbsVal(eva);
          if ((res < 0.0) or else (res < qd_abseva))
           then res := hihi_part(qd_abseva);
          end if;
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    else
      for i in 1..grid.a'last loop
        sd_abseva := Maximal_Error(q,grid.a(i).all);
        if sd_abseva > res
         then res := sd_abseva;
        end if;
      end loop;
      eva := Eval(q,Sample_Point(grid.spt).v);
      qd_abseva := AbsVal(eva);
      if res < qd_abseva
       then res := hihi_part(qd_abseva);
      end if;
    end if;
    return res;
  end Maximal_Error;

-- DESTRUCTORS :

  procedure Clear ( q : in out Newton_Interpolator1 ) is

    procedure free is new unchecked_deallocation
      (Newton_Interpolator1_Rep,Newton_Interpolator1);

  begin
    if q /= null
     then free(q);
    end if;
  end Clear;

  procedure Clear ( q : in out Newton_Form_Evaluator1 ) is

    procedure free is new unchecked_deallocation
      (Newton_Form_Evaluator1_Rep,Newton_Form_Evaluator1);

  begin
    if q /= null
     then Clear(q.dvd); free(q);
    end if;
  end Clear;

  procedure Clear ( q : in out Newton_Taylor_Form1 ) is

    procedure free is new unchecked_deallocation
      (Newton_Taylor_Form1_Rep,Newton_Taylor_Form1);

  begin
    if q /= null
     then free(q);
    end if;
  end Clear;

  procedure Clear ( q : in out Newton_Taylor_Form ) is

    procedure free is new unchecked_deallocation
      (Newton_Taylor_Form_Rep,Newton_Taylor_Form);

  begin
    if q /= null then
      if q.k > 1 then
        for i in q.qp'range loop
          Clear(q.qp(i));
        end loop;
      end if;
      free(q);
    end if;
  end Clear;

end QuadDobl_Divided_Differences;
