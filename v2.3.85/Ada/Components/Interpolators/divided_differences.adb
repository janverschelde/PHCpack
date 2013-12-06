with unchecked_deallocation;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Number_Tools;     use Multprec_Complex_Number_Tools;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Sample_Points;                     use Sample_Points;

package body Divided_Differences is

-- DATA STRUCTURES :

  type Newton_Interpolator_Rep ( n,d,d1 : integer32 ) is record
    strot : Standard_Complex_Vectors.Vector(1..n); 
    mprot : Multprec_Complex_Vectors.Vector(1..n); 
                   -- normal to all parallel hyperplanes in grid
    stpts : Standard_Complex_Matrices.Matrix(0..d,0..d1);      -- data points
    mppts : Multprec_Complex_Matrices.Matrix(0..d,0..d1);      -- data points
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

-- DIAGNOSTICS :

  procedure Difference_pts ( file : in file_type;
                             stpts : Standard_Complex_Matrices.Matrix;
                             mppts : Multprec_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes and prints the difference between the two tables.

    diff : Standard_Complex_Numbers.Complex_Number;
    absdif,max : double_float := 0.0;
    use Standard_Complex_Numbers;

  begin
    put_line(file,"Differences in pts tables :");
    for i in stpts'range(1) loop
      for j in stpts'range(1) loop
        diff := stpts(i,j) - Round(mppts(i,j));
        absdif := AbsVal(diff);
        put(file,absdif,2,3,3);
        if absdif > max
         then max := absdif;
        end if;
      end loop;
      new_line(file);
    end loop;
    put(file,"Maximal difference : "); put(file,max,2,3,3); new_line(file);
  end Difference_pts;

-- NUMERICAL EVALUATION OF THE TABLE OF DIVIDED DIFFERENCES :

  function Rotate_and_Project
             ( v,x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Performs a linear coordinate transformation so that the
  --   first component of the vector on return equals v*x.
  --   The general formula for this linear transformation is
  --     y(1) := v(x'range)*x;
  --     y(i) := v(i)*x(1) - v(1)*x(i), for i in 2..x'last.
  --   Only the first two coordinates of y are needed and returned.

    use Standard_Complex_Numbers,Standard_Complex_Vectors;

    res : Vector(1..2);

  begin
    res(1) := v(x'range)*x;
    res(2) := v(2)*x(1) - v(1)*x(2);
    return res;
  end Rotate_and_Project;

  function Rotate_and_Project
             ( v,x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Performs a linear coordinate transformation so that the
  --   first component of the vector on return equals v*x.
  --   The general formula for this linear transformation is
  --     y(1) := v(x'range)*x;
  --     y(i) := v(i)*x(1) - v(1)*x(i), for i in 2..x'last.
  --   Only the first two coordinates of y are needed and returned.

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;

    res : Vector(1..2);
    acc : Complex_Number;

  begin
    res(1) := v(x'range)*x;
    res(2) := v(2)*x(1);
    acc := v(1)*x(2);
    Sub(res(2),acc);
    Clear(acc);
    return res;
  end Rotate_and_Project;

  function Rotate_and_Project2
             ( v,x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns only the second coordinate of the rotation.

    use Standard_Complex_Numbers;

  begin
    return v(2)*x(1) - v(1)*x(2);
  end Rotate_and_Project2;

  function Rotate_and_Project2
             ( v,x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns only second coordinate of the rotation.

    use Multprec_Complex_Numbers;

    res,acc : Complex_Number;

  begin
    res := v(2)*x(1);
    acc := v(1)*x(2);
    Sub(res,acc);
    Clear(acc);
    return res;
  end Rotate_and_Project2;

  procedure Numerical_Divided_Differences
               ( file : in file_type;
                 stdvd : in out Standard_Complex_Vectors.Vector;
                 mpdvd : in out Multprec_Complex_Vectors.Vector;
                 stx2 : in Standard_Complex_Numbers.Complex_Number;
                 mpx2 : in Multprec_Complex_Numbers.Complex_Number;
                 stpts : in Standard_Complex_Matrices.Matrix;
                 mppts : in Multprec_Complex_Matrices.Matrix;
                 i,j,d : in integer32 ) is

  -- DESCRIPTION :
  --   Updates the vector of divided differences.

  -- ON ENTRY :
  --   dvd       current last row in the table of divided differences;
  --   x2        second coordinate of the data point;
  --   pts       compact matrix form of the interpolation grid;
  --   i         index to current row in table of divided differences;
  --   j         index to current column in table and dvd;
  --   d         degree and size of the table.

  -- ON RETURN :
  --   dvd       dvd(j) will contained updated divided differences.

    use Standard_Complex_Numbers,Multprec_Complex_Numbers;
    acc : Standard_Complex_Numbers.Complex_Number;
    quot,frac : Multprec_Complex_Numbers.Complex_Number;
    acc1,acc2,dif1,dif2 : Multprec_Complex_Numbers.Complex_Number; 
    i0 : boolean;

  begin
    if i = j then
      if i = 0 then
        -- stdvd(j) := Create(1);
        -- mpdvd(j) := Create(1);
        i0 := true;
      else
        stdvd(j) := stpts(i,d+1);
        -- Copy(mppts(i,d+1),mpdvd(j));
        i0 := false;
      end if;
      -- acc1 := Create(1);
      -- acc2 := Create(1);
      acc1 := mpx2 - mppts(i,1);
      acc2 := mppts(0,d+1) - mppts(i,1);
      if i0 then
        stdvd(j) := (stx2 - stpts(i,1))/(stpts(0,d+1) - stpts(i,1));
      else
        acc := (stx2 - stpts(i,1))/(stpts(0,d+1) - stpts(i,1));
        stdvd(j) := stdvd(j)*acc;
      end if;
      for k in 2..d loop
        acc := (stx2 - stpts(i,k))/(stpts(0,d+1) - stpts(i,k));
        stdvd(j) := stdvd(j)*acc;
        -- quot := mpx2 - mppts(i,k);
        -- frac := mppts(0,d+1) - mppts(i,k);
        -- Div(quot,frac);
        -- Mul(mpdvd(j),quot);
        -- Clear(quot); Clear(frac);
        dif1 := mpx2 - mppts(i,k);
        dif2 := mppts(0,d+1) - mppts(i,k);
        Mul(acc1,dif1); Clear(dif1);
        Mul(acc2,dif2); Clear(dif2);
      end loop;
      -- Mul(mpdvd(j),acc1); Clear(acc1);
      Copy(acc1,mpdvd(j));  Clear(acc1);
      if not i0
       then Mul(mpdvd(j),mppts(i,d+1)); Clear(acc1);
      end if;
      Div(mpdvd(j),acc2); Clear(acc2);
    else
      put(file,"ST numerator   : "); put(file,stdvd(j) - stdvd(j+1));
      new_line(file);
      stdvd(j) := (stdvd(j) - stdvd(j+1))/(stpts(j,0) - stpts(i,0));
      quot := mpdvd(j) - mpdvd(j+1);
      frac := mppts(j,0) - mppts(i,0);
      put(file,"MP numerator   : "); put(file,quot); new_line(file);
      put(file,"ST denominator : "); put(file,stpts(j,0) - stpts(i,0));
      new_line(file);
      put(file,"MP denominator : "); put(file,frac); new_line(file);
      -- put(file,"Size numerator   (real part): ");
      -- put(file,Size_Fraction(REAL_PART(quot)),1); new_line(file);
      -- put(file,"Size numerator   (imag part): ");
      -- put(file,Size_Fraction(IMAG_PART(quot)),1); new_line(file);
      -- put(file,"Size denominator (real part): ");
      -- put(file,Size_Fraction(REAL_PART(frac)),1); new_line(file);
      -- put(file,"Size denominator (imag part): ");
      -- put(file,Size_Fraction(IMAG_PART(frac)),1); new_line(file);
      Div(quot,frac);
      Copy(quot,mpdvd(j));
      Clear(quot); Clear(frac);
    end if;
  end Numerical_Divided_Differences;

-- TARGET PROCEDURES :

  function Create ( file : file_type; size : natural32;
                    stgrid : Array_of_Standard_Sample_Lists;
                    mpgrid : Array_of_Multprec_Sample_Lists )
                  return Newton_Interpolator is

    c : constant Standard_Complex_Numbers.Complex_Number := Random1;

  begin
    return Create(file,size,stgrid,mpgrid,c);
  end Create;

  function Create ( file : file_type; size : natural32;
                    stgrid : Array_of_Standard_Sample_Lists;
                    mpgrid : Array_of_Multprec_Sample_Lists;
                    c : Standard_Complex_Numbers.Complex_Number )
                  return Newton_Interpolator is

    res : Newton_Interpolator;
    n : constant integer32 := Number_of_Variables(Head_Of(stgrid(0)));
    d : constant integer32 := stgrid'last;
    sthyp : constant Standard_Complex_VecVecs.VecVec
          := Hyperplane_Sections(Head_Of(stgrid(0)));
    mphyp : constant Multprec_Complex_VecVecs.VecVec
          := Hyperplane_Sections(Head_Of(mpgrid(0)));
    res_rep : Newton_Interpolator_Rep(n,d,d+1);
    sttmp : Standard_Sample_List;
    mptmp : Multprec_Sample_List;
    stspt : Standard_Sample;
    mpspt : Multprec_Sample;
    strpx : Standard_Complex_Vectors.Vector(1..2);
    mpc : Multprec_Complex_Numbers.Complex_Number := Create(c);

    use Standard_Complex_Numbers,Multprec_Complex_Numbers;
    first : boolean;
    acc : Standard_Complex_Numbers.Complex_Number;
    quot,frac : Multprec_Complex_Numbers.Complex_Number;

  begin
    Set_Size(mpc,size);
    res_rep.strot := sthyp(sthyp'first)(1..n);
    Multprec_Complex_Vectors.Copy(mphyp(mphyp'first)(1..n),res_rep.mprot);
    for i in 0..d loop          -- process all samples on same slice #i
      sttmp := stgrid(i);      mptmp := mpgrid(i);
      stspt := Head_Of(sttmp); mpspt := Head_Of(mptmp);
      strpx := Rotate_and_Project(res_rep.strot,Sample_Point(stspt).v);
      res_rep.stpts(i,0) := strpx(1);            -- x coordinate
      res_rep.stpts(i,1) := strpx(2);            -- y coordinate
      declare
        mprpx : constant Multprec_Complex_Vectors.Vector(1..2)
              := Rotate_and_Project(res_rep.mprot,Sample_Point(mpspt).v);
      begin
        res_rep.mppts(i,0) := mprpx(1);
        res_rep.mppts(i,1) := mprpx(2);
      end;
      sttmp := Tail_Of(sttmp); mptmp := Tail_Of(mptmp);
      for j in 2..d loop
        stspt := Head_Of(sttmp); mpspt := Head_Of(mptmp);
        res_rep.stpts(i,j)
          := Rotate_and_Project2(res_rep.strot,Sample_Point(stspt).v);
        res_rep.mppts(i,j)
          := Rotate_and_Project2(res_rep.mprot,Sample_Point(mpspt).v);
        sttmp := Tail_Of(sttmp); mptmp := Tail_Of(mptmp);
      end loop;
    end loop;
    res_rep.stpts(0,d+1) := c;
    res_rep.mppts(0,d+1) := mpc;
    for i in 1..d loop
      first := true;
      for j in 1..d loop                       -- compute p(x_i,c)
        if first then
          res_rep.stpts(i,d+1)
            := (c - res_rep.stpts(i,j))/(c - res_rep.stpts(0,j));
          quot := mpc - res_rep.mppts(i,j);
          frac := mpc - res_rep.mppts(0,j);
          res_rep.mppts(i,d+1) := quot/frac;
        else
          first := false;
          acc := (c - res_rep.stpts(i,j))/(c - res_rep.stpts(0,j));
          res_rep.stpts(i,d+1) := res_rep.stpts(i,d+1)*acc;
          quot := mpc - res_rep.mppts(i,j);
          frac := mpc - res_rep.mppts(0,j);
          Div(quot,frac);
          Mul(res_rep.mppts(i,d+1),quot);
        end if;
        Clear(quot); Clear(frac);
      end loop;
    end loop;
    put_line(file,"Sizes of the numbers in pts : ");
    for i in res_rep.mppts'range(1) loop
      for j in res_rep.mppts'range(2) loop
        put(file,Size_Fraction(REAL_PART(res_rep.mppts(i,j))),2);
      end loop;
      new_line(file);
    end loop;
    Difference_pts(file,res_rep.stpts,res_rep.mppts);
    res := new Newton_Interpolator_Rep'(res_rep);
    return res;
  end Create;

-- EVALUATOR :

  function Maximal_Difference
                 ( stdvd : Standard_Complex_Vectors.Vector;
                   mpdvd : Multprec_Complex_Vectors.Vector )
                 return double_float is

  -- DESCRIPTION :
  --   Returns the largest absolute value of all differences between
  --   the vectors stdvd and mpdvd.

    max : double_float;
    use Standard_Complex_Numbers;
    dif : Complex_Number;

  begin
    dif := stdvd(stdvd'first) - Round(mpdvd(mpdvd'first));
    max := AbsVal(dif);
    for i in stdvd'first+1..stdvd'last loop
      dif := stdvd(i) - Round(mpdvd(i));
      if AbsVal(dif) > max
       then max := AbsVal(dif);
      end if;
    end loop;
    return max;
  end Maximal_Difference;

  procedure Eval ( file : in file_type; q : in Newton_Interpolator;
                   stx : in Standard_Complex_Vectors.Vector;
                   mpx : in Multprec_Complex_Vectors.Vector;
                   steva : out Standard_Complex_Numbers.Complex_Number;
                   mpeva : out Multprec_Complex_Numbers.Complex_Number ) is

    stdvd : Standard_Complex_Vectors.Vector(0..q.d);
    mpdvd : Multprec_Complex_Vectors.Vector(0..q.d);
    strpx : Standard_Complex_Vectors.Vector(1..2)
          := Rotate_and_Project(q.strot,stx);
    mprpx : Multprec_Complex_Vectors.Vector(1..2)
          := Rotate_and_Project(q.mprot,mpx);

    use Standard_Complex_Numbers,Multprec_Complex_Numbers;
    acc : Multprec_Complex_Numbers.Complex_Number;
    dif : Standard_Complex_Numbers.Complex_Number;
    absdif : double_float;

  begin
    put_line(file,"The table of divided differences : ");
    for i in 0..stdvd'last loop
      for j in reverse 0..i loop
        Numerical_Divided_Differences
          (file,stdvd,mpdvd,strpx(2),mprpx(2),q.stpts,q.mppts,i,j,q.d);
        put(file,"ST dvd(");
        put(file,i,1); put(file,","); put(file,j,1); put(file,") : ");
        put(file,stdvd(j)); new_line(file);
        put(file,"MP dvd(");
        put(file,i,1); put(file,","); put(file,j,1); put(file,") : ");
        put(file,mpdvd(j)); new_line(file);
        dif := stdvd(j) - Round(mpdvd(j));
        put(file,"dvd difference : "); put(file,dif); new_line(file);
        absdif := AbsVal(dif);
        if absdif > 0.01 then
          put(file,"AbsVal(dvd difference) : ");
          put(file,absdif,3); put_line(file," > 0.01 BUG!!!");
        end if;
      end loop;
    end loop;
    put(file,"The largest difference ST-MP : ");
    put(file,Maximal_Difference(stdvd,mpdvd),3); new_line(file);
    steva := stdvd(0);
    Copy(mpdvd(0),mpeva);
    put(file,"ST eva(0) : "); put(file,steva); new_line(file);
    put(file,"MP eva(0) : "); put(file,mpeva); new_line(file);
    for i in 1..q.d loop
      steva := steva*(strpx(1) - q.stpts(i,0));
      steva := steva + stdvd(i);
      acc := mprpx(1) - q.mppts(i,0);
      Mul(mpeva,acc);
      Clear(acc);
      Add(mpeva,mpdvd(i));
      put(file,"ST eva("); put(file,i,1); put(file,") : ");
      put(file,steva); new_line(file);
      put(file,"MP eva("); put(file,i,1); put(file,") : ");
      put(file,mpeva); new_line(file);
      dif := steva - Round(mpeva);
      put(file,"eva difference : "); put(file,dif); new_line(file);
      absdif := AbsVal(dif);
      if absdif > 0.01 then
        put(file,"AbsVal(eva difference) : "); put(file,absdif,3);
        put_line(file," > 0.01 BUG!!!");
      end if;
    end loop;
    Multprec_Complex_Vectors.Clear(mpdvd);
    Multprec_Complex_Vectors.Clear(mprpx);
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( q : in out Newton_Interpolator ) is

    procedure free is new unchecked_deallocation
      (Newton_Interpolator_Rep,Newton_Interpolator);

  begin
    if q /= null
     then free(q);
    end if;
  end Clear;

end Divided_Differences;
