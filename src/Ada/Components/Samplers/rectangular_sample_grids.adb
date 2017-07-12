with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers_Polar;      use Standard_Complex_Numbers_Polar;
with Standard_Complex_Norms_Equals;       use Standard_Complex_Norms_Equals;
with Standard_Random_Numbers;             use Standard_Random_Numbers;
with Extended_Random_Numbers;             use Extended_Random_Numbers;
with Multprec_Complex_Number_Tools;       use Multprec_Complex_Number_Tools;
with Multprec_Complex_Numbers_Polar;      use Multprec_Complex_Numbers_Polar;
with Multprec_Complex_Norms_Equals;       use Multprec_Complex_Norms_Equals;

package body Rectangular_Sample_Grids is

-- WARNING :
--   roots of unity conflicted with multi-precision divided differences!

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean ) is
  begin
    Sample_Points.Set_Polynomial_Type(laurent);
    Sample_Point_Lists.Set_Polynomial_Type(laurent);
  end Set_Polynomial_Type;

-- ROOTS OF UNITY :

  function Standard_Roots_of_Unity
             ( d : natural32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with the d complex roots of unity.

    res : Standard_Complex_Vectors.Vector(1..integer32(d));
    use Standard_Complex_Numbers;
    one : constant Complex_Number := Create(1.0);

  begin
    for i in 1..d loop
      res(integer32(i)) := Root(one,d,i);
    end loop;
    return res;
  end Standard_Roots_of_Unity;

  function Multprec_Roots_of_Unity
             ( d,size : natural32 ) return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with the d complex roots of unity.

    res : Multprec_Complex_Vectors.Vector(1..integer32(d));
    use Multprec_Complex_Numbers;
    one : Complex_Number := Multprec_Complex_Numbers.Create(integer(1));

  begin
    Set_Size(one,size);
    for i in 1..d loop
      res(integer32(i)) := Root(one,d,i);
    end loop;
    Clear(one);
    return res;
  end Multprec_Roots_of_Unity;

-- UTILITIES :

  function Extended_Random
             ( v : Standard_Complex_Vectors.Vector; size : natural32 )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);
    tol : constant double_float := 1.0E-12;
    mpvi : Multprec_Complex_Numbers.Complex_Number;

  begin
    for i in res'range loop
      mpvi := Create(v(i));
      if Standard_Complex_Numbers.AbsVal(v(i)) < tol
       then res(i) := mpvi;
       else res(i) := Extended_Random(mpvi,size);
            Multprec_Complex_Numbers.Clear(mpvi);
      end if;
    end loop;
    return res;
  end Extended_Random;

  function Extended_Random
             ( v : Standard_Complex_VecVecs.VecVec; size : natural32 )
             return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      res(i) := new Multprec_Complex_Vectors.Vector'
                      (Extended_Random(v(i).all,size));
    end loop;
    return res;
  end Extended_Random;

-- CREATORS :

  function Create1 ( sps : Standard_Sample_List; m : natural32 )
                   return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..integer32(m));
    sli : constant Standard_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    roots : constant Standard_Complex_Vectors.Vector(1..integer32(m))
          := Standard_Roots_of_Unity(m);

  begin
    res(0) := sps;
    for i in 1..res'last loop
      declare
        newsli : Standard_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : Standard_Sample_List;
      begin
        for j in sli'range loop
          newsli(j) := new Standard_Complex_Vectors.Vector'(sli(j).all);
        end loop;
        newsli(1)(0) := roots(i);
        Sample(sps,newsli,newsps,newsps_last);
        res(i) := newsps;
      end;   
    end loop;
    return res;
  end Create1;

  function Create1 ( sps : Standard_Sample_List; m,size : natural32 )
                   return Array_of_Multprec_Sample_Lists is

    res : Array_of_Multprec_Sample_Lists(0..integer32(m));
    sli : constant Standard_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    extsli : constant Multprec_Complex_VecVecs.VecVec(sli'range)
           := Extended_Random(sli,size);
   -- stroots : Standard_Complex_Vectors.Vector(1..m)
   --         := Standard_Roots_of_Unity(m);
   -- mproots : Multprec_Complex_Vectors.Vector(1..m)
   --         := Multprec_Roots_of_Unity(m,size);
    grid0_last : Multprec_Sample_List;
    sps1 : Standard_Sample_List := sps;
    mpc : Multprec_Complex_Numbers.Complex_Number;

  begin
    Refine_on_Slices(sps1,extsli,res(0),grid0_last);
   -- Refine(sps1,res(0),grid0_last);
    for i in 1..integer32(Length_Of(sps)) loop
      declare
        newstsli : Standard_Complex_VecVecs.VecVec(sli'range);
        newmpsli : Multprec_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : Multprec_Sample_List;
      begin
        for j in sli'range loop
          newstsli(j) := new Standard_Complex_Vectors.Vector'(sli(j).all);
          newmpsli(j) := new Multprec_Complex_Vectors.Vector(extsli(j)'range);
          Multprec_Complex_Vectors.Copy(extsli(j).all,newmpsli(j).all);
        end loop;
        newstsli(1)(0) := Random1;
        mpc := Create(newstsli(1)(0));
        Multprec_Complex_Numbers.Clear(newmpsli(1)(0));
        newmpsli(1)(0) := Extended_Random(mpc,size);
        Multprec_Complex_Numbers.Clear(mpc);
        Sample_on_Slices(sps1,newstsli,newmpsli,newsps,newsps_last);
       -- Sample(sps1,newstsli,newsps,newsps_last);
        res(i) := newsps;
      end;
    end loop;
   -- Multprec_Complex_Vectors.Clear(mproots);
    return res;
  end Create1;

  procedure Create1 ( sps : in Standard_Sample_List; m,size : natural32;
                      stgrid : out Array_of_Standard_Sample_Lists;
                      mpgrid : out Array_of_Multprec_Sample_Lists ) is

    sli : constant Standard_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    extsli : constant Multprec_Complex_VecVecs.VecVec(sli'range)
           := Extended_Random(sli,size);
    stroots : constant Standard_Complex_Vectors.Vector(1..integer32(m))
            := Standard_Roots_of_Unity(m);
    mproots : constant Multprec_Complex_Vectors.Vector(1..integer32(m))
            := Multprec_Roots_of_Unity(m,size);
    mpgrid0_last : Multprec_Sample_List;

  begin
    stgrid(0) := sps;
   -- Refine(stgrid(0),mpgrid(0),mpgrid0_last);
    Refine_on_Slices(stgrid(0),extsli,mpgrid(0),mpgrid0_last);
    for i in 1..integer32(m) loop
      declare
        newstsli : Standard_Complex_VecVecs.VecVec(sli'range);
        newmpsli : Multprec_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : Standard_Sample_List;
        refsps,refsps_last : Multprec_Sample_List;
      begin
        for j in sli'range loop
          newstsli(j) := new Standard_Complex_Vectors.Vector'(sli(j).all);
          newstsli(j)(0) := stroots(i); 
          newmpsli(j) := new Multprec_Complex_Vectors.Vector(extsli(j)'range);
          Multprec_Complex_Vectors.Copy(extsli(j).all,newmpsli(j).all);
          Multprec_Complex_Numbers.Copy(mproots(i),newmpsli(j)(0)); 
        end loop;
        Sample(sps,newstsli,newsps,newsps_last);
        stgrid(i) := newsps;
       -- Refine(newsps,refsps,refsps_last);
        Refine_on_Slices(newsps,newmpsli,refsps,refsps_last);
        mpgrid(i) := refsps;
      end;
    end loop;
  end Create1;

  function Triangular_Create1 ( sps : Standard_Sample_List; m : natural32 )
                              return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..integer32(m));
    sli : constant Standard_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    roots : constant Standard_Complex_Vectors.Vector(1..integer32(m))
          := Standard_Roots_of_Unity(m);
    wrksps : Standard_Sample_List;

  begin
    res(0) := sps;
    wrksps := sps;
    for i in 1..integer32(m) loop
      declare
        newsli : Standard_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : Standard_Sample_List;
      begin
        for j in sli'range loop
          newsli(j) := new Standard_Complex_Vectors.Vector'(sli(j).all);
        end loop;
        newsli(1)(0) := roots(i);
        Sample(wrksps,newsli,newsps,newsps_last);
        res(i) := newsps;
      end;
      wrksps := Tail_Of(wrksps);
    end loop;
    return res;
  end Triangular_Create1;

  function Triangular_Create1
             ( sps : Standard_Sample_List; m,size : natural32 )
             return Array_of_Multprec_Sample_Lists is

    res : Array_of_Multprec_Sample_Lists(0..integer32(m));
    sli : constant Standard_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    wrksps : Standard_Sample_List;
    extsli : constant Multprec_Complex_VecVecs.VecVec(sli'range)
           := Extended_Random(sli,size);
    grid0_last : Multprec_Sample_List;
    sps1 : Standard_Sample_List := sps;
    mpc : Multprec_Complex_Numbers.Complex_Number;

  begin
    Refine_on_Slices(sps1,extsli,res(0),grid0_last);
    wrksps := sps1;
    for i in 1..integer32(Length_Of(sps)) loop
      declare
        newstsli : Standard_Complex_VecVecs.VecVec(sli'range);
        newmpsli : Multprec_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : Multprec_Sample_List;
      begin
        for j in sli'range loop
          newstsli(j) := new Standard_Complex_Vectors.Vector'(sli(j).all);
          newmpsli(j) := new Multprec_Complex_Vectors.Vector(extsli(j)'range);
          Multprec_Complex_Vectors.Copy(extsli(j).all,newmpsli(j).all);
        end loop;
        newstsli(1)(0) := Random1;
        mpc := Create(newstsli(1)(0));
        Multprec_Complex_Numbers.Clear(newmpsli(1)(0));
        newmpsli(1)(0) := Extended_Random(mpc,size);
        Multprec_Complex_Numbers.Clear(mpc);
        Sample_on_Slices(wrksps,newstsli,newmpsli,newsps,newsps_last);
        res(i) := newsps;
      end;
      wrksps := Tail_Of(wrksps);
    end loop;
    return res;
  end Triangular_Create1;

-- DIAGNOSTICS :

  function Errors ( grid : Array_of_Standard_Sample_Lists ) 
                  return Standard_Floating_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Standard_Floating_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : Standard_Sample_List;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := 0.0;
      end loop;
    end loop;
    for i in grid'range loop
      tmp := grid(i);
      for j in 1..integer32(len) loop
        res(i,j) := Sample_Point(Head_Of(tmp)).err;
        tmp := Tail_Of(tmp);
        exit when Is_Null(tmp);
      end loop;
    end loop;
    return res;
  end Errors;

  function Errors ( grid : Array_of_Multprec_Sample_Lists ) 
                  return Multprec_Floating_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Multprec_Floating_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : Multprec_Sample_List;

  begin
    for i in grid'range loop
      tmp := grid(i);
      for j in 1..integer32(len) loop
        res(i,j) := Sample_Point(Head_Of(tmp)).err;
        tmp := Tail_Of(tmp);
        exit when Is_Null(tmp);
      end loop;
    end loop;
    return res;
  end Errors;

  function Maximal_Error ( grid_errors : Standard_Floating_Matrices.Matrix ) 
                         return double_float is

    max : double_float
        := grid_errors(grid_errors'first(1),grid_errors'first(2));

  begin
    for i in grid_errors'range(1) loop
      for j in grid_errors'range(2) loop
        if grid_errors(i,j) > max
         then max := grid_errors(i,j);
        end if;
      end loop;
    end loop;
    return max;
  end Maximal_Error;

  function Maximal_Error ( grid_errors : Multprec_Floating_Matrices.Matrix ) 
                         return Floating_Number is

    max : Floating_Number;

  begin
    Copy(grid_errors(grid_errors'first(1),grid_errors'first(2)),max);
    for i in grid_errors'range(1) loop
      for j in grid_errors'range(2) loop
        if grid_errors(i,j) > max
         then Copy(grid_errors(i,j),max);
        end if;
      end loop;
    end loop;
    return max;
  end Maximal_Error;

  function Maximal_Error ( grid : Array_of_Standard_Sample_Lists )
                         return double_float is
  begin
    return Maximal_Error(Errors(grid));
  end Maximal_Error;

  function Maximal_Error ( grid : Array_of_Multprec_Sample_Lists )
                         return Floating_Number is
  begin
    return Maximal_Error(Errors(grid));
  end Maximal_Error;

  function Distance ( spt1,spt2 : Standard_Sample ) return double_float is

    use Standard_Complex_Vectors;

  begin
    return Max_Norm(Sample_Point(spt1).v - Sample_Point(spt2).v);
  end Distance;

  function Distance ( spt1,spt2 : Multprec_Sample ) return Floating_Number is

    res : Floating_Number;
    use Multprec_Complex_Vectors;
    diff : Vector(1..Number_of_Variables(spt1));

  begin
    diff := Sample_Point(spt1).v - Sample_Point(spt2).v;
    res := Max_Norm(diff);
    Clear(diff);
    return res;
  end Distance;

  function Distance ( sps : Standard_Sample_List; i : natural32;
                      spt : Standard_Sample ) return double_float is

    res,dst : double_float;
    tmp : Standard_Sample_List := sps;
    first : boolean := true;

  begin
    for j in 1..Length_Of(sps) loop
      if j /= i then
        dst := Distance(Head_Of(tmp),spt);
        if first then
          res := dst; first := false;
        elsif dst < res then
          res := dst;
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Distance;

  function Distance ( sps : Multprec_Sample_List; i : natural32;
                      spt : Multprec_Sample ) return Floating_Number is

    res,dst : Floating_Number;
    tmp : Multprec_Sample_List := sps;
    first : boolean := true;

  begin
    for j in 1..Length_Of(sps) loop
      if j /= i then
        dst := Distance(Head_Of(tmp),spt);
        if first then
          Copy(dst,res); first := false;
        elsif dst < res then
          Copy(dst,res);
        end if;
        Clear(dst);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Distance;

  function Distances ( grid : Array_of_Standard_Sample_Lists )
                     return Standard_Floating_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Standard_Floating_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : Standard_Sample_List;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := 1.0E+50;
      end loop;
    end loop;
    for i in grid'range loop
      tmp := grid(i);
      if not Is_Null(Tail_Of(tmp)) then    -- implies Length_Of(tmp) > 1
        for j in 1..integer32(len) loop
          res(i,j) := Distance(grid(i),natural32(j),Head_Of(tmp));
          tmp := Tail_Of(tmp);
          exit when Is_Null(tmp);
        end loop;
      end if;
    end loop; 
    return res;
  end Distances;

  function Distances ( grid : Array_of_Multprec_Sample_Lists )
                     return Multprec_Floating_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Multprec_Floating_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : Multprec_Sample_List;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(1.0E+50);
      end loop;
    end loop;
    for i in grid'range loop
      tmp := grid(i);
      if not Is_Null(Tail_Of(tmp)) then     -- implies Length_Of(tmp) > 1
        for j in 1..integer32(len) loop
          res(i,j) := Distance(grid(i),natural32(j),Head_Of(tmp));
          tmp := Tail_Of(tmp);
          exit when Is_Null(tmp);
        end loop;
      end if;
    end loop; 
    return res;
  end Distances;

  function Minimal_Distance ( grid_dist : Standard_Floating_Matrices.Matrix )
                            return double_float is

    min : double_float
        := grid_dist(grid_dist'first(1),grid_dist'first(2));

  begin
    for i in grid_dist'range(1) loop
      for j in grid_dist'range(2) loop
        if grid_dist(i,j) < min
         then min := grid_dist(i,j);
        end if;
      end loop;
    end loop;
    return min;
  end Minimal_Distance;

  function Minimal_Distance ( grid_dist : Multprec_Floating_Matrices.Matrix )
                            return Floating_Number is

    min : Floating_Number;

  begin
    Copy(grid_dist(grid_dist'first(1),grid_dist'first(2)),min);
    for i in grid_dist'range(1) loop
      for j in grid_dist'range(2) loop
        if grid_dist(i,j) < min
         then Copy(grid_dist(i,j),min);
        end if;
      end loop;
    end loop;
    return min;
  end Minimal_Distance;

  function Minimal_Distance ( grid : Array_of_Standard_Sample_Lists )
                            return double_float is
  begin
    return Minimal_Distance(Distances(grid));
  end Minimal_Distance;

  function Minimal_Distance ( grid : Array_of_Multprec_Sample_Lists )
                            return Floating_Number is

    res : Floating_Number;
    dist : Multprec_Floating_Matrices.Matrix
             (grid'range,1..integer32(Length_Of(grid(grid'first))))
         := Distances(grid);

  begin
    res := Minimal_Distance(dist);
    Multprec_Floating_Matrices.Clear(dist);
    return res;
  end Minimal_Distance;

-- SELECTORS :

  function Abscisses
             ( grid : Array_of_Standard_Sample_Lists; i : natural32 )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(0..integer32(i));

    use Standard_Complex_Numbers;

  begin
    for j in 0..integer32(i) loop
      declare
        hyp : constant Standard_Complex_VecVecs.VecVec
            := Hyperplane_Sections(Head_Of(grid(j)));
      begin
        res(j) := -hyp(hyp'first)(0);
      end;
    end loop;
    return res;
  end Abscisses;

  function Abscisses
             ( grid : Array_of_Multprec_Sample_Lists; i : natural32 )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(0..integer32(i));

    use Multprec_Complex_Numbers;

  begin
    for j in 0..integer32(i) loop
      declare
        hyp : constant Multprec_Complex_VecVecs.VecVec
            := Hyperplane_Sections(Head_Of(grid(j)));
      begin
        res(j) := -hyp(hyp'first)(0);
      end;
    end loop;
    return res;
  end Abscisses;

  function Extract_Samples
                ( grid : Array_of_Standard_Sample_Lists )
                return Standard_Complex_VecVecs.VecVec is

    d : constant integer32 := grid'last;
    res : Standard_Complex_VecVecs.VecVec(1..d*(d+1));
    tmp : Standard_Sample_List;
    spt : Standard_Sample;
    ind : integer32 := 1;
    rpx : Standard_Complex_Vectors.Vector(1..2);

    use Standard_Complex_Numbers,Standard_Complex_Vectors;

  begin
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        spt := Head_Of(tmp);
        rpx(1) := Sample_Point(spt).v(1);
        rpx(2) := Sample_Point(spt).v(2); 
        res(ind) := new Standard_Complex_Vectors.Vector'(rpx);
        tmp := Tail_Of(tmp);
        ind := ind+d;
        if ind > res'last then
          ind := ind mod res'last;
          ind := ind+1;
        end if;
      end loop;
    end loop;
    return res;
  end Extract_Samples;

  function Extract_Samples
                ( grid : Array_of_Multprec_Sample_Lists )
                return Multprec_Complex_VecVecs.VecVec is

    d : constant integer32 := grid'last;
    res : Multprec_Complex_VecVecs.VecVec(1..d*(d+1));
    tmp : Multprec_Sample_List;
    spt : Multprec_Sample;
    ind : integer32 := 1;

    use Multprec_Complex_Numbers,Standard_Complex_Vectors;

  begin
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        spt := Head_Of(tmp);
        declare
          rpx : Multprec_Complex_Vectors.Vector(1..2);
        begin
          Copy(Sample_Point(spt).v(1),rpx(1));
          Copy(Sample_Point(spt).v(2),rpx(2)); 
          res(ind) := new Multprec_Complex_Vectors.Vector'(rpx);
        end;
        tmp := Tail_Of(tmp);
        ind := ind+d;
        if ind > res'last then
          ind := ind mod res'last;
          ind := ind+1;
        end if;
      end loop;
    end loop;
    return res;
  end Extract_Samples;

  function Rotate_and_Project ( v,x : Standard_Complex_Vectors.Vector )
                              return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(1..2);

    use Standard_Complex_Numbers,Standard_Complex_Vectors;

  begin
    res(1) := v(x'range)*x;
    res(2) := v(2)*x(1) - v(1)*x(2);
    return res;
  end Rotate_and_Project;

  function Rotate_and_Project ( v,x : Multprec_Complex_Vectors.Vector )
                              return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(1..2);
    acc : Multprec_Complex_Numbers.Complex_Number;

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;

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

    use Standard_Complex_Numbers;

  begin
    return v(2)*x(1) - v(1)*x(2);
  end Rotate_and_Project2;

  function Rotate_and_Project2
               ( v,x : Multprec_Complex_Vectors.Vector )
               return Multprec_Complex_Numbers.Complex_Number is

    res,acc : Multprec_Complex_Numbers.Complex_Number;

    use Multprec_Complex_Numbers;

  begin
    res := v(2)*x(1);
    acc := v(1)*x(2);
    Sub(res,acc);
    Clear(acc);
    return res;
  end Rotate_and_Project2;

  function Inverse_Rotate ( v,z : Standard_Complex_Vectors.Vector )
                          return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers,Standard_Complex_Vectors;

    res : Vector(1..2);
    det : Complex_Number := v(1)*v(1) + v(2)*v(2);

  begin
    res(1) := (z(2)*v(2) + z(1)*v(1))/det;
    res(2) := (v(2)*z(1) - v(1)*z(2))/det;
    return res;
  end Inverse_Rotate;

  function Inverse_Rotate ( v,z : Multprec_Complex_Vectors.Vector )
                          return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(1..2);
    acc,det : Multprec_Complex_Numbers.Complex_Number;

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;

  begin
    det := v(1)*v(1);
    acc := v(2)*v(2);
    Add(det,acc); Clear(acc);
    res(1) := z(2)*v(2);
    acc := z(1)*v(1);
    Add(res(1),acc); Clear(acc);
    Div(res(1),det);
    res(2) := v(2)*z(1);
    acc := v(1)*z(2);
    Sub(res(2),acc); Clear(acc);
    Div(res(2),det);
    Clear(det);
    return res;
  end Inverse_Rotate;

  function Rotate_Samples ( d,k : natural32;
                            rot : Standard_Complex_Vectors.Vector;
                            grid : Array_of_Standard_Sample_Lists )
                          return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..integer32(d),0..integer32(k));
    tmp : Standard_Sample_List;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Standard_Complex_Numbers.Create(0.0);
      end loop;
    end loop;
    for j in 0..integer32(k) loop
      tmp := grid(j);
      for i in 1..integer32(d) loop
        res(i,j) := Rotate_and_Project2(rot,Sample_Point(Head_Of(tmp)).v);
       -- res(i,j) := Sample_Point(Head_Of(tmp)).v(2);
        tmp := Tail_Of(tmp);
        exit when Is_Null(tmp);
      end loop;
    end loop;
    return res;
  end Rotate_Samples;

  function Rotate_Samples ( d,k : natural32;
                            rot : Multprec_Complex_Vectors.Vector;
                            grid : Array_of_Multprec_Sample_Lists )
                          return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(1..integer32(d),0..integer32(k));
    tmp : Multprec_Sample_List;

  begin
    for j in 0..integer32(k) loop
      tmp := grid(j);
      for i in 1..integer32(d) loop
        res(i,j) := Rotate_and_Project2(rot,Sample_Point(Head_Of(tmp)).v);
        tmp := Tail_Of(tmp);
        exit when Is_Null(tmp);
      end loop;
    end loop;
    return res;
  end Rotate_Samples;

end Rectangular_Sample_Grids;
