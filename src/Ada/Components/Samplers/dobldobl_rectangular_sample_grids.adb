with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with DoblDobl_Complex_Numbers_Polar;      use DoblDobl_Complex_Numbers_Polar;
with DoblDobl_Complex_Vector_Norms;       use DoblDobl_Complex_Vector_Norms;

package body DoblDobl_Rectangular_Sample_Grids is

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean ) is
  begin
    DoblDobl_Sample_Points.Set_Polynomial_Type(laurent);
    DoblDobl_Sample_Lists.Set_Polynomial_Type(laurent);
  end Set_Polynomial_Type;

-- ROOTS OF UNITY :

  function DoblDobl_Roots_of_Unity
             ( d : natural32 ) return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with the d complex roots of unity.

    res : DoblDobl_Complex_Vectors.Vector(1..integer32(d));
    use DoblDobl_Complex_Numbers;
    dd_one : constant double_double := create(1.0);
    one : constant Complex_Number := Create(dd_one);

  begin
    for i in 1..d loop
      res(integer32(i)) := Root(one,d,i);
    end loop;
    return res;
  end DoblDobl_Roots_of_Unity;

-- CREATORS :

  function Create1 ( sps : DoblDobl_Sample_List; m : natural32 )
                   return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..integer32(m));
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    roots : constant DoblDobl_Complex_Vectors.Vector(1..integer32(m))
          := DoblDobl_Roots_of_Unity(m);

  begin
    res(0) := sps;
    for i in 1..res'last loop
      declare
        newsli : DoblDobl_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : DoblDobl_Sample_List;
      begin
        for j in sli'range loop
          newsli(j) := new DoblDobl_Complex_Vectors.Vector'(sli(j).all);
        end loop;
        newsli(1)(0) := roots(i);
        Sample(sps,newsli,newsps,newsps_last);
        res(i) := newsps;
      end;   
    end loop;
    return res;
  end Create1;

  function Triangular_Create1 ( sps : DoblDobl_Sample_List; m : natural32 )
                              return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..integer32(m));
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    roots : constant DoblDobl_Complex_Vectors.Vector(1..integer32(m))
          := DoblDobl_Roots_of_Unity(m);
    wrksps : DoblDobl_Sample_List;

  begin
    res(0) := sps;
    wrksps := sps;
    for i in 1..integer32(m) loop
      declare
        newsli : DoblDobl_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : DoblDobl_Sample_List;
      begin
        for j in sli'range loop
          newsli(j) := new DoblDobl_Complex_Vectors.Vector'(sli(j).all);
        end loop;
        newsli(1)(0) := roots(i);
        Sample(wrksps,newsli,newsps,newsps_last);
        res(i) := newsps;
      end;
      wrksps := Tail_Of(wrksps);
    end loop;
    return res;
  end Triangular_Create1;

-- DIAGNOSTICS :

  function Errors ( grid : Array_of_DoblDobl_Sample_Lists ) 
                  return Double_Double_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Double_Double_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : DoblDobl_Sample_List;
    zero : constant double_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := zero;
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

  function Maximal_Error ( grid_errors : Double_Double_Matrices.Matrix ) 
                         return double_double is

    max : double_double
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

  function Maximal_Error ( grid : Array_of_DoblDobl_Sample_Lists )
                         return double_double is
  begin
    return Maximal_Error(Errors(grid));
  end Maximal_Error;

  function Distance ( spt1,spt2 : DoblDobl_Sample ) return double_double is

    use DoblDobl_Complex_Vectors;

  begin
    return Max_Norm(Sample_Point(spt1).v - Sample_Point(spt2).v);
  end Distance;

  function Distance ( sps : DoblDobl_Sample_List; i : natural32;
                      spt : DoblDobl_Sample ) return double_double is

    res,dst : double_double;
    tmp : DoblDobl_Sample_List := sps;
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

  function Distances ( grid : Array_of_DoblDobl_Sample_Lists )
                     return Double_Double_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Double_Double_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : DoblDobl_Sample_List;
    init : constant double_double := create(1.0E+50);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := init;
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

  function Minimal_Distance ( grid_dist : Double_Double_Matrices.Matrix )
                            return double_double is

    min : double_double
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

  function Minimal_Distance ( grid : Array_of_DoblDobl_Sample_Lists )
                            return double_double is
  begin
    return Minimal_Distance(Distances(grid));
  end Minimal_Distance;

-- SELECTORS :

  function Abscisses
             ( grid : Array_of_DoblDobl_Sample_Lists; i : natural32 )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(0..integer32(i));

    use DoblDobl_Complex_Numbers;

  begin
    for j in 0..integer32(i) loop
      declare
        hyp : constant DoblDobl_Complex_VecVecs.VecVec
            := Hyperplane_Sections(Head_Of(grid(j)));
      begin
        res(j) := -hyp(hyp'first)(0);
      end;
    end loop;
    return res;
  end Abscisses;

  function Extract_Samples
                ( grid : Array_of_DoblDobl_Sample_Lists )
                return DoblDobl_Complex_VecVecs.VecVec is

    d : constant integer32 := grid'last;
    res : DoblDobl_Complex_VecVecs.VecVec(1..d*(d+1));
    tmp : DoblDobl_Sample_List;
    spt : DoblDobl_Sample;
    ind : integer32 := 1;
    rpx : DoblDobl_Complex_Vectors.Vector(1..2);

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Vectors;

  begin
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        spt := Head_Of(tmp);
        rpx(1) := Sample_Point(spt).v(1);
        rpx(2) := Sample_Point(spt).v(2); 
        res(ind) := new DoblDobl_Complex_Vectors.Vector'(rpx);
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

  function Rotate_and_Project ( v,x : DoblDobl_Complex_Vectors.Vector )
                              return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(1..2);

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Vectors;

  begin
    res(1) := v(x'range)*x;
    res(2) := v(2)*x(1) - v(1)*x(2);
    return res;
  end Rotate_and_Project;

  function Rotate_and_Project2
               ( v,x : DoblDobl_Complex_Vectors.Vector )
               return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

  begin
    return v(2)*x(1) - v(1)*x(2);
  end Rotate_and_Project2;

  function Inverse_Rotate ( v,z : DoblDobl_Complex_Vectors.Vector )
                          return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Vectors;

    res : Vector(1..2);
    det : Complex_Number := v(1)*v(1) + v(2)*v(2);

  begin
    res(1) := (z(2)*v(2) + z(1)*v(1))/det;
    res(2) := (v(2)*z(1) - v(1)*z(2))/det;
    return res;
  end Inverse_Rotate;

  function Rotate_Samples ( d,k : natural32;
                            rot : DoblDobl_Complex_Vectors.Vector;
                            grid : Array_of_DoblDobl_Sample_Lists )
                          return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(1..integer32(d),0..integer32(k));
    tmp : DoblDobl_Sample_List;
    zero : constant double_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DoblDobl_Complex_Numbers.Create(zero);
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

end DoblDobl_Rectangular_Sample_Grids;
