with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers_Polar;      use QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Vector_Norms;       use QuadDobl_Complex_Vector_Norms;

package body QuadDobl_Rectangular_Sample_Grids is

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean ) is
  begin
    QuadDobl_Sample_Points.Set_Polynomial_Type(laurent);
    QuadDobl_Sample_Lists.Set_Polynomial_Type(laurent);
  end Set_Polynomial_Type;

-- ROOTS OF UNITY :

  function QuadDobl_Roots_of_Unity
             ( d : natural32 ) return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with the d complex roots of unity.

    res : QuadDobl_Complex_Vectors.Vector(1..integer32(d));
    use QuadDobl_Complex_Numbers;
    qd_one : constant quad_double := create(1.0);
    one : constant Complex_Number := Create(qd_one);

  begin
    for i in 1..d loop
      res(integer32(i)) := Root(one,d,i);
    end loop;
    return res;
  end QuadDobl_Roots_of_Unity;

-- CREATORS :

  function Create1 ( sps : QuadDobl_Sample_List; m : natural32 )
                   return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..integer32(m));
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    roots : constant QuadDobl_Complex_Vectors.Vector(1..integer32(m))
          := QuadDobl_Roots_of_Unity(m);

  begin
    res(0) := sps;
    for i in 1..res'last loop
      declare
        newsli : QuadDobl_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : QuadDobl_Sample_List;
      begin
        for j in sli'range loop
          newsli(j) := new QuadDobl_Complex_Vectors.Vector'(sli(j).all);
        end loop;
        newsli(1)(0) := roots(i);
        Sample(sps,newsli,newsps,newsps_last);
        res(i) := newsps;
      end;   
    end loop;
    return res;
  end Create1;

  function Triangular_Create1 ( sps : QuadDobl_Sample_List; m : natural32 )
                              return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..integer32(m));
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Hyperplane_Sections(Head_Of(sps));
    roots : constant QuadDobl_Complex_Vectors.Vector(1..integer32(m))
          := QuadDobl_Roots_of_Unity(m);
    wrksps : QuadDobl_Sample_List;

  begin
    res(0) := sps;
    wrksps := sps;
    for i in 1..integer32(m) loop
      declare
        newsli : QuadDobl_Complex_VecVecs.VecVec(sli'range);
        newsps,newsps_last : QuadDobl_Sample_List;
      begin
        for j in sli'range loop
          newsli(j) := new QuadDobl_Complex_Vectors.Vector'(sli(j).all);
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

  function Errors ( grid : Array_of_QuadDobl_Sample_Lists ) 
                  return Quad_Double_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Quad_Double_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : QuadDobl_Sample_List;
    zero : constant quad_double := create(0.0);

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

  function Maximal_Error ( grid_errors : Quad_Double_Matrices.Matrix ) 
                         return quad_double is

    max : quad_double
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

  function Maximal_Error ( grid : Array_of_QuadDobl_Sample_Lists )
                         return quad_double is
  begin
    return Maximal_Error(Errors(grid));
  end Maximal_Error;

  function Distance ( spt1,spt2 : QuadDobl_Sample ) return quad_double is

    use QuadDobl_Complex_Vectors;

  begin
    return Max_Norm(Sample_Point(spt1).v - Sample_Point(spt2).v);
  end Distance;

  function Distance ( sps : QuadDobl_Sample_List; i : natural32;
                      spt : QuadDobl_Sample ) return quad_double is

    res,dst : quad_double;
    tmp : QuadDobl_Sample_List := sps;
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

  function Distances ( grid : Array_of_QuadDobl_Sample_Lists )
                     return Quad_Double_Matrices.Matrix is

    len : constant natural32 := Length_Of(grid(grid'first));
    res : Quad_Double_Matrices.Matrix(grid'range,1..integer32(len));
    tmp : QuadDobl_Sample_List;
    init : constant quad_double := create(1.0E+50);

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

  function Minimal_Distance ( grid_dist : Quad_Double_Matrices.Matrix )
                            return quad_double is

    min : quad_double
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

  function Minimal_Distance ( grid : Array_of_QuadDobl_Sample_Lists )
                            return quad_double is
  begin
    return Minimal_Distance(Distances(grid));
  end Minimal_Distance;

-- SELECTORS :

  function Abscisses
             ( grid : Array_of_QuadDobl_Sample_Lists; i : natural32 )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(0..integer32(i));

    use QuadDobl_Complex_Numbers;

  begin
    for j in 0..integer32(i) loop
      declare
        hyp : constant QuadDobl_Complex_VecVecs.VecVec
            := Hyperplane_Sections(Head_Of(grid(j)));
      begin
        res(j) := -hyp(hyp'first)(0);
      end;
    end loop;
    return res;
  end Abscisses;

  function Extract_Samples
                ( grid : Array_of_QuadDobl_Sample_Lists )
                return QuadDobl_Complex_VecVecs.VecVec is

    d : constant integer32 := grid'last;
    res : QuadDobl_Complex_VecVecs.VecVec(1..d*(d+1));
    tmp : QuadDobl_Sample_List;
    spt : QuadDobl_Sample;
    ind : integer32 := 1;
    rpx : QuadDobl_Complex_Vectors.Vector(1..2);

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Vectors;

  begin
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        spt := Head_Of(tmp);
        rpx(1) := Sample_Point(spt).v(1);
        rpx(2) := Sample_Point(spt).v(2); 
        res(ind) := new QuadDobl_Complex_Vectors.Vector'(rpx);
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

  function Rotate_and_Project ( v,x : QuadDobl_Complex_Vectors.Vector )
                              return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(1..2);

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Vectors;

  begin
    res(1) := v(x'range)*x;
    res(2) := v(2)*x(1) - v(1)*x(2);
    return res;
  end Rotate_and_Project;

  function Rotate_and_Project2
               ( v,x : QuadDobl_Complex_Vectors.Vector )
               return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

  begin
    return v(2)*x(1) - v(1)*x(2);
  end Rotate_and_Project2;

  function Inverse_Rotate ( v,z : QuadDobl_Complex_Vectors.Vector )
                          return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Vectors;

    res : Vector(1..2);
    det : Complex_Number := v(1)*v(1) + v(2)*v(2);

  begin
    res(1) := (z(2)*v(2) + z(1)*v(1))/det;
    res(2) := (v(2)*z(1) - v(1)*z(2))/det;
    return res;
  end Inverse_Rotate;

  function Rotate_Samples ( d,k : natural32;
                            rot : QuadDobl_Complex_Vectors.Vector;
                            grid : Array_of_QuadDobl_Sample_Lists )
                          return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(1..integer32(d),0..integer32(k));
    tmp : QuadDobl_Sample_List;
    zero : constant quad_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Complex_Numbers.Create(zero);
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

end QuadDobl_Rectangular_Sample_Grids;
