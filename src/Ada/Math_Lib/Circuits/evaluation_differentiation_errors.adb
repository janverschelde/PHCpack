with Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Matrices;

package body Evaluation_Differentiation_Errors is

  function Difference ( s : Standard_Complex_Series.Link_to_Series;
                        c : Standard_Complex_Vectors.Link_to_Vector )
                      return double_float is

    use Standard_Integer_Numbers;
    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_Series;

    res : double_float := 0.0;
    avl : double_float;
    val : Complex_Number;

  begin
    if s /= null and c /= null then
      for i in c'range loop
        exit when (i > s.cff'last);
        val := s.cff(i) - c(i);
        avl := AbsVal(val);
        res := res + avl;
      end loop;
    end if;
    return res;
  end Difference;

  function Difference ( s : DoblDobl_Complex_Series.Link_to_Series;
                        c : DoblDobl_Complex_Vectors.Link_to_Vector )
                      return double_double is

    use Standard_Integer_Numbers;
    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Series;

    res : double_double := create(0.0);
    avl : double_double;
    val : Complex_Number;

  begin
    if s /= null and c /= null then
      for i in c'range loop
        exit when (i > s.cff'last);
        val := s.cff(i) - c(i);
        avl := AbsVal(val);
        res := res + avl;
      end loop;
    end if;
    return res;
  end Difference;

  function Difference ( s : QuadDobl_Complex_Series.Link_to_Series;
                        c : QuadDobl_Complex_Vectors.Link_to_Vector )
                      return quad_double is

    use Standard_Integer_Numbers;
    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Series;

    res : quad_double := create(0.0);
    avl : quad_double;
    val : Complex_Number;

  begin
    if s /= null and c /= null then
      for i in c'range loop
        exit when (i > s.cff'last);
        val := s.cff(i) - c(i);
        avl := AbsVal(val);
        res := res + avl;
      end loop;
    end if;
    return res;
  end Difference;

  function Difference ( x : Standard_Complex_Vectors.Link_to_Vector;
                        y : Standard_Complex_Vectors.Link_to_Vector )
                      return double_float is

    use Standard_Complex_Numbers;

    res : double_float := 0.0;
    avl : double_float;
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

  function Difference ( x : DoblDobl_Complex_Vectors.Link_to_Vector;
                        y : DoblDobl_Complex_Vectors.Link_to_Vector )
                      return double_double is

    use DoblDobl_Complex_Numbers;

    res : double_double := create(0.0);
    avl : double_double;
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

  function Difference ( x : QuadDobl_Complex_Vectors.Link_to_Vector;
                        y : QuadDobl_Complex_Vectors.Link_to_Vector )
                      return quad_double is

    use QuadDobl_Complex_Numbers;

    res : quad_double := create(0.0);
    avl : quad_double;
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

  function Difference ( s : Standard_Complex_Series_Vectors.Vector;
                        c : Standard_Complex_VecVecs.VecVec )
                      return double_float is

    res : double_float := 0.0;
    val : double_float;

  begin
    for i in s'range loop
      val := Difference(s(i),c(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( s : DoblDobl_Complex_Series_Vectors.Vector;
                        c : DoblDobl_Complex_VecVecs.VecVec )
                      return double_double is

    res : double_double := create(0.0);
    val : double_double;

  begin
    for i in s'range loop
      val := Difference(s(i),c(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( s : QuadDobl_Complex_Series_Vectors.Vector;
                        c : QuadDobl_Complex_VecVecs.VecVec )
                      return quad_double is

    res : quad_double := create(0.0);
    val : quad_double;

  begin
    for i in s'range loop
      val := Difference(s(i),c(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( x : Standard_Complex_VecVecs.VecVec;
                        y : Standard_Complex_VecVecs.VecVec )
                      return double_float is

    res : double_float := 0.0;
    val : double_float;

  begin
    for i in x'range loop
      val := Difference(x(i),y(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( x : DoblDobl_Complex_VecVecs.VecVec;
                        y : DoblDobl_Complex_VecVecs.VecVec )
                      return double_double is

    res : double_double := create(0.0);
    val : double_double;

  begin
    for i in x'range loop
      val := Difference(x(i),y(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( x : QuadDobl_Complex_VecVecs.VecVec;
                        y : QuadDobl_Complex_VecVecs.VecVec )
                      return quad_double is

    res : quad_double := create(0.0);
    val : quad_double;

  begin
    for i in x'range loop
      val := Difference(x(i),y(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( jm : Standard_Complex_Series_Matrices.Matrix;
                        vm : Standard_Complex_VecMats.VecMat )
                      return double_float is

    use Standard_Complex_Numbers;

    res : double_float := 0.0;
    avl : double_float;
    val : Complex_Number;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        for k in jm(i,j).cff'range loop
          declare
            mat : constant Standard_Complex_Matrices.Link_to_Matrix := vm(k);
          begin
            val := jm(i,j).cff(k) - mat(i,j);
            avl := AbsVal(val);
            res := res + avl;
          end;
        end loop;
      end loop;
    end loop;
    return res;
  end Difference;

  function Difference ( jm : DoblDobl_Complex_Series_Matrices.Matrix;
                        vm : DoblDobl_Complex_VecMats.VecMat )
                      return double_double is

    use DoblDobl_Complex_Numbers;

    res : double_double := create(0.0);
    avl : double_double;
    val : Complex_Number;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        for k in jm(i,j).cff'range loop
          declare
            mat : constant DoblDobl_Complex_Matrices.Link_to_Matrix := vm(k);
          begin
            val := jm(i,j).cff(k) - mat(i,j);
            avl := AbsVal(val);
            res := res + avl;
          end;
        end loop;
      end loop;
    end loop;
    return res;
  end Difference;

  function Difference ( jm : QuadDobl_Complex_Series_Matrices.Matrix;
                        vm : QuadDobl_Complex_VecMats.VecMat )
                      return quad_double is

    use QuadDobl_Complex_Numbers;

    res : quad_double := create(0.0);
    avl : quad_double;
    val : Complex_Number;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        for k in jm(i,j).cff'range loop
          declare
            mat : constant QuadDobl_Complex_Matrices.Link_to_Matrix := vm(k);
          begin
            val := jm(i,j).cff(k) - mat(i,j);
            avl := AbsVal(val);
            res := res + avl;
          end;
        end loop;
      end loop;
    end loop;
    return res;
  end Difference;

  function Difference ( vm1 : Standard_Complex_VecMats.VecMat;
                        vm2 : Standard_Complex_VecMats.VecMat )
                      return double_float is

    use Standard_Complex_Numbers;

    res : double_float := 0.0;
    avl : double_float;
    val : Complex_Number;

  begin
    for k in vm1'range loop
      declare
        mat1 : constant Standard_Complex_Matrices.Link_to_Matrix := vm1(k);
        mat2 : constant Standard_Complex_Matrices.Link_to_Matrix := vm2(k);
      begin
        for i in mat1'range(1) loop
          for j in mat2'range(2) loop
            val := mat1(i,j) - mat2(i,j);
            avl := AbsVal(val);
            res := res + avl;
          end loop;
        end loop;
      end;
    end loop;
    return res;
  end Difference;

  function Difference ( vm1 : DoblDobl_Complex_VecMats.VecMat;
                        vm2 : DoblDobl_Complex_VecMats.VecMat )
                      return double_double is

    use DoblDobl_Complex_Numbers;

    res : double_double := create(0.0);
    avl : double_double;
    val : Complex_Number;

  begin
    for k in vm1'range loop
      declare
        mat1 : constant DoblDobl_Complex_Matrices.Link_to_Matrix := vm1(k);
        mat2 : constant DoblDobl_Complex_Matrices.Link_to_Matrix := vm2(k);
      begin
        for i in mat1'range(1) loop
          for j in mat2'range(2) loop
            val := mat1(i,j) - mat2(i,j);
            avl := AbsVal(val);
            res := res + avl;
          end loop;
        end loop;
      end;
    end loop;
    return res;
  end Difference;

  function Difference ( vm1 : QuadDobl_Complex_VecMats.VecMat;
                        vm2 : QuadDobl_Complex_VecMats.VecMat )
                      return quad_double is

    use QuadDobl_Complex_Numbers;

    res : quad_double := create(0.0);
    avl : quad_double;
    val : Complex_Number;

  begin
    for k in vm1'range loop
      declare
        mat1 : constant QuadDobl_Complex_Matrices.Link_to_Matrix := vm1(k);
        mat2 : constant QuadDobl_Complex_Matrices.Link_to_Matrix := vm2(k);
      begin
        for i in mat1'range(1) loop
          for j in mat2'range(2) loop
            val := mat1(i,j) - mat2(i,j);
            avl := AbsVal(val);
            res := res + avl;
          end loop;
        end loop;
      end;
    end loop;
    return res;
  end Difference;

  function Sum_of_Errors
             ( x,y : in Standard_Complex_Vectors.Vector )
             return double_float is

    use Standard_Complex_numbers;

    res : double_float := 0.0;
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      res := res + AbsVal(val);
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( x,y : in DoblDobl_Complex_Vectors.Vector )
             return double_double is

    use DoblDobl_Complex_numbers;

    res : double_double := create(0.0);
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      res := res + AbsVal(val);
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( x,y : in QuadDobl_Complex_Vectors.Vector )
             return quad_double is

    use QuadDobl_Complex_numbers;

    res : quad_double := create(0.0);
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      res := res + AbsVal(val);
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( A,B : in Standard_Complex_Matrices.Matrix )
             return double_float is

    use Standard_Complex_numbers;

    res : double_float := 0.0;
    val : Complex_Number;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        val := A(i,j) - B(i,j);
        res := res + AbsVal(val);
      end loop;
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( A,B : in DoblDobl_Complex_Matrices.Matrix )
             return double_double is

    use DoblDobl_Complex_numbers;

    res : double_double := create(0.0);
    val : Complex_Number;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        val := A(i,j) - B(i,j);
        res := res + AbsVal(val);
      end loop;
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( A,B : in QuadDobl_Complex_Matrices.Matrix )
             return quad_double is

    use QuadDobl_Complex_numbers;

    res : quad_double := create(0.0);
    val : Complex_Number;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        val := A(i,j) - B(i,j);
        res := res + AbsVal(val);
      end loop;
    end loop;
    return res;
  end Sum_of_Errors;

end Evaluation_Differentiation_Errors;
