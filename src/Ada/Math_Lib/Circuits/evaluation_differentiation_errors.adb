with Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;

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

  function Difference ( s : TripDobl_Complex_Series.Link_to_Series;
                        c : TripDobl_Complex_Vectors.Link_to_Vector )
                      return triple_double is

    use Standard_Integer_Numbers;
    use TripDobl_Complex_Numbers;
    use TripDobl_Complex_Vectors;
    use TripDobl_Complex_Series;

    res : triple_double := create(0.0);
    avl : triple_double;
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

  function Difference ( s : PentDobl_Complex_Series.Link_to_Series;
                        c : PentDobl_Complex_Vectors.Link_to_Vector )
                      return penta_double is

    use Standard_Integer_Numbers;
    use PentDobl_Complex_Numbers;
    use PentDobl_Complex_Vectors;
    use PentDobl_Complex_Series;

    res : penta_double := create(0.0);
    avl : penta_double;
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

  function Difference ( s : OctoDobl_Complex_Series.Link_to_Series;
                        c : OctoDobl_Complex_Vectors.Link_to_Vector )
                      return octo_double is

    use Standard_Integer_Numbers;
    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Vectors;
    use OctoDobl_Complex_Series;

    res : octo_double := create(0.0);
    avl : octo_double;
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

  function Difference ( s : DecaDobl_Complex_Series.Link_to_Series;
                        c : DecaDobl_Complex_Vectors.Link_to_Vector )
                      return deca_double is

    use Standard_Integer_Numbers;
    use DecaDobl_Complex_Numbers;
    use DecaDobl_Complex_Vectors;
    use DecaDobl_Complex_Series;

    res : deca_double := create(0.0);
    avl : deca_double;
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

  function Difference ( x : TripDobl_Complex_Vectors.Link_to_Vector;
                        y : TripDobl_Complex_Vectors.Link_to_Vector )
                      return triple_double is

    use TripDobl_Complex_Numbers;

    res : triple_double := create(0.0);
    avl : triple_double;
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

  function Difference ( x : PentDobl_Complex_Vectors.Link_to_Vector;
                        y : PentDobl_Complex_Vectors.Link_to_Vector )
                      return penta_double is

    use PentDobl_Complex_Numbers;

    res : penta_double := create(0.0);
    avl : penta_double;
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

  function Difference ( x : OctoDobl_Complex_Vectors.Link_to_Vector;
                        y : OctoDobl_Complex_Vectors.Link_to_Vector )
                      return octo_double is

    use OctoDobl_Complex_Numbers;

    res : octo_double := create(0.0);
    avl : octo_double;
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

  function Difference ( x : DecaDobl_Complex_Vectors.Link_to_Vector;
                        y : DecaDobl_Complex_Vectors.Link_to_Vector )
                      return deca_double is

    use DecaDobl_Complex_Numbers;

    res : deca_double := create(0.0);
    avl : deca_double;
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

  function Difference ( s : TripDobl_Complex_Series_Vectors.Vector;
                        c : TripDobl_Complex_VecVecs.VecVec )
                      return triple_double is

    res : triple_double := create(0.0);
    val : triple_double;

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

  function Difference ( s : PentDobl_Complex_Series_Vectors.Vector;
                        c : PentDobl_Complex_VecVecs.VecVec )
                      return penta_double is

    res : penta_double := create(0.0);
    val : penta_double;

  begin
    for i in s'range loop
      val := Difference(s(i),c(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( s : OctoDobl_Complex_Series_Vectors.Vector;
                        c : OctoDobl_Complex_VecVecs.VecVec )
                      return octo_double is

    res : octo_double := create(0.0);
    val : octo_double;

  begin
    for i in s'range loop
      val := Difference(s(i),c(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( s : DecaDobl_Complex_Series_Vectors.Vector;
                        c : DecaDobl_Complex_VecVecs.VecVec )
                      return deca_double is

    res : deca_double := create(0.0);
    val : deca_double;

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

  function Difference ( x : TripDobl_Complex_VecVecs.VecVec;
                        y : TripDobl_Complex_VecVecs.VecVec )
                      return triple_double is

    res : triple_double := create(0.0);
    val : triple_double;

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

  function Difference ( x : PentDobl_Complex_VecVecs.VecVec;
                        y : PentDobl_Complex_VecVecs.VecVec )
                      return penta_double is

    res : penta_double := create(0.0);
    val : penta_double;

  begin
    for i in x'range loop
      val := Difference(x(i),y(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( x : OctoDobl_Complex_VecVecs.VecVec;
                        y : OctoDobl_Complex_VecVecs.VecVec )
                      return octo_double is

    res : octo_double := create(0.0);
    val : octo_double;

  begin
    for i in x'range loop
      val := Difference(x(i),y(i));
      res := res + val;
    end loop;
    return res;
  end Difference;

  function Difference ( x : DecaDobl_Complex_VecVecs.VecVec;
                        y : DecaDobl_Complex_VecVecs.VecVec )
                      return deca_double is

    res : deca_double := create(0.0);
    val : deca_double;

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

  function Difference ( jm : TripDobl_Complex_Series_Matrices.Matrix;
                        vm : TripDobl_Complex_VecMats.VecMat )
                      return triple_double is

    use TripDobl_Complex_Numbers;

    res : triple_double := create(0.0);
    avl : triple_double;
    val : Complex_Number;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        for k in jm(i,j).cff'range loop
          declare
            mat : constant TripDobl_Complex_Matrices.Link_to_Matrix := vm(k);
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

  function Difference ( jm : PentDobl_Complex_Series_Matrices.Matrix;
                        vm : PentDobl_Complex_VecMats.VecMat )
                      return penta_double is

    use PentDobl_Complex_Numbers;

    res : penta_double := create(0.0);
    avl : penta_double;
    val : Complex_Number;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        for k in jm(i,j).cff'range loop
          declare
            mat : constant PentDobl_Complex_Matrices.Link_to_Matrix := vm(k);
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

  function Difference ( jm : OctoDobl_Complex_Series_Matrices.Matrix;
                        vm : OctoDobl_Complex_VecMats.VecMat )
                      return octo_double is

    use OctoDobl_Complex_Numbers;

    res : octo_double := create(0.0);
    avl : octo_double;
    val : Complex_Number;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        for k in jm(i,j).cff'range loop
          declare
            mat : constant OctoDobl_Complex_Matrices.Link_to_Matrix := vm(k);
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

  function Difference ( jm : DecaDobl_Complex_Series_Matrices.Matrix;
                        vm : DecaDobl_Complex_VecMats.VecMat )
                      return deca_double is

    use DecaDobl_Complex_Numbers;

    res : deca_double := create(0.0);
    avl : deca_double;
    val : Complex_Number;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        for k in jm(i,j).cff'range loop
          declare
            mat : constant DecaDobl_Complex_Matrices.Link_to_Matrix := vm(k);
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

  function Difference ( vm1 : TripDobl_Complex_VecMats.VecMat;
                        vm2 : TripDobl_Complex_VecMats.VecMat )
                      return triple_double is

    use TripDobl_Complex_Numbers;

    res : triple_double := create(0.0);
    avl : triple_double;
    val : Complex_Number;

  begin
    for k in vm1'range loop
      declare
        mat1 : constant TripDobl_Complex_Matrices.Link_to_Matrix := vm1(k);
        mat2 : constant TripDobl_Complex_Matrices.Link_to_Matrix := vm2(k);
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

  function Difference ( vm1 : PentDobl_Complex_VecMats.VecMat;
                        vm2 : PentDobl_Complex_VecMats.VecMat )
                      return penta_double is

    use PentDobl_Complex_Numbers;

    res : penta_double := create(0.0);
    avl : penta_double;
    val : Complex_Number;

  begin
    for k in vm1'range loop
      declare
        mat1 : constant PentDobl_Complex_Matrices.Link_to_Matrix := vm1(k);
        mat2 : constant PentDobl_Complex_Matrices.Link_to_Matrix := vm2(k);
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

  function Difference ( vm1 : OctoDobl_Complex_VecMats.VecMat;
                        vm2 : OctoDobl_Complex_VecMats.VecMat )
                      return octo_double is

    use OctoDobl_Complex_Numbers;

    res : octo_double := create(0.0);
    avl : octo_double;
    val : Complex_Number;

  begin
    for k in vm1'range loop
      declare
        mat1 : constant OctoDobl_Complex_Matrices.Link_to_Matrix := vm1(k);
        mat2 : constant OctoDobl_Complex_Matrices.Link_to_Matrix := vm2(k);
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

  function Difference ( vm1 : DecaDobl_Complex_VecMats.VecMat;
                        vm2 : DecaDobl_Complex_VecMats.VecMat )
                      return deca_double is

    use DecaDobl_Complex_Numbers;

    res : deca_double := create(0.0);
    avl : deca_double;
    val : Complex_Number;

  begin
    for k in vm1'range loop
      declare
        mat1 : constant DecaDobl_Complex_Matrices.Link_to_Matrix := vm1(k);
        mat2 : constant DecaDobl_Complex_Matrices.Link_to_Matrix := vm2(k);
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
             ( x,y : in TripDobl_Complex_Vectors.Vector )
             return triple_double is

    use TripDobl_Complex_numbers;

    res : triple_double := create(0.0);
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
             ( x,y : in PentDobl_Complex_Vectors.Vector )
             return penta_double is

    use PentDobl_Complex_numbers;

    res : penta_double := create(0.0);
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      res := res + AbsVal(val);
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( x,y : in OctoDobl_Complex_Vectors.Vector )
             return octo_double is

    use OctoDobl_Complex_numbers;

    res : octo_double := create(0.0);
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      res := res + AbsVal(val);
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( x,y : in DecaDobl_Complex_Vectors.Vector )
             return deca_double is

    use DecaDobl_Complex_numbers;

    res : deca_double := create(0.0);
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
             ( A,B : in TripDobl_Complex_Matrices.Matrix )
             return triple_double is

    use TripDobl_Complex_numbers;

    res : triple_double := create(0.0);
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

  function Sum_of_Errors
             ( A,B : in PentDobl_Complex_Matrices.Matrix )
             return penta_double is

    use PentDobl_Complex_numbers;

    res : penta_double := create(0.0);
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
             ( A,B : in OctoDobl_Complex_Matrices.Matrix )
             return octo_double is

    use OctoDobl_Complex_numbers;

    res : octo_double := create(0.0);
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
             ( A,B : in DecaDobl_Complex_Matrices.Matrix )
             return deca_double is

    use DecaDobl_Complex_numbers;

    res : deca_double := create(0.0);
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
