with TripDobl_Complex_Numbers;

package body TripDobl_Vector_Splitters is

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return TripDobl_Complex_Vectors.Link_to_Vector is

    zero : constant TripDobl_Complex_Numbers.Complex_Number
         := TripDobl_Complex_Numbers.Create(integer(0));
    cff : constant TripDobl_Complex_Vectors.Vector(0..deg) := (0..deg => zero);
    res : constant TripDobl_Complex_Vectors.Link_to_Vector
        := new TripDobl_Complex_Vectors.Vector'(cff);

  begin
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return TripDobl_Complex_VecVecs.VecVec is

    res : TripDobl_Complex_VecVecs.VecVec(1..dim);

  begin
    for k in res'range loop
      res(k) := Allocate_Complex_Coefficients(deg);
    end loop;
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return TripDobl_Complex_VecVecs.Link_to_VecVec is

    res : TripDobl_Complex_VecVecs.Link_to_VecVec;
    cff : constant TripDobl_Complex_VecVecs.VecVec(1..dim)
        := Allocate_Complex_Coefficients(dim,deg);

  begin
    res := new TripDobl_Complex_VecVecs.VecVec'(cff);
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return TripDobl_Complex_VecVecs.VecVec is

    res : TripDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in res'range loop
      declare
        v : constant TripDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => TripDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new TripDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return TripDobl_Complex_VecVecs.Link_to_VecVec is

    res : TripDobl_Complex_VecVecs.Link_to_VecVec;
    vv : TripDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in vv'range loop
      declare
        v : constant TripDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => TripDobl_Complex_Numbers.Create(integer(0)));
      begin
        vv(i) := new TripDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    res := new TripDobl_Complex_VecVecs.VecVec'(vv);
    return res;
  end Allocate;

end TripDobl_Vector_Splitters;
