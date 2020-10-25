with OctoDobl_Complex_Numbers;

package body OctoDobl_Vector_Splitters is

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return OctoDobl_Complex_Vectors.Link_to_Vector is

    zero : constant OctoDobl_Complex_Numbers.Complex_Number
         := OctoDobl_Complex_Numbers.Create(integer(0));
    cff : constant OctoDobl_Complex_Vectors.Vector(0..deg) := (0..deg => zero);
    res : constant OctoDobl_Complex_Vectors.Link_to_Vector
        := new OctoDobl_Complex_Vectors.Vector'(cff);

  begin
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return OctoDobl_Complex_VecVecs.VecVec is

    res : OctoDobl_Complex_VecVecs.VecVec(1..dim);

  begin
    for k in res'range loop
      res(k) := Allocate_Complex_Coefficients(deg);
    end loop;
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return OctoDobl_Complex_VecVecs.Link_to_VecVec is

    res : OctoDobl_Complex_VecVecs.Link_to_VecVec;
    cff : constant OctoDobl_Complex_VecVecs.VecVec(1..dim)
        := Allocate_Complex_Coefficients(dim,deg);

  begin
    res := new OctoDobl_Complex_VecVecs.VecVec'(cff);
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return OctoDobl_Complex_VecVecs.VecVec is

    res : OctoDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in res'range loop
      declare
        v : constant OctoDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => OctoDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new OctoDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return OctoDobl_Complex_VecVecs.Link_to_VecVec is

    res : OctoDobl_Complex_VecVecs.Link_to_VecVec;
    vv : OctoDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in vv'range loop
      declare
        v : constant OctoDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => OctoDobl_Complex_Numbers.Create(integer(0)));
      begin
        vv(i) := new OctoDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    res := new OctoDobl_Complex_VecVecs.VecVec'(vv);
    return res;
  end Allocate;

end OctoDobl_Vector_Splitters;
