with Standard_Integer_Vectors;

package body Transformation_of_Supports is

  function Transform ( support : Lists_of_Integer_Vectors.List;
                       transfo : Standard_Integer_Matrices.Matrix )
                     return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Matrices,Lists_of_Integer_Vectors;

    res,res_last : List;
    tmp : List := support;
    lv : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      Append(res,res_last,transfo*lv.all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Transform;

  function Transform ( p : Standard_Complex_Laurentials.Poly;
                       transfo : Standard_Integer_Matrices.Matrix )
                     return Standard_Complex_Laurentials.Poly is

    use Standard_Integer_Matrices,Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      nt : Term;
      v : constant Standard_Integer_Vectors.Vector := transfo*t.dg.all;

    begin
      nt.cf := t.cf;
      nt.dg := new Standard_Integer_Vectors.Vector'(v);
      Add(res,nt); Clear(nt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Transform;

  function Transform ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                       transfo : Standard_Integer_Matrices.Matrix )
                     return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laur_Systems;

    res : Laur_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Transform(p(i),transfo);
    end loop;
    return res;
  end Transform;

end Transformation_of_Supports;
