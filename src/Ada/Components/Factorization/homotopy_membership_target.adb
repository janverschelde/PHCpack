with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with Planes_and_Polynomials;
with Witness_Sets;

package body Homotopy_Membership_Target is

  function Adjusted_Slices 
             ( sli : Standard_Complex_VecVecs.VecVec;
               sol : Standard_Complex_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(sli'range);

    use Standard_Complex_Numbers;

  begin
    for i in res'range loop
      res(i) := new Standard_Complex_Vectors.Vector'(sli(i).all);
      res(i)(0) := -res(i)(1)*sol(1);
      for j in 2..sol'last loop
        res(i)(0) := res(i)(0) - res(i)(j)*sol(j);
      end loop;
    end loop;
    return res;
  end Adjusted_Slices;

  function Adjusted_Slices 
             ( sli : DoblDobl_Complex_VecVecs.VecVec;
               sol : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(sli'range);

    use DoblDobl_Complex_Numbers;

  begin
    for i in res'range loop
      res(i) := new DoblDobl_Complex_Vectors.Vector'(sli(i).all);
      res(i)(0) := -res(i)(1)*sol(1);
      for j in 2..sol'last loop
        res(i)(0) := res(i)(0) - res(i)(j)*sol(j);
      end loop;
    end loop;
    return res;
  end Adjusted_Slices;

  function Adjusted_Slices 
             ( sli : QuadDobl_Complex_VecVecs.VecVec;
               sol : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(sli'range);

    use QuadDobl_Complex_Numbers;

  begin
    for i in res'range loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'(sli(i).all);
      res(i)(0) := -res(i)(1)*sol(1);
      for j in 2..sol'last loop
        res(i)(0) := res(i)(0) - res(i)(j)*sol(j);
      end loop;
    end loop;
    return res;
  end Adjusted_Slices;

  function Adjusted_Target
             ( ep : Standard_Complex_Poly_Systems.Poly_Sys;
               dim : natural32;
               pnt : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    idm : constant integer32 := integer32(dim);
    res : Standard_Complex_Poly_Systems.Poly_Sys(ep'range);
    sli : Standard_Complex_VecVecs.VecVec(1..idm)
        := Witness_Sets.Slices(ep,dim);
    adj : Standard_Complex_VecVecs.VecVec(1..idm)
        := Adjusted_Slices(sli,pnt);
    offset : constant integer32 := ep'last-idm;

  begin
    for i in ep'first..offset loop
      Standard_Complex_Polynomials.Copy(ep(i),res(i));
    end loop;
    for i in 1..idm loop
      res(offset+i) := Planes_and_Polynomials.Hyperplane(adj(i).all);
    end loop;
    Standard_Complex_VecVecs.Clear(sli);
    Standard_Complex_VecVecs.Clear(adj);
    return res;
  end Adjusted_Target;

  function Adjusted_Target
             ( ep : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               dim : natural32;
               pnt : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    idm : constant integer32 := integer32(dim);
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(ep'range);
    sli : DoblDobl_Complex_VecVecs.VecVec(1..idm)
        := Witness_Sets.Slices(ep,dim);
    adj : DoblDobl_Complex_VecVecs.VecVec(1..idm)
        := Adjusted_Slices(sli,pnt);
    offset : constant integer32 := ep'last-idm;

  begin
    for i in ep'first..offset loop
      DoblDobl_Complex_Polynomials.Copy(ep(i),res(i));
    end loop;
    for i in 1..idm loop
      res(offset+i) := Planes_and_Polynomials.Hyperplane(adj(i).all);
    end loop;
    DoblDobl_Complex_VecVecs.Clear(sli);
    DoblDobl_Complex_VecVecs.Clear(adj);
    return res;
  end Adjusted_Target;

  function Adjusted_Target
             ( ep : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               dim : natural32;
               pnt : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    idm : constant integer32 := integer32(dim);
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(ep'range);
    sli : QuadDobl_Complex_VecVecs.VecVec(1..idm)
        := Witness_Sets.Slices(ep,dim);
    adj : QuadDobl_Complex_VecVecs.VecVec(1..idm)
        := Adjusted_Slices(sli,pnt);
    offset : constant integer32 := ep'last-idm;

  begin
    for i in ep'first..offset loop
      QuadDobl_Complex_Polynomials.Copy(ep(i),res(i));
    end loop;
    for i in 1..idm loop
      res(offset+i) := Planes_and_Polynomials.Hyperplane(adj(i).all);
    end loop;
    QuadDobl_Complex_VecVecs.Clear(sli);
    QuadDobl_Complex_VecVecs.Clear(adj);
    return res;
  end Adjusted_Target;

  function Adjusted_Target
             ( ep : Standard_Complex_Laur_Systems.Laur_Sys;
               dim : natural32;
               pnt : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    idm : constant integer32 := integer32(dim);
    res : Standard_Complex_Laur_Systems.Laur_Sys(ep'range);
    sli : Standard_Complex_VecVecs.VecVec(1..idm)
        := Witness_Sets.Slices(ep,dim);
    adj : Standard_Complex_VecVecs.VecVec(1..idm)
        := Adjusted_Slices(sli,pnt);
    offset : constant integer32 := ep'last-idm;

  begin
    for i in ep'first..offset loop
      Standard_Complex_Laurentials.Copy(ep(i),res(i));
    end loop;
    for i in 1..idm loop
      res(offset+i) := Planes_and_Polynomials.Hyperplane(adj(i).all);
    end loop;
    Standard_Complex_VecVecs.Clear(sli);
    Standard_Complex_VecVecs.Clear(adj);
    return res;
  end Adjusted_Target;

  function Adjusted_Target
             ( ep : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               dim : natural32;
               pnt : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    idm : constant integer32 := integer32(dim);
    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(ep'range);
    sli : DoblDobl_Complex_VecVecs.VecVec(1..idm)
        := Witness_Sets.Slices(ep,dim);
    adj : DoblDobl_Complex_VecVecs.VecVec(1..idm)
        := Adjusted_Slices(sli,pnt);
    offset : constant integer32 := ep'last-idm;

  begin
    for i in ep'first..offset loop
      DoblDobl_Complex_Laurentials.Copy(ep(i),res(i));
    end loop;
    for i in 1..idm loop
      res(offset+i) := Planes_and_Polynomials.Hyperplane(adj(i).all);
    end loop;
    DoblDobl_Complex_VecVecs.Clear(sli);
    DoblDobl_Complex_VecVecs.Clear(adj);
    return res;
  end Adjusted_Target;

  function Adjusted_Target
             ( ep : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               dim : natural32;
               pnt : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    idm : constant integer32 := integer32(dim);
    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(ep'range);
    sli : QuadDobl_Complex_VecVecs.VecVec(1..idm)
        := Witness_Sets.Slices(ep,dim);
    adj : QuadDobl_Complex_VecVecs.VecVec(1..idm)
        := Adjusted_Slices(sli,pnt);
    offset : constant integer32 := ep'last-idm;

  begin
    for i in ep'first..offset loop
      QuadDobl_Complex_Laurentials.Copy(ep(i),res(i));
    end loop;
    for i in 1..idm loop
      res(offset+i) := Planes_and_Polynomials.Hyperplane(adj(i).all);
    end loop;
    QuadDobl_Complex_VecVecs.Clear(sli);
    QuadDobl_Complex_VecVecs.Clear(adj);
    return res;
  end Adjusted_Target;

end Homotopy_Membership_Target;
