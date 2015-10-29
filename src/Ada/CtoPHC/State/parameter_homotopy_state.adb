with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Parameter_Systems;

package body Parameter_Homotopy_State is

-- INTERNAL DATA :

  nbequ,nbvar,nbpar : integer32;
  paridx : Standard_Integer_Vectors.Link_to_Vector;
  st_startv : Standard_Complex_Vectors.Link_to_Vector;
  dd_startv : DoblDobl_Complex_Vectors.Link_to_Vector;
  qd_startv : QuadDobl_Complex_Vectors.Link_to_Vector;
  st_target : Standard_Complex_Vectors.Link_to_Vector;
  dd_target : DoblDobl_Complex_Vectors.Link_to_Vector;
  qd_target : QuadDobl_Complex_Vectors.Link_to_Vector;

-- CONSTRUCTORS :

  procedure Set_Number_of_Equations ( n : in integer32 ) is
  begin
    nbequ := n;
  end Set_Number_of_Equations;

  procedure Set_Number_of_Variables ( n : in integer32 ) is
  begin
    nbvar := n;
  end Set_Number_of_Variables;

  procedure Set_Number_of_Parameters ( n : in integer32 ) is
  begin
    nbpar := n;
  end Set_Number_of_Parameters;

  procedure Set_Indices ( idx : in Standard_Integer_Vectors.Vector ) is
  begin
    paridx := new Standard_Integer_Vectors.Vector'(idx);
    Standard_Parameter_Systems.Sort(paridx.all);
  end Set_Indices;

  procedure Set_Start ( v : in Standard_Complex_Vectors.Vector ) is
  begin
    st_startv := new Standard_Complex_Vectors.Vector'(v);
  end Set_Start;

  procedure Set_Start ( v : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    dd_startv := new DoblDobl_Complex_Vectors.Vector'(v);
  end Set_Start;

  procedure Set_Start ( v : in QuadDobl_Complex_Vectors.Vector ) is
  begin
    qd_startv := new QuadDobl_Complex_Vectors.Vector'(v);
  end Set_Start;

  procedure Set_Target ( v : in Standard_Complex_Vectors.Vector ) is
  begin
    st_target := new Standard_Complex_Vectors.Vector'(v);
  end Set_Target;

  procedure Set_Target ( v : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    dd_target := new DoblDobl_Complex_Vectors.Vector'(v);
  end Set_Target;

  procedure Set_Target ( v : in QuadDobl_Complex_Vectors.Vector ) is
  begin
    qd_target := new QuadDobl_Complex_Vectors.Vector'(v);
  end Set_Target;

-- SELECTORS :

  function Get_Number_of_Equations return integer32 is
  begin
    return nbequ;
  end Get_Number_of_Equations;

  function Get_Number_of_Variables return integer32 is
  begin
    return nbvar;
  end Get_Number_of_Variables;

  function Get_Number_of_Parameters return integer32 is
  begin
    return nbpar;
  end Get_Number_of_Parameters;

  function Get_Indices return Standard_Integer_Vectors.Link_to_Vector is
  begin
    return paridx;
  end Get_Indices;

  function Get_Start return Standard_Complex_Vectors.Link_to_Vector is
  begin
    return st_startv;
  end Get_Start;

  function Get_Start return DoblDobl_Complex_Vectors.Link_to_Vector is
  begin
    return dd_startv;
  end Get_Start;

  function Get_Start return QuadDobl_Complex_Vectors.Link_to_Vector is
  begin
    return qd_startv;
  end Get_Start;

  function Get_Target return Standard_Complex_Vectors.Link_to_Vector is
  begin
    return st_target;
  end Get_Target;

  function Get_Target return DoblDobl_Complex_Vectors.Link_to_Vector is
  begin
    return dd_target;
  end Get_Target;

  function Get_Target return QuadDobl_Complex_Vectors.Link_to_Vector is
  begin
    return qd_target;
  end Get_Target;

-- DESTRUCTOR :

  procedure Clear is

    use Standard_Integer_Vectors;
    use Standard_Complex_Vectors;
    use DoblDobl_Complex_Vectors;
    use QuadDobl_Complex_Vectors;

  begin
    nbequ := 0;
    nbvar := 0;
    nbpar := 0;
    if paridx /= null
     then Standard_Integer_Vectors.Clear(paridx);
    end if;
    if st_startv /= null
     then Standard_Complex_Vectors.Clear(st_startv);
    end if;
    if dd_startv /= null
     then DoblDobl_Complex_Vectors.Clear(dd_startv);
    end if;
    if qd_startv /= null
     then QuadDobl_Complex_Vectors.Clear(qd_startv);
    end if;
    if st_target /= null
     then Standard_Complex_Vectors.Clear(st_target);
    end if;
    if dd_target /= null
     then DoblDobl_Complex_Vectors.Clear(dd_target);
    end if;
    if qd_target /= null
     then QuadDobl_Complex_Vectors.Clear(qd_target);
    end if;
  end Clear;

begin
  nbequ := 0;
  nbvar := 0;
  nbpar := 0;
end Parameter_Homotopy_State;
