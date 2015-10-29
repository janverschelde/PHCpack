with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Parameter_Systems;

package body Parameter_Homotopy_State is

-- INTERNAL DATA :

  nbequ,nbvar,nbpar : integer32;
  paridx : Standard_Integer_Vectors.Link_to_Vector;

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

-- DESTRUCTOR :

  procedure Clear is

    use Standard_Integer_Vectors;

  begin
    nbequ := 0;
    nbvar := 0;
    nbpar := 0;
    if paridx /= null
     then Standard_Integer_Vectors.Clear(paridx);
    end if;
  end Clear;

begin
  nbequ := 0;
  nbvar := 0;
  nbpar := 0;
end Parameter_Homotopy_State;
