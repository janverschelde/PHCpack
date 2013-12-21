package body Interpolation_Point_Lists is

-- CREATORS :

  procedure Create ( spl : in Standard_Sample_List;
                     prospl : in Standard_Complex_VecVecs.VecVec;
                     first,last : in out Standard_Sample_Node_List ) is

    tmp : Standard_Sample_List := spl;

  begin
    for i in prospl'range loop
      Append(first,last,Create(Head_Of(tmp),prospl(i).all));
      tmp := Tail_Of(tmp);
    end loop;
  end Create;

  procedure Create ( spl : in Multprec_Sample_List;
                     prospl : in Multprec_Complex_VecVecs.VecVec;
                     first,last : in out Multprec_Sample_Node_List ) is

    tmp : Multprec_Sample_List := spl;

  begin
    for i in prospl'range loop
      Append(first,last,Create(Head_Of(tmp),prospl(i).all));
      tmp := Tail_Of(tmp);
    end loop;
  end Create;

-- SELECTORS :

-- DESTRUCTORS :

  procedure Shallow_Clear ( l : in out Standard_Sample_Node_List ) is

    tmp : Standard_Sample_Node_List := l;

  begin
    while not Is_Null(tmp) loop
      declare
        sn : Standard_Sample_Node := Head_Of(tmp);
      begin
        Shallow_Clear(sn);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Standard_Sample_Nodes.Clear
      (Lists_of_Standard_Sample_Nodes.List(l));
  end Shallow_Clear;

  procedure Shallow_Clear ( l : in out Multprec_Sample_Node_List ) is

    tmp : Multprec_Sample_Node_List := l;

  begin
    while not Is_Null(tmp) loop
      declare
        sn : Multprec_Sample_Node := Head_Of(tmp);
      begin
        Shallow_Clear(sn);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Multprec_Sample_Nodes.Clear
      (Lists_of_Multprec_Sample_Nodes.List(l));
  end Shallow_Clear;

  procedure Deep_Clear ( l : in out Standard_Sample_Node_List ) is

    tmp : Standard_Sample_Node_List := l;

  begin
    while not Is_Null(tmp) loop
      declare
        sn : Standard_Sample_Node := Head_Of(tmp);
      begin
        Deep_Clear(sn);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Standard_Sample_Nodes.Clear
      (Lists_of_Standard_Sample_Nodes.List(l));
  end Deep_Clear;

  procedure Deep_Clear ( l : in out Multprec_Sample_Node_List ) is

    tmp : Multprec_Sample_Node_List := l;

  begin
    while not Is_Null(tmp) loop
      declare
        sn : Multprec_Sample_Node := Head_Of(tmp);
      begin
        Deep_Clear(sn);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Multprec_Sample_Nodes.Clear
      (Lists_of_Multprec_Sample_Nodes.List(l));
  end Deep_Clear;

end Interpolation_Point_Lists;
