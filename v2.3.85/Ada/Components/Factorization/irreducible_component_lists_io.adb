with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Irreducible_Components;             use Irreducible_Components;
with Irreducible_Components_io;          use Irreducible_Components_io;

package body Irreducible_Component_Lists_io is


  procedure get ( first,last : in out Standard_Irreducible_Component_List ) is
  begin
    get(Standard_Input,first,last);
  end get;

  procedure get ( file : in file_type;
                  first,last : in out Standard_Irreducible_Component_List ) is
 
    m : natural32 := 0;
    
  begin
    get(file,m);
    for i in 1..m loop
      declare
        c : Standard_Irreducible_Component;
      begin
        get(file,c);
        Append(first,last,c);
      end;
    end loop;
  end get;

  procedure get ( first,last : in out Multprec_Irreducible_Component_List ) is
  begin
    get(Standard_Input,first,last);
  end get;

  procedure get ( file : in file_type;
                  first,last : in out Multprec_Irreducible_Component_List ) is
 
    m : natural32 := 0;
    
  begin
    get(file,m);
    for i in 1..m loop
      declare
        c : Multprec_Irreducible_Component;
      begin
        get(file,c);
        Append(first,last,c);
      end;
    end loop;
  end get;

  procedure put ( L : in Standard_Irreducible_Component_List ) is
  begin
    put(Standard_Output,l);
  end put;

  procedure put ( file : file_type;
                  L : in Standard_Irreducible_Component_List ) is

    m : constant natural32 := Length_Of(L);
    tmp : Standard_Irreducible_Component_List := L;

  begin
    put(file,m,1); new_line(file);
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( L : in Multprec_Irreducible_Component_List ) is
  begin
    put(Standard_Output,L);
  end put;

  procedure put ( file : file_type;
                  L : in Multprec_Irreducible_Component_List ) is

    m : constant natural32 := Length_Of(L);
    tmp : Multprec_Irreducible_Component_List := L;

  begin
    put(file,m,1); new_line(file);
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

-- MINIMAL DATA :

  procedure get_labels 
               ( first,last : in out Standard_Irreducible_Component_List ) is
  begin
    get_labels(Standard_Input,first,last);
  end get_labels;

  procedure get_labels 
               ( file : in file_type;
                 first,last : in out Standard_Irreducible_Component_List ) is

    len : natural32 := 0;

  begin
    get(file,len);
    for i in 1..len loop
      declare
        c : Standard_Irreducible_Component;
      begin
        get_labels(file,c);
        Append(first,last,c);
      end;
    end loop;
  end get_labels;

  procedure get_labels 
               ( first,last : in out Multprec_Irreducible_Component_List ) is
  begin
    get_labels(Standard_Input,first,last);
  end get_labels;

  procedure get_labels 
               ( file : in file_type;
                 first,last : in out Multprec_Irreducible_Component_List ) is

    len : natural32 := 0;

  begin
    get(file,len);
    for i in 1..len loop
      declare
        c : Multprec_Irreducible_Component;
      begin
        get_labels(file,c);
        Append(first,last,c);
      end;
    end loop;
  end get_labels;

  procedure put_labels ( L : in Standard_Irreducible_Component_List ) is
  begin
    put(Standard_Output,L);
  end put_labels;

  procedure put_labels ( file : in file_type;
                         L : in Standard_Irreducible_Component_List ) is

    tmp : Standard_Irreducible_Component_List := L;

  begin
    put(file,Length_Of(L),1); new_line(file);
    while not Is_Null(tmp) loop
      put_labels(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put_labels;

end Irreducible_Component_Lists_io;
