with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Brackets_io;                       use Brackets_io;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors; 
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io; 

package body Localization_Posets_io is

  procedure put ( top,bottom : in Bracket; roco : in natural32 ) is
  begin
    put(Standard_Output,top,bottom,roco);
  end put;

  procedure put ( file : in file_type;
                  top,bottom : in Bracket; roco : in natural32 ) is
  begin
    put(file,"("); put(file,top);
    put(file,","); put(file,bottom); put(file,","); put(file,roco,1);
    put(file,")");
  end put;

  procedure put ( file : in file_type;
                  top,bottom : in Bracket;
                  roco : in natural32; label : in integer32;
                  child_labels : in Link_to_Vector ) is
  begin
    put(file,"("); put(file,label,1); put(file,",");
    put(file,top); put(file,",");
    put(file,bottom); put(file,","); put(file,roco,1);
    put(file,",{");
    if child_labels /= null
     then put(file,child_labels.all);
    end if;
    put(file,"}");
    put(file,")");
  end put;

  procedure put ( root : in Node; lvl : in integer32 ) is
  begin
    put(Standard_Output,root,lvl);
  end put;

  procedure put ( file : in file_type; root : in Node; lvl : in integer32 ) is

    lroot : constant Link_to_Node := new Node'(root);
    lvlnd : constant Link_to_Node := Find_Node(lroot,lvl);

    procedure Write_Node ( nd : in Node; continue : out boolean ) is
    begin
      put(file,nd.top,nd.bottom,natural32(nd.roco));
      continue := true;
    end Write_Node;
    procedure Write_Siblings is new Enumerate_Siblings(Write_Node);

  begin
    if lvlnd /= null
     then Write_Siblings(lvlnd.all);
    end if;
  end put;

  procedure put ( poset : in Node ) is
  begin
    put(Standard_Output,poset);
  end put;

  procedure put ( file : in file_type; poset : in Node ) is

    np : natural32;

  begin
    if poset.level < 10
     then np := 1;
     else np := 2;
    end if;
    for i in 0..poset.level loop
      put(file,"n = "); put(file,i,np); put(file," : ");
      put(file,poset,i); new_line(file);
    end loop;
  end put;

  procedure put ( poset : in Array_of_Nodes ) is
  begin
    put(Standard_Output,poset);
  end put;

  procedure put ( file : in file_type; poset : in Array_of_Nodes ) is

    np : natural32;

    procedure Write_Node ( nd : in Node; continue : out boolean ) is
    begin
      put(file,nd.top,nd.bottom,natural32(nd.roco));
      continue := true;
    end Write_Node;
    procedure Write_Siblings is new Enumerate_Siblings(Write_Node);

  begin
    if poset'last < 10
     then np := 1;
     else np := 2;
    end if;
    for i in poset'range loop
      put(file,"n = "); put(file,i,np); put(file," : ");
     -- put(file,"l = "); put(file,Number_of_Siblings(poset(i)),np);
     -- put(file," : ");
      if poset(i) /= null
       then Write_Siblings(poset(i).all);
      end if;
      new_line(file);
    end loop;
  end put;

  procedure put ( poset : in Array_of_Array_of_Nodes ) is
  begin
    put(Standard_Output,poset);
  end put;

  procedure put ( file : in file_type; poset : in Array_of_Array_of_Nodes ) is

    np : natural32;
    lnd : Link_to_Node;

  begin
    if poset'last < 10
     then np := 1;
     else np := 2;
    end if;
    for i in poset'range loop
      put(file,"n = "); put(file,i,np); put(file," : ");
      if poset(i) /= null then
        for j in poset(i)'range loop
          lnd := poset(i)(j);
          put(file,lnd.top,lnd.bottom,natural32(lnd.roco),
              lnd.label,lnd.child_labels);
        end loop;
      end if;
      new_line(file);
    end loop;
  end put;

  procedure put_roco ( poset : in Array_of_Array_of_Nodes ) is
  begin
    put_roco(Standard_Output,poset);
  end put_roco;

  procedure put_roco
              ( file : in file_type; poset : in Array_of_Array_of_Nodes ) is

    np : natural32;
    lnd : Link_to_Node;

  begin
    if poset'last < 10
     then np := 1;
     else np := 2;
    end if;
    for i in poset'range loop
      put(file,"n = "); put(file,i,np); put(file," : ");
      if poset(i) /= null then
        for j in poset(i)'range loop
          lnd := poset(i)(j);
          if lnd = null
           then put(file," 0");
           else put(file," "); put(file,lnd.roco,1);
          end if;
        end loop;
      end if;
      new_line(file);
    end loop;
  end put_roco;

end Localization_Posets_io;
