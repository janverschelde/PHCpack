with unchecked_deallocation;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package body DoblDobl_Derivative_Trees is

-- CONSTRUCTORS :

  function Initialize ( p : Poly ) return Node is

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Node(n);

  begin
    Copy(p,res.p);
    for i in res.d'range loop
      res.d(i) := null;
    end loop;
    return res;
  end Initialize;

  function Initialize ( p : Poly_Sys ) return Array_of_Nodes is

    res : Array_of_Nodes(p'range);

  begin
    for i in p'range loop
      res(i) := new Node'(Initialize(p(i)));
    end loop;
    return res;
  end Initialize;

  procedure Create ( nd : in out Node ) is
  begin
    if nd.p /= Null_Poly then
      for i in 1..nd.n loop 
        if Degree(nd.p,i) <=  0 then
          nd.d(i) := null;
        else
          declare
            dp : constant Poly := Diff(nd.p,i);
            child : Node(nd.n);
          begin
            child.p := dp;
            Create(child);
            nd.d(i) := new Node'(child);
          end;
        end if;
      end loop;
    end if;
  end Create;

  function Create ( p : Poly ) return Node is

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Node(n);

  begin
    Copy(p,res.p);
    Create(res);
    return res;
  end Create;

  function Create ( p : Poly_Sys ) return Array_of_Nodes is

    res : Array_of_Nodes(p'range);

  begin
    for i in p'range loop
      res(i) := new Node'(Create(p(i)));
    end loop;
    return res;
  end Create;

-- SELECTORS :

  procedure Derivative ( nd : in out Node; v : in Vector; dp : out Poly ) is

    allzero : boolean := true;

  begin
    for i in v'range loop
      if v(i) > 0 then
        allzero := false;
        if Degree(nd.p,i) <= 0 then
          dp := Null_Poly;
        else
          declare
            nv : Vector(v'range) := v;
          begin
            nv(i) := nv(i)-1;
            if nd.d(i) = null then
              declare
                child : Node(nd.n);
              begin
                child.p := Diff(nd.p,i);
                nd.d(i) := new Node'(child);
              end;
            end if;
            Derivative(nd.d(i).all,nv,dp);
          end;
        end if;
      end if;
      exit when not allzero;
    end loop;
    if allzero
     then Copy(nd.p,dp);
    end if;
  end Derivative;

-- ENUMERATORS :

  procedure Enumerate ( nd : in Node ) is

    continue : boolean := true;

    procedure Enumerate ( current : in Node ) is
    begin
      Process(current,continue);
      if continue then
        for i in current.d'range loop
          if current.d(i) /= null
           then Enumerate(current.d(i).all);
          end if;
          exit when not continue;
        end loop;
      end if;
    end Enumerate;

  begin
    Enumerate(nd);
  end Enumerate;

  procedure Enumerate_Nodes ( nd : in Array_of_Nodes ) is

    continue : boolean := true;

    procedure Enumerate ( current : in Array_of_Nodes ) is
    begin
      for i in current'range loop
        if current(i) /= null then
          Process(current(i).all,continue);
          if continue
           then Enumerate(current(i).d);
          end if;
        end if;
        exit when not continue;
      end loop;
    end Enumerate;

  begin
    Enumerate(nd);
  end Enumerate_Nodes;

-- DESTRUCTORS :

  procedure Clear ( nd : in out Node ) is
  begin
    Clear(nd.p);
    Clear(nd.d);
  end Clear;

  procedure Clear ( nd : in out Link_to_Node ) is

    procedure free is new unchecked_deallocation(Node,Link_to_Node);

  begin
    if nd /= null
     then Clear(nd.all); free(nd);
    end if;
  end Clear;

  procedure Clear ( nd : in out Array_of_Nodes ) is
  begin
    for i in nd'range loop
      Clear(nd(i));
    end loop;
  end Clear;

end DoblDobl_Derivative_Trees;
