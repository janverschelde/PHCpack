with unchecked_deallocation;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with Monomial_Hashing;                  use Monomial_Hashing;

package body QuadDobl_Jacobian_Trees is

-- CONSTRUCTORS :

  function Initialize ( a : Jaco_Mat ) return Node is

    res : Node(a'last(2));

  begin
    res.a := new Jaco_Mat'(a);
    for i in res.d'range loop
      res.d(i) := null;
    end loop;
    return res;
  end Initialize;

  function Initialize ( a : Jaco_Mat ) return Eval_Node is

    res : Eval_Node(a'last(2));

  begin
    res.a := new Jaco_Mat'(a);
    res.f := new Eval_Jaco_Mat'(Create(a));
    for i in res.d'range loop
      res.d(i) := null;
    end loop;
    return res;
  end Initialize;

  procedure Create ( nd : in out Node ) is
  begin
    if nd.a /= null then
      for i in 1..nd.n loop 
        if ((nd.a = null) or else (Degree(nd.a.all,i) <= 0)) then
          nd.d(i) := null;
        else
          declare
            da : constant Link_to_Jaco_Mat := Diff(nd.a,i);
            child : Node(nd.n);
          begin
            child.a := da;
            Create(child);
            nd.d(i) := new Node'(child);
          end;
        end if;
      end loop;
    end if;
  end Create;

  procedure Create ( nd : in out Eval_Node ) is
  begin
    if nd.a /= null then
      for i in 1..nd.n loop 
        if ((nd.a = null) or else (Degree(nd.a.all,i) <= 0)) then
          nd.d(i) := null;
        else
          declare
            da : constant Link_to_Jaco_Mat := Diff(nd.a,i);
            child : Eval_Node(nd.n);
          begin
            child.a := da;
            child.f := new Eval_Jaco_Mat'(Create(da.all));
            Create(child);
            nd.d(i) := new Eval_Node'(child);
          end;
        end if;
      end loop;
    end if;
  end Create;

  function Create ( a : Jaco_Mat ) return Node is

    res : Node(a'last(2));

  begin
    res.a := new Jaco_Mat'(a);
    Create(res);
    return res;
  end Create;

  function Create ( a : Jaco_Mat ) return Eval_Node is

    res : Eval_Node(a'last(2));

  begin
    res.a := new Jaco_Mat'(a);
    res.f := new Eval_Jaco_Mat'(Create(a));
    Create(res);
    return res;
  end Create;

  function Diff ( a : Jaco_Mat; k : integer32 ) return Jaco_Mat is

    res : Jaco_Mat(a'range(1),a'range(2));

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        res(i,j) := Diff(a(i,j),k);
      end loop;
    end loop;
    return res;
  end Diff;

  function Diff ( a : Link_to_Jaco_Mat; k : integer32 ) 
                return Link_to_Jaco_Mat is

    res : Link_to_Jaco_Mat := null;

  begin
    if a /= null
     then res := new Jaco_Mat'(Diff(a.all,k));
    end if;
    return res;
  end Diff;

-- SELECTORS :

  function Degree ( a : Jaco_Mat; k : integer32 ) return integer32 is

    res : integer32 := -1;
    deg : integer32;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        deg := Degree(a(i,j),k);
        if deg > res
         then res := deg;
        end if;
      end loop;
    end loop;
    return res;
  end Degree;

  procedure Derivative ( nd : in out Node;
                         v : in Standard_Natural_Vectors.Vector;
                         da : out Link_to_Jaco_Mat ) is

    allzero : boolean := true;

  begin
    for i in v'range loop
      if v(i) > 0 then
        allzero := false;
        if ((nd.a = null) or else (Degree(nd.a.all,i) <= 0)) then
          da := null;
        else
          declare
            nv : Standard_Natural_Vectors.Vector(v'range) := v;
          begin
            nv(i) := nv(i)-1;
            if nd.d(i) = null then
              declare
                child : Node(nd.n);
              begin
                child.a := Diff(nd.a,i);
                nd.d(i) := new Node'(child);
              end;
            end if;
            Derivative(nd.d(i).all,nv,da);
          end;
        end if;
      end if;
      exit when not allzero;
    end loop;
    if allzero
     then da := nd.a;
    end if;
  end Derivative;

  procedure Derivative ( nd : in out Eval_Node;
                         v : in Standard_Natural_Vectors.Vector;
                         da : out Link_to_Jaco_Mat;
                         df : out Link_to_Eval_Jaco_Mat ) is

    allzero : boolean := true;

  begin
    for i in v'range loop
      if v(i) > 0 then
        allzero := false;
        if ((nd.a = null) or else (Degree(nd.a.all,i) <= 0)) then
          da := null;
          df := null;
        else
          declare
            nv : Standard_Natural_Vectors.Vector(v'range) := v;
          begin
            nv(i) := nv(i)-1;
            if nd.d(i) = null then
              declare
                child : Eval_Node(nd.n);
              begin
                child.a := Diff(nd.a,i);
                child.f := new Eval_Jaco_Mat'(Create(child.a.all));
                nd.d(i) := new Eval_Node'(child);
              end;
            end if;
            Derivative(nd.d(i).all,nv,da,df);
          end;
        end if;
      end if;
      exit when not allzero;
    end loop;
    if allzero
     then da := nd.a;
          df := nd.f;
    end if;
  end Derivative;

  procedure Derivative ( nd : in out Node; k : in integer32 ) is
  begin
    if k > 0 then
      for i in nd.d'range loop
        if ((nd.a /= null) and then (Degree(nd.a.all,i) > 0)) then
          declare
            child : Node(nd.n);
          begin
            child.a := Diff(nd.a,i);
            nd.d(i) := new Node'(child);
            if k > 1 
             then Derivative(nd.d(i).all,k-1);
            end if;
          end;
        end if;
      end loop;
    end if;
  end Derivative;

  procedure Derivative ( nd : in out Eval_Node; k : in integer32 ) is
  begin
    if k > 0 then
      for i in nd.d'range loop
        if ((nd.a /= null) and then (Degree(nd.a.all,i) > 0)) then
          declare
            child : Eval_Node(nd.n);
          begin
            child.a := Diff(nd.a,i);
            child.f := new Eval_Jaco_Mat'(Create(child.a.all));
            nd.d(i) := new Eval_Node'(child);
            if k > 1 
             then Derivative(nd.d(i).all,k-1);
            end if;
          end;
        end if;
      end loop;
    end if;
  end Derivative;

  procedure Dimensions ( v : in Link_to_VecMat; rows,cols : out integer32 ) is

    A : Link_to_Matrix;

  begin
    if v /= null then
      for i in v'range loop
        A := v(i);
        if A /= null then
          rows := A'last(1);
          cols := A'last(2);
          return;
        end if;
      end loop;
    end if;
    rows := 0;
    cols := 0;
  end Dimensions;

-- ENUMERATORS :

  procedure Enumerate_Nodes ( nd : in Node ) is

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
  end Enumerate_Nodes;

  procedure Enumerate_Eval_Nodes ( nd : in Eval_Node ) is

    continue : boolean := true;

    procedure Enumerate ( current : in Eval_Node ) is
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
  end Enumerate_Eval_Nodes;

-- EVALUATORS :

  procedure Create_Remember_Derivatives
               ( a : in Jaco_Mat; k : in integer32;
                 nd : out Link_to_Eval_Node ) is

    n : constant integer32 := a'last(2);

    procedure Create_Derivative
                ( d : in Standard_Natural_Vectors.Vector ) is

      da : Link_to_Jaco_Mat;
      df : Link_to_Eval_Jaco_Mat;

    begin
      Derivative(nd.all,d,da,df);
    end Create_Derivative;

    procedure Create_Remember_Table is 
      new Enumerate_Monomials(Create_Derivative);

  begin
    nd := new Eval_Node'(Initialize(a));
    Create_Remember_Table(natural32(k),natural32(n));
  end Create_Remember_Derivatives;

  function Evaluate_Jacobian_Remember_Table
               ( nd : Link_to_Eval_Node; k,n : integer32;
                 x : QuadDobl_Complex_Vectors.Vector )
               return Jacobian_Remember_Table is        

    res : Jacobian_Remember_Table(0..k);
    deg,ind : integer32;

    procedure Evaluate_Derivative
                ( m : in Standard_Natural_Vectors.Vector ) is

      da : Link_to_Jaco_Mat;
      df : Link_to_Eval_Jaco_Mat := null;

    begin
      ind := ind + 1;
      Derivative(nd.all,m,da,df);
      if df = null
       then res(deg)(ind) := null;
       else res(deg)(ind) := new Matrix'(Eval(df.all,x));
      end if;
    end Evaluate_Derivative;
    procedure Evaluate_Derivatives is
      new Enumerate_Monomials(Evaluate_Derivative);
     
  begin
    for i in 0..k loop
      declare
        mci : constant natural32 := Monomial_Count(natural32(i),natural32(n));
        vmi : VecMat(1..integer32(mci));
      begin
        res(i) := new VecMat'(vmi);
      end;
      deg := i;
      ind := 0;
      Evaluate_Derivatives(natural32(i),natural32(n));
    end loop;
    return res;
  end Evaluate_Jacobian_Remember_Table;

-- DESTRUCTORS :

  procedure Clear ( nd : in out Node ) is
  begin
    Clear(nd.a);
    Clear(nd.d);
  end Clear;

  procedure Clear ( nd : in out Eval_Node ) is
  begin
    Clear(nd.a);
    Clear(nd.f);
    Clear(nd.d);
  end Clear;

  procedure Clear ( nd : in out Link_to_Node ) is

    procedure free is new unchecked_deallocation(Node,Link_to_Node);

  begin
    if nd /= null then
      Clear(nd.all);
      free(nd);
    end if;
  end Clear;

  procedure Clear ( nd : in out Link_to_Eval_Node ) is

    procedure free is new unchecked_deallocation(Eval_Node,Link_to_Eval_Node);

  begin
    if nd /= null then
      Clear(nd.all);
      free(nd);
    end if;
  end Clear;

  procedure Clear ( nd : in out Array_of_Nodes ) is
  begin
    for i in nd'range loop
      Clear(nd(i));
    end loop;
  end Clear;

  procedure Clear ( nd : in out Array_of_Eval_Nodes ) is
  begin
    for i in nd'range loop
      Clear(nd(i));
    end loop;
  end Clear;

  procedure Clear ( jrt : in out Jacobian_Remember_Table ) is
  begin
    for i in jrt'range loop
      QuadDobl_Complex_VecMats.Deep_Clear(jrt(i));
    end loop;
  end Clear;

end QuadDobl_Jacobian_Trees;
