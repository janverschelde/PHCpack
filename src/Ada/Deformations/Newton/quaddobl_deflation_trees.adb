with unchecked_deallocation;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with QuadDobl_Deflate_Singularities;    use QuadDobl_Deflate_Singularities;

package body QuadDobl_Deflation_Trees is

-- CREATORS :

  function Create_Root ( p : Poly_Sys ) return Node is

    ne : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'last)));
    res : Node(ne,nv);

  begin
    res.d := 0;
    Copy(p,res.s);
    res.f := Create(p);
    res.jm := Create(p);
    res.jf := Create(res.jm);
    return res;
  end Create_Root;

  procedure Create_Child ( nd : in out Node; m : in integer32 ) is
 
    df : constant Poly_Sys := Deflate(nd.s,nd.jm,natural32(m));
    ne : constant integer32 := df'last;
    nv : constant integer32 := nd.nv+m;
    child : Node(ne,nv);
 
  begin
    child.d := nd.d+1;
    child.s := df;
    child.f := Create(df);
    child.jm := Create(df);
    child.jf := Create(child.jm);
    nd.children(m) := new Node'(child);
  end Create_Child;

-- DESTRUCTORS :

  procedure Clear ( nd : in out Node ) is
  begin
    Clear(nd.s);
    Clear(nd.f);
    Clear(nd.jm);
    Clear(nd.jf);
    Clear(nd.children);
    Clear(nd.sols);
  end Clear;

  procedure Clear ( nd : in out Link_to_Node ) is

    procedure free is
      new unchecked_deallocation(Node,Link_to_Node);
 
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
 
end QuadDobl_Deflation_Trees;
