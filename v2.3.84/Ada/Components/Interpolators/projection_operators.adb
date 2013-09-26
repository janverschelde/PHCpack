with unchecked_deallocation;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
--with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Standard_Linear_Projections;        use Standard_Linear_Projections;
with Multprec_Linear_Projections;        use Multprec_Linear_Projections;
with Standard_Central_Projections;       use Standard_Central_Projections;
with Multprec_Central_Projections;       use Multprec_Central_Projections;

package body Projection_Operators is

-- DATA STRUCTURES :
--   The base points are stored as sample points to allow later refinement.
--   With each base point corresponds a hyperplane that does not contain
--   the base point.  We also store the projected base points.

  type Standard_Projector_Rep ( n,nb : integer32 ) is record
    hyps : Standard_Complex_VecVecs.VecVec(1..n);
    base,base_last : Standard_Sample_List;
    len : integer32;
    prbase : Standard_Complex_VecVecs.VecVec(1..nb);
    bahyps : Standard_Complex_VecVecs.VecVec(1..nb); 
  end record;

  type Multprec_Projector_Rep ( n,nb : integer32 ) is record
    hyps : Multprec_Complex_VecVecs.VecVec(1..n);
    base,base_last : Multprec_Sample_List;
    len : integer32;
    prbase : Multprec_Complex_VecVecs.VecVec(1..nb);
    bahyps : Multprec_Complex_VecVecs.VecVec(1..nb);
  end record;

-- CREATORS :

  function Create ( hyps : Standard_Complex_VecVecs.VecVec )
                  return Standard_Projector is

    res : Standard_Projector;
    k : constant integer32 := hyps'last;
    n : constant integer32 := hyps(hyps'first)'last; 
    res_rep : Standard_Projector_Rep(k,n);

  begin
    res_rep.hyps := hyps;
    res_rep.len := 0;
    res := new Standard_Projector_Rep'(res_rep);
    return res;
  end Create;

  function Create ( hyps : Multprec_Complex_VecVecs.VecVec )
                  return Multprec_Projector is

    res : Multprec_Projector;
    k : constant integer32 := hyps'last;
    n : constant integer32 := hyps(hyps'first)'last;
    res_rep : Multprec_Projector_Rep(k,n);

  begin
    res_rep.hyps := hyps;
    res_rep.len := 0;
    res := new Multprec_Projector_Rep'(res_rep);
    return res;
  end Create;

  function Create ( k,n : natural32 ) return Standard_Projector is

    hyps : Standard_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in 1..integer32(k) loop
      hyps(i) := new Standard_Complex_Vectors.Vector'
                       (Random_Vector(0,integer32(n)));
    end loop;
    return Create(hyps);
  end Create;

  function Create ( k,n,size : natural32 ) return Multprec_Projector is

    hyps : Multprec_Complex_VecVecs.VecVec(1..integer32(k));
    ranv : Standard_Complex_Vectors.Vector(0..integer32(n));

  begin
    for i in 1..integer32(k) loop
   -- hyps(i) := new Multprec_Complex_Vectors.Vector'(Random_Vector(0,n,size));
      ranv := Random_Vector(0,integer32(n));
      hyps(i) := new Multprec_Complex_Vectors.Vector'(Create(ranv));
    end loop;
    return Create(hyps);
  end Create;

  procedure Update ( p : in out Standard_Projector;
                     basept : in Standard_Sample;
                     basehyp : in Standard_Complex_Vectors.Vector ) is

    basevc : constant Standard_Complex_Vectors.Vector
           := Sample_Point(basept).v;
    projpt : Standard_Complex_Vectors.Vector(basevc'range);

  begin
    if p /= null then
      Append(p.base,p.base_last,basept);
      p.len := p.len+1;
      p.bahyps(p.len) := new Standard_Complex_Vectors.Vector'(basehyp);
      if p.len = 1 then
        p.prbase(p.len) := new Standard_Complex_Vectors.Vector'(basevc);
      else
        projpt := Intersect(p.bahyps(1).all,p.prbase(1).all,basevc,basevc'last);
        for i in 2..p.len-1 loop
          projpt := Intersect(p.bahyps(i).all,p.prbase(i).all,
                              projpt,projpt'last);
        end loop;
        p.prbase(p.len) := new Standard_Complex_Vectors.Vector'(projpt);
      end if;
    end if;
  end Update;

  procedure Update ( p : in out Multprec_Projector;
                     basept : in Multprec_Sample;
                     basehyp : in Multprec_Complex_Vectors.Vector ) is

    basevc : constant Multprec_Complex_Vectors.Vector
           := Sample_Point(basept).v;
    wrk,projpt : Multprec_Complex_Vectors.Vector(basevc'range);

  begin
    if p /= null then
      Append(p.base,p.base_last,basept);
      p.len := p.len+1;
      p.bahyps(p.len) := new Multprec_Complex_Vectors.Vector'(basehyp);
      if p.len = 1 then
        p.prbase(p.len) := new Multprec_Complex_Vectors.Vector'(basevc);
      else
        projpt := Intersect(p.bahyps(1).all,p.prbase(1).all,basevc,basevc'last);
        for i in 2..p.len-1 loop
          Multprec_Complex_Vectors.Copy(projpt,wrk);
          Multprec_Complex_Vectors.Clear(projpt);
          projpt := Intersect(p.bahyps(i).all,p.prbase(i).all,wrk,wrk'last);
        end loop;
        p.prbase(p.len) := new Multprec_Complex_Vectors.Vector'(projpt);
        Multprec_Complex_Vectors.Clear(wrk);
      end if;
    end if;
  end Update;

-- SELECTORS :

  function Empty ( p : Standard_Projector ) return boolean is
  begin
    return (p = null);
  end Empty;

  function Empty ( p : Multprec_Projector ) return boolean is
  begin
    return (p = null);
  end Empty;

  function Origin_Dimension ( p : Standard_Projector ) return natural32 is
  begin
    if p = null
     then return 0;
     else return natural32(p.hyps(1)'last);
    end if;
  end Origin_Dimension;

  function Origin_Dimension ( p : Multprec_Projector ) return natural32 is
  begin
    if p = null
     then return 0;
     else return natural32(p.hyps(1)'last);
    end if;
  end Origin_Dimension;

  function Target_Dimension ( p : Standard_Projector ) return natural32 is
  begin
    if p = null
     then return 0;
     else return natural32(p.n);
    end if;
  end Target_Dimension;

  function Target_Dimension ( p : Multprec_Projector ) return natural32 is
  begin
    if p = null
     then return 0;
     else return natural32(p.n);
    end if;
  end Target_Dimension;

  function Hyperplanes ( p : Standard_Projector )
                       return Standard_Complex_VecVecs.VecVec is
  begin
    return p.hyps;
  end Hyperplanes;

  function Hyperplanes ( p : Multprec_Projector )
                       return Multprec_Complex_VecVecs.VecVec is
  begin
    return p.hyps;
  end Hyperplanes;

  function Is_Central ( p : Standard_Projector ) return boolean is
  begin
    return not Is_Null(p.base);
  end Is_Central;

  function Is_Central ( p : Multprec_Projector ) return boolean is
  begin
    return not Is_Null(p.base);
  end Is_Central;

  function Centrality ( p : Standard_Projector ) return natural32 is
  begin
    return natural32(p.len);
  end Centrality;

  function Centrality ( p : Multprec_Projector ) return natural32 is
  begin
    return natural32(p.len);
  end Centrality;

  function Base_Points ( p : Standard_Projector )
                       return Standard_Sample_List is

    res : Standard_Sample_List;

  begin
    if p /= null
     then res := p.base;
    end if;
    return res;
  end Base_Points;

  function Base_Points ( p : Multprec_Projector )
                       return Multprec_Sample_List is

    res : Multprec_Sample_List;

  begin
    if p /= null
     then res := p.base;
    end if;
    return res;
  end Base_Points;

-- PROJECTORS ON BASIC DATA :

  function Project ( p : Standard_Projector;
                     x : Standard_Complex_Vectors.Vector )
                   return Standard_Complex_Vectors.Vector is
  begin
    if p = null then
      return x;
    elsif Is_Null(p.base) then
      return Evaluate(p.hyps,x,p.n);
    else
      return Intersect(p.bahyps(1..p.len),p.prbase(1..p.len),x,p.n);
    end if;
  end Project;

  function Project ( p : Standard_Projector;
                     x : Standard_Complex_VecVecs.VecVec )
                   return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(x'range);

  begin
    if p = null then
      return x;
    elsif Is_Null(p.base) then
      return Evaluate(p.hyps,x,p.n);
    else
      for i in x'range loop
        res(i) := new Standard_Complex_Vectors.Vector'
          (Intersect(p.bahyps(1..p.len),p.prbase(1..p.len),x(i).all,p.n));
      end loop;
      return res;
    end if;
  end Project;

  function Project ( p : Multprec_Projector;
                     x : Standard_Complex_Vectors.Vector )
                   return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(1..p.n);
    mpx : Multprec_Complex_Vectors.Vector(x'range) := Create(x);

  begin
    if p = null then
      return mpx;
    elsif Is_Null(p.base) then
      return Evaluate(p.hyps,x,p.n);
    else
      res := Intersect(p.bahyps(1..p.len),p.prbase(1..p.len),mpx,p.n);
      Multprec_Complex_Vectors.Clear(mpx);
      return res;
    end if;
  end Project;

  function Project ( p : Multprec_Projector;
                     x : Multprec_Complex_Vectors.Vector )
                   return Multprec_Complex_Vectors.Vector is
  begin
    if p = null then
      return x;
    elsif Is_Null(p.base) then
      return Evaluate(p.hyps,x,p.n);
    else
      return Intersect(p.bahyps(1..p.len),p.prbase(1..p.len),x,p.n);
    end if;
  end Project;

  function Project ( p : Multprec_Projector;
                     x : Standard_Complex_VecVecs.VecVec )
                   return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(x'range);
    mpx : Multprec_Complex_Vectors.Vector(x(x'first)'range);

  begin
    if p = null then
      return Create(x);
    elsif Is_Null(p.base) then
      return Evaluate(p.hyps,x,p.n);
    else
      for i in x'range loop
        mpx := Create(x(i).all);
        res(i) := new Multprec_Complex_Vectors.Vector'
          (Intersect(p.bahyps(1..p.len),p.prbase(1..p.len),mpx,p.n));
        Multprec_Complex_Vectors.Clear(mpx);
      end loop;
      return res;
    end if;
  end Project;

  function Project ( p : Multprec_Projector;
                     x : Multprec_Complex_VecVecs.VecVec )
                   return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(x'range);

  begin
    if p = null then
      return x;
    elsif Is_Null(p.base) then
      return Evaluate(p.hyps,x,p.n);
    else
      for i in x'range loop
        res(i) := new Multprec_Complex_Vectors.Vector'
          (Intersect(p.bahyps(1..p.len),p.prbase(1..p.len),x(i).all,p.n));
      end loop;
      return res;
    end if;
  end Project;

-- PROJECTORS ON ENCAPSULATED DATA :

  function Project ( p : Standard_Projector; s : Standard_Sample )
                   return Standard_Sample_Node is

    res : Standard_Sample_Node;

  begin
    res := Create(s,Project(p,Sample_Point(s).v));
    return res;
  end Project;

  function Project ( p : Multprec_Projector; s : Multprec_Sample )
                   return Multprec_Sample_Node is

    res : Multprec_Sample_Node;

  begin
    res := Create(s,Project(p,Sample_Point(s).v));
    return res;
  end Project;

  procedure Project ( p : in Standard_Projector; s : in Standard_Sample_List;
                      sn,sn_last : in out Standard_Sample_Node_List ) is

    tmp : Standard_Sample_List := s;

  begin
    while not Is_Null(tmp) loop
      Append(sn,sn_last,Project(p,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
  end Project;

  procedure Project ( p : in Multprec_Projector; s : in Multprec_Sample_List;
                      sn,sn_last : in out Multprec_Sample_Node_List ) is

    tmp : Multprec_Sample_List := s;

  begin
    while not Is_Null(tmp) loop
      Append(sn,sn_last,Project(p,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
  end Project;

-- DESTRUCTORS :

  procedure free is
    new unchecked_deallocation(Standard_Projector_Rep,Standard_Projector);
  procedure free is
    new unchecked_deallocation(Multprec_Projector_Rep,Multprec_Projector);

  procedure Shallow_Clear ( p : in out Standard_Projector ) is
  begin
    if p /= null
     then free(p);
    end if;
  end Shallow_Clear;

  procedure Shallow_Clear ( p : in out Multprec_Projector ) is
  begin
    if p /= null
     then free(p);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( p : in out Standard_Projector ) is
  begin
    if p /= null
     then Deep_Clear(p.base);
          free(p);
    end if;
  end Deep_Clear;

  procedure Deep_Clear ( p : in out Multprec_Projector ) is
  begin
    if p /= null
     then Deep_Clear(p.base);
          free(p);
    end if;
  end Deep_Clear;

end Projection_Operators;
