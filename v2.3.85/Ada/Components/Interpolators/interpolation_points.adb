with unchecked_deallocation;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Interpolation_Points is

-- DATA STRUCTURES :

  type Standard_Sample_Node_Rep ( k : integer32 ) is record
    v : Standard_Complex_Vectors.Vector(1..k);
    s : Standard_Sample;
  end record;
  type Multprec_Sample_Node_Rep ( k : integer32 ) is record
    v : Multprec_Complex_Vectors.Vector(1..k);
    s : Multprec_Sample;
  end record;

-- CREATORS :

  function Create ( spt : Standard_Sample;
                    prospt : Standard_Complex_Vectors.Vector )
                  return Standard_Sample_Node is

    res : Standard_Sample_Node;
    res_rep : Standard_Sample_Node_Rep(prospt'last);

  begin
    res_rep.v := prospt;
    res_rep.s := spt;
    res := new Standard_Sample_Node_Rep'(res_rep);
    return res;
  end Create;

  function Create ( spt : Multprec_Sample;
                    prospt : Multprec_Complex_Vectors.Vector )
                  return Multprec_Sample_Node is

    res : Multprec_Sample_Node;
    res_rep : Multprec_Sample_Node_Rep(prospt'last);

  begin
    res_rep.v := prospt;
    res_rep.s := spt;
    res := new Multprec_Sample_Node_Rep'(res_rep);
    return res;
  end Create;

  procedure Update ( snd : in out Standard_Sample_Node;
                     prospt : in Standard_Complex_Vectors.Vector ) is
  begin
    snd.v := prospt;
  end Update;

  procedure Update ( snd : in out Multprec_Sample_Node;
                     prospt : in Multprec_Complex_Vectors.Vector ) is
  begin
    Multprec_Complex_Vectors.Clear(snd.v);
    snd.v := prospt;
  end Update;

-- SELECTORS :

  function Empty ( sn : Standard_Sample_Node ) return boolean is
  begin
    return (sn = null);
  end Empty;

  function Empty ( sn : Multprec_Sample_Node ) return boolean is
  begin
    return (sn = null);
  end Empty;

  function Sample_Node ( sn : Standard_Sample_Node )
                       return Standard_Complex_Vectors.Vector is
  begin
    return sn.v;
  end Sample_Node;

  function Sample_Node ( sn : Multprec_Sample_Node )
                       return Multprec_Complex_Vectors.Vector is
  begin
    return sn.v;
  end Sample_node;

  function Sample_Point ( sn : Standard_Sample_Node ) return Standard_Sample is
  begin
    return sn.s;
  end Sample_Point;

  function Sample_Point ( sn : Multprec_Sample_Node ) return Multprec_Sample is
  begin
    return sn.s;
  end Sample_Point;

  function Sample_Vector ( sn : Standard_Sample_Node )
                         return Standard_Complex_Vectors.Vector is
  begin
    return Sample_Point(sn.s).v;
  end Sample_Vector;

  function Sample_Vector ( sn : Multprec_Sample_Node )
                         return Multprec_Complex_Vectors.Vector is
  begin
    return Sample_Point(sn.s).v;
  end Sample_Vector;

-- DESTRUCTORS :

  procedure Shallow_Clear ( sn : in out Standard_Sample_Node ) is

    procedure free is
      new unchecked_deallocation(Standard_Sample_Node_Rep,
                                 Standard_Sample_Node);

  begin
    if sn /= null
     then free(sn);
    end if;
  end Shallow_Clear;
 
  procedure Shallow_Clear ( sn : in out Multprec_Sample_Node ) is

    procedure free is
      new unchecked_deallocation(Multprec_Sample_Node_Rep,
                                 Multprec_Sample_Node);

  begin
    if sn /= null
     then Multprec_Complex_Vectors.Clear(sn.v);
          free(sn);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( sn : in out Standard_Sample_Node ) is
  begin
    if sn /= null
     then Deep_Clear(sn.s);
          Shallow_Clear(sn);
    end if;
  end Deep_Clear;

  procedure Deep_Clear ( sn : in out Multprec_Sample_Node ) is
  begin
    if sn /= null
     then Deep_Clear(sn.s);
          Shallow_Clear(sn);
    end if;
  end Deep_Clear;

end Interpolation_Points;
