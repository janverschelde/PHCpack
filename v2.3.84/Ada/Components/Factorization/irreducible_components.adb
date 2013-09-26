with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;

package body Irreducible_Components is

-- IMPLEMENTION NOTE :
--   The main responsability of this data structure is to link the
--   interpolation filter with the span of the component.
--   In particular, the component test on interpolation filters
--   should only be done on restricted points for nonempty spans.

-- CREATORS :

  function Create ( f : Standard_Filter )
                  return Standard_Irreducible_Component is

    res : Standard_Irreducible_Component;

  begin
    res.f := f;
    return res;
  end Create;

  function Create ( f : Standard_Filter; s : Standard_Span )
                  return Standard_Irreducible_Component is

    res : Standard_Irreducible_Component;

  begin
    res.f := f;
    res.s := s;
    return res;
  end Create;

  function Create ( f : Multprec_Filter )
                  return Multprec_Irreducible_Component is

    res : Multprec_Irreducible_Component;

  begin
    res.f := f;
    return res;
  end Create;

  function Create ( f : Multprec_Filter; s : Multprec_Span )
                  return Multprec_Irreducible_Component is

    res : Multprec_Irreducible_Component;

  begin
    res.f := f;
    res.s := s;
    return res;
  end Create;

  procedure Initialize_Labels ( c : in out Standard_Irreducible_Component;
                                d : in integer32 ) is
  begin
    c.g := new Standard_Natural_Vectors.Vector'(1..d => 0);
  end Initialize_Labels;

  procedure Initialize_Labels ( c : in out Standard_Irreducible_Component;
                                lab : in Standard_Natural_Vectors.Vector ) is
  begin
    c.g := new Standard_Natural_Vectors.Vector'(lab);
  end Initialize_Labels;

  procedure Initialize_Labels ( c : in out Multprec_Irreducible_Component;
                                d : in integer32 ) is
  begin
    c.g := new Standard_Natural_Vectors.Vector'(1..d => 0);
  end Initialize_Labels;

  procedure Initialize_Labels ( c : in out Multprec_Irreducible_Component;
                                lab : in Standard_Natural_Vectors.Vector ) is
  begin
    c.g := new Standard_Natural_Vectors.Vector'(lab);
  end Initialize_Labels;

  procedure Add_Label ( v : in Standard_Natural_Vectors.Link_to_Vector;
                        i : in natural32; ind : out integer32 ) is

    use Standard_Natural_Vectors;

  begin
    if v = null then
      ind := 0;
    else
      ind := v'last+1;
      for j in v'range loop
        if v(j) = 0 then
          ind := j;
          v(j) := i;
          exit;
        end if;
      end loop;
    end if;
  end Add_Label;

  procedure Add_Label ( c : in out Standard_Irreducible_Component;
                        i : in natural32; ind : out integer32 ) is
  begin
    Add_Label(c.g,i,ind);
  end Add_Label;

  procedure Add_Label ( c : in out Multprec_Irreducible_Component;
                        i : in natural32; ind : out integer32 ) is
  begin
    Add_Label(c.g,i,ind);
  end Add_Label;

  procedure Add_Point ( c : in out Standard_Irreducible_Component;
                        spt : in Standard_Sample ) is
  begin
    Append(c.p,c.plast,spt);
  end Add_Point;

  procedure Add_Point ( c : in out Multprec_Irreducible_Component;
                        spt : in Standard_Sample ) is
  begin
    Append(c.p,c.plast,spt);
  end Add_Point;

  procedure Add_Points ( c : in out Standard_Irreducible_Component;
                         sps,sps_last : in Standard_Sample_List ) is
  begin
    c.p := sps;
    c.plast := sps_last;
  end Add_Points;

  procedure Add_Points ( c : in out Multprec_Irreducible_Component;
                         sps,sps_last : in Standard_Sample_List ) is
  begin
    c.p := sps;
    c.plast := sps_last;
  end Add_Points;

  procedure Select_Labeled_Points 
               ( c : in out Standard_Irreducible_Component;
                 sps : in Standard_Sample_List ) is

    lab : integer32;
    ind : natural32;
    tmp : Standard_Sample_List := sps;
    use Standard_Natural_Vectors;

  begin
    if c.g /= null then
      lab := c.g'first;
      ind := 1;
      while not Is_Null(tmp) loop
        if c.g(lab) = ind then
          Append(c.p,c.plast,Head_Of(tmp));
          lab := lab + 1;
        end if;
        exit when lab > c.g'last;
        tmp := Tail_Of(tmp);
        ind := ind + 1;
      end loop;
    end if;
  end Select_Labeled_Points;

  procedure Select_Labeled_Points 
               ( c : in out Multprec_Irreducible_Component;
                 sps : in Standard_Sample_List ) is

    ind : natural32;
    lab : integer32;
    tmp : Standard_Sample_List := sps;
    use Standard_Natural_Vectors;

  begin
    if c.g /= null then
      lab := c.g'first;
      ind := 1;
      while not Is_Null(tmp) loop
        if c.g(lab) = ind then
          Append(c.p,c.plast,Head_Of(tmp));
          lab := lab + 1;
        end if;
        exit when lab > c.g'last;
        tmp := Tail_Of(tmp);
        ind := ind + 1;
      end loop;
    end if;
  end Select_Labeled_Points;

  procedure Add_Filter ( c : in out Standard_Irreducible_Component;
                         f : in Standard_Filter ) is
  begin
    c.f := f;
  end Add_Filter;

-- SELECTORS :

  function Dimension ( c : Standard_Irreducible_Component ) return natural32 is
  begin
    return Dimension(c.f);
  end Dimension;

  function Dimension ( c : Multprec_Irreducible_Component ) return natural32 is
  begin
    return Dimension(c.f);
  end Dimension;

  function Degree ( c : Standard_Irreducible_Component ) return natural32 is

    df : constant natural32 := Degree(c.f);
    use Standard_Natural_Vectors;

  begin
    if df = 0 then
      if c.g = null
       then return 0;
       else return c.g'length;
      end if;
    else
      return df + Centrality(c.f);
    end if;
  end Degree;

  function Degree ( c : Multprec_Irreducible_Component ) return natural32 is

    df : constant natural32 := Degree(c.f);
    use Standard_Natural_Vectors;

  begin
    if df = 0 then
      if c.g = null
       then return 0;
       else return c.g'length;
      end if;
    else
      return df + Centrality(c.f);
    end if;
  end Degree;

  function Filter ( c : Standard_Irreducible_Component )
                  return Standard_Filter is
  begin
    return c.f;
  end Filter;

  function Filter ( c : Multprec_Irreducible_Component )
                  return Multprec_Filter is
  begin
    return c.f;
  end Filter;

  function Span ( c : Standard_Irreducible_Component ) return Standard_Span is
  begin
    return c.s;
  end Span;

  function Span ( c : Multprec_Irreducible_Component ) return Multprec_Span is
  begin
    return c.s;
  end Span;

  function Labels ( c : Standard_Irreducible_Component )
                  return Standard_Natural_Vectors.Link_to_Vector is
  begin
    return c.g;
  end Labels;

  function Labels ( c : Multprec_Irreducible_Component )
                  return Standard_Natural_Vectors.Link_to_Vector is
  begin
    return c.g;
  end Labels;

  function Points ( c : Standard_Irreducible_Component )
                  return Standard_Sample_List is
  begin
    return c.p;
  end Points;

  function Points ( c : Multprec_Irreducible_Component )
                  return Standard_Sample_List is
  begin
    return c.p;
  end Points;

-- MEMBERSHIP TESTS :

  function On_Component ( c : Standard_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean is
  begin
    if Empty(c.s) then
      return On_Component(c.f,tol,x);
    elsif not In_Span(c.s,tol,x) then
      return false;
    else
      return On_Component(c.f,tol,Restrict(c.s,0,x));
    end if;
  end On_Component;

  function On_Component ( file : file_type;
                          c : Standard_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean is
  begin
    if Empty(c.s) then
      return On_Component(file,c.f,tol,x);
    elsif not In_Span(file,c.s,tol,x) then
      return false;
    else
      return On_Component(file,c.f,tol,Restrict(c.s,0,x));
    end if;
  end On_Component;

  function On_Component ( c : Multprec_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean is
  begin
    if Empty(c.s) then
      return On_Component(c.f,tol,x);
    elsif not In_Span(c.s,tol,x) then
      return false;
    else
      return On_Component(c.f,tol,Restrict(c.s,0,x));
    end if;
  end On_Component;

  function On_Component ( file : file_type;
                          c : Multprec_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean is
  begin
    if Empty(c.s) then
      return On_Component(file,c.f,tol,x);
    elsif not In_Span(file,c.s,tol,x) then
      return false;
    else
      return On_Component(file,c.f,tol,Restrict(c.s,0,x));
    end if;
  end On_Component;

  function On_Component ( c : Multprec_Irreducible_Component;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return boolean is
  begin
    if Empty(c.s) then
      return On_Component(c.f,tol,x);
    elsif not In_Span(c.s,tol,x) then
      return false;
    else
      return On_Component(c.f,tol,Restrict(c.s,0,x));
    end if;
  end On_Component;

  function On_Component ( file : file_type;
                          c : Multprec_Irreducible_Component;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return boolean is
  begin
    if Empty(c.s) then
      return On_Component(file,c.f,tol,x);
    elsif not In_Span(file,c.s,tol,x) then
      return false;
    else
      return On_Component(file,c.f,tol,Restrict(c.s,0,x));
    end if;
  end On_Component;

  function Homotopy_Membership_Test
              ( c : Standard_Irreducible_Component;
                x : Standard_Complex_Vectors.Vector; tol : double_float )
              return boolean is

    sps,sps_last : Standard_Sample_List;
    isin : natural32;

  begin
    Membership_Test(c.p,x,tol,isin,sps,sps_last);
    Deep_Clear(sps);
    return (isin > 0);
  end Homotopy_Membership_Test;

  function Homotopy_Membership_Test
              ( file : file_type; c : Standard_Irreducible_Component;
                x : Standard_Complex_Vectors.Vector; tol : double_float )
              return boolean is

    sps,sps_last : Standard_Sample_List;
    isin : natural32;

  begin
    Membership_Test(c.p,x,tol,isin,sps,sps_last);
    Deep_Clear(sps);
    if isin > 0 then
      put(file,"Test point corresponds to generic point ");
      put(file,c.g(integer32(isin)),1); put_line(file,".");
      return true;
    else
      put_line(file,"Test point does not correspond to any generic point.");
      return false;
    end if;
  end Homotopy_Membership_Test;

-- DESTRUCTORS :

  procedure Clear ( c : in out Standard_Irreducible_Component ) is
  begin
    Deep_Clear(c.f);
    Clear(c.s);
  end Clear;

  procedure Clear ( c : in out Multprec_Irreducible_Component ) is
  begin
    Deep_Clear(c.f);
    Clear(c.s);
  end Clear;

end Irreducible_Components;
