with Standard_Integer_Vectors;
with Standard_Floating_Matrices;
with Lists_of_Integer_Vectors;
with Floating_Face_Enumerators;          use Floating_Face_Enumerators;
with Facet_Vertex_Enumeration;           use Facet_Vertex_Enumeration;

package body Floating_Faces_of_Polytope is

-- AUXILIAIRIES :

  function Create_Edge ( pts : VecVec; i,j : integer32 ) return Face is

  -- DESCRIPTION :
  --   Creates the edge spanned by pts(i) and pts(j).

    res : constant Face(0..1) := new VecVec(0..1);

  begin
    res(0) := new Vector'(pts(i).all);
    res(1) := new Vector'(pts(j).all);
    return res;
  end Create_Edge;

  function Create_Face ( pts : VecVec;
                         f : Standard_Integer_Vectors.Vector ) return Face is

  -- DESCRIPTION :
  --   Returns vector of points pts(f(i)) that span the face.

    res : constant Face(f'range) := new VecVec(f'range);

  begin
    for i in f'range loop
      res(i) := new Vector'(pts(f(i)).all);
    end loop;
    return res;
  end Create_Face;

  procedure Move_to_Front ( pts : in out VecVec;
                            x : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector x is move to the front of the vector pts.

  begin
    if pts(pts'first).all /= x then
      for i in pts'first+1..pts'last loop
        if pts(i).all = x then
          pts(i).all := pts(pts'first).all;
          pts(pts'first).all := x;
          return;
        end if;
      end loop;
    end if;
  end Move_to_Front;

-- CONSTRUCTORS :

  function Create ( k,n : integer32; p : List; tol : double_float )
                  return Faces is

    res : Faces;

  begin
    if k > n then
      return res;
    else
      declare
        m : constant integer32 := integer32(Length_Of(p));
        pts : constant VecVec(1..m) := Shallow_Create(p);
        res_last : Faces := res;
      begin
        if k = 1 then
          declare
            procedure Append_Edge ( i,j : in integer32; cont : out boolean ) is
              f : constant Face := Create_Edge(pts,i,j);
            begin
              Append(res,res_last,f); cont := true;
            end Append_Edge;
            procedure Enum_Edges is new Enumerate_Edges(Append_Edge);
          begin
            Enum_Edges(pts,tol);
          end;
        else
          declare
            procedure Append_Face ( fa : in Standard_Integer_Vectors.Vector;
                                    cont : out boolean ) is
              f : constant Face := Create_Face(pts,fa);
            begin
              Append(res,res_last,f); cont := true;
            end Append_Face;
            procedure Enum_Faces is new Enumerate_Faces(Append_Face); 
          begin
            Enum_Faces(k,pts,tol);
          end;
        end if;
        return res;
      end;
    end if;
  end Create;

  function Create ( k,n : integer32; p : List; x : Vector; tol : double_float )
                  return Faces is

    res : Faces;

  begin
    if k > n then
      return res;
    else
      declare
        m : constant integer32 := integer32(Length_Of(p));
        pts : VecVec(1..m) := Shallow_Create(p);
        res_last : Faces := res;
      begin
        Move_to_Front(pts,x);
        if k = 1 then
          declare
            procedure Append_Edge ( i,j : in integer32; cont : out boolean ) is
              f : Face;
            begin
              if i = pts'first then
                f := Create_Edge(pts,i,j);
                Append(res,res_last,f);
                cont := true;
              else
                cont := false;
              end if;
            end Append_Edge;
            procedure Enum_Edges is new Enumerate_Edges(Append_Edge);
          begin
            Enum_Edges(pts,tol);
          end;
        else
          declare
            procedure Append_Face ( fa : in Standard_Integer_Vectors.Vector;
                                    cont : out boolean ) is
              f : Face;
            begin
              if fa(fa'first) = pts'first then
                f := Create_Face(pts,fa);
                Append(res,res_last,f);
                cont := true;
              else 
                cont := false;
              end if;
            end Append_Face;
            procedure Enum_Faces is new Enumerate_Faces(Append_Face);
          begin
            Enum_Faces(k,pts,tol);
          end;
        end if;
        return res;
      end;
    end if;
  end Create;

  function Create_Lower ( k,n : integer32; p : List; tol : double_float )
                        return Faces is

    res : Faces;

  begin
    if k > n then
      return res;
    else
      declare
        m : constant integer32 := integer32(Length_Of(p));
        pts : constant VecVec(1..m) := Shallow_Create(p);
        res_last : Faces := res;
      begin
        if k = 1 then
          declare
            procedure Append_Edge ( i,j : in integer32; cont : out boolean ) is
              f : constant Face := Create_Edge(pts,i,j);
            begin
              Append(res,res_last,f); cont := true;
            end Append_Edge;
            procedure Enum_Edges is new Enumerate_Lower_Edges(Append_Edge);
          begin
            Enum_Edges(pts,tol);
          end;
        else
          declare
            procedure Append_Face ( fa : in Standard_Integer_Vectors.Vector;
                                    cont : out boolean ) is
              f : constant Face := Create_Face(pts,fa);
            begin
              Append(res,res_last,f); cont := true;
            end Append_Face;
            procedure Enum_Faces is new Enumerate_Lower_Faces(Append_Face);
          begin
            Enum_Faces(k,pts,tol);
          end;
        end if;
        return res;
      end;
    end if;
  end Create_Lower;

  function Create_Lower ( k,n : integer32; p : List; x : Vector;
                          tol : double_float ) return Faces is

    res : Faces;

  begin
    if k > n then
      return res;
    else
      declare
        m : constant integer32 := integer32(Length_Of(p));
        pts : VecVec(1..m) := Shallow_Create(p);
        res_last : Faces := res;
      begin
        Move_to_Front(pts,x);
        if k = 1 then
          declare
            procedure Append_Edge ( i,j : in integer32; cont : out boolean ) is
              f : Face := Create_Edge(pts,i,j);
            begin
              if i = pts'first then
                f := Create_Edge(pts,i,j);
                Append(res,res_last,f);
                cont := true;
              else
                cont := false;
              end if;
            end Append_Edge;
            procedure Enum_Edges is new Enumerate_Lower_Edges(Append_Edge);
          begin
            Enum_Edges(pts,tol);
          end;
        else
          declare
            procedure Append_Face ( fa : in Standard_Integer_Vectors.Vector;
                                    cont : out boolean ) is
              f : Face;
            begin
              if fa(fa'first) = pts'first then
                f := Create_Face(pts,fa);
                Append(res,res_last,f);
                cont := true;
              else
                cont := false;
              end if;
            end Append_Face;
            procedure Enum_Faces is new Enumerate_Lower_Faces(Append_Face);
          begin
            Enum_Faces(k,pts,tol);
          end;
        end if;
        return res;
      end;
    end if;
  end Create_Lower;

  function Create_Lower_Facets
             ( n : integer32; p : List; tol : double_float ) return Faces is

    res,res_last : Faces;
    m : constant integer32 := integer32(Length_Of(p));
    pts : constant VecVec(1..m) := Shallow_Create(p); 
    ptsmat : Standard_Floating_Matrices.Matrix(1..n,1..m);
    facets,tmp : Lists_of_Integer_Vectors.List;

  begin
    for i in 1..n loop
      for j in 1..m loop
        ptsmat(i,j) := pts(j)(i);
      end loop;
    end loop;
    facets := Enumerate_Lower_Facet_Labels(ptsmat,tol);
    tmp := facets;
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      declare
        f : constant Face
          := Create_Face(pts,Lists_of_Integer_Vectors.Head_Of(tmp).all);
      begin
        Append(res,res_last,f);
      end;
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    Lists_of_Integer_Vectors.Clear(facets);
    return res;
  end Create_Lower_Facets;

  procedure Construct ( first : in out Faces; fs : in Faces ) is

    tmp : Faces := fs;

  begin
    while not Is_Null(tmp) loop
      Construct(Head_Of(tmp),first);
      tmp := Tail_Of(tmp);
    end loop;
  end Construct;

-- SELECTORS :

  function Is_Equal ( f1,f2 : Face ) return boolean is

    found : boolean;

  begin
    for i in f1'range loop
      found := false;
      for j in f2'range loop
        found := Equal(f1(i).all,f2(j).all);
        exit when found;
      end loop;
      if not found
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Is_In ( f : Face; x : Vector ) return boolean is
  begin
    for i in f'range loop
      if f(i).all = x
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( fs : Faces; f : Face ) return boolean is

    tmp : Faces := fs;

  begin
    while not Is_Null(tmp) loop
      if Is_Equal(f,Head_Of(tmp))
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

-- DESTRUCTORS :

  procedure Deep_Clear ( f : in out Face ) is
  begin
    if f /= null
     then for i in f'range loop
            Clear(f(i));
          end loop;
    end if;
  end Deep_Clear;

  procedure Shallow_Clear ( f : in out Face ) is
  begin
    if f /= null
     then Clear(f.all);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( fa : in out Face_Array ) is
  begin
    for i in fa'range loop
      Deep_Clear(fa(i));
    end loop;
  end Deep_Clear;

  procedure Shallow_Clear ( fa : in out Face_Array ) is
  begin
    for i in fa'range loop
      Shallow_Clear(fa(i));
    end loop;
  end Shallow_Clear;

  procedure Deep_Clear ( fs : in out Faces ) is

    tmp : Faces := fs;

  begin
    while not Is_Null(tmp) loop
      declare
	f : Face := Head_Of(tmp);
      begin
	Deep_Clear(f);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Faces.Clear(Lists_of_Faces.List(fs));
  end Deep_Clear;

  procedure Shallow_Clear ( fs : in out Faces ) is
  begin
    Lists_of_Faces.Clear(Lists_of_Faces.List(fs));
  end Shallow_Clear;

  procedure Deep_Clear ( afs : in out Array_of_Faces ) is
  begin
    for i in afs'range loop
      Deep_Clear(afs(i));
    end loop;
  end Deep_Clear;

  procedure Shallow_Clear ( afs : in out Array_of_Faces ) is
  begin
    for i in afs'range loop
      Shallow_Clear(afs(i));
    end loop;
  end Shallow_Clear;

end Floating_Faces_of_Polytope;
