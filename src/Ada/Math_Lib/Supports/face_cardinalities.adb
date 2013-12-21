with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Integer_Face_Enumerators;
with Floating_Face_Enumerators;

package body Face_Cardinalities is

  function fvector ( pts : Standard_Integer_VecVecs.VecVec ) return Vector is

  -- ALGORITHM : plain enumeration of vertices, edges, k-faces...

    use Integer_Face_Enumerators;

    n : constant integer32 := pts(pts'first).all'length;
    f : Vector(-1..n);
 
    procedure Count_Vertex ( i : in integer32; cont : out boolean ) is
    begin
      f(0) := f(0) + 1;
      cont := true;
    end Count_Vertex;
    procedure Count_Vertices is new Enumerate_Vertices(Count_Vertex);

    procedure Count_Edge ( i,j : in integer32; cont : out boolean ) is
    begin
      f(1) := f(1) + 1;
      cont := true;
    end Count_Edge;
    procedure Count_Edges is new Enumerate_Edges(Count_Edge);

    procedure Count_Face ( face : in Vector; cont : out boolean ) is
    begin
      f(face'length-1) := f(face'length-1) + 1;
      cont := true;
    end Count_Face;
    procedure Count_Faces is new Enumerate_Faces(Count_Face);

  begin
    f(-1) := 1;
    f(0..n) := (0..n => 0);
    Count_Vertices(pts);
    Count_Edges(pts);
    for i in 2..(n-1) loop
      Count_Faces(i,pts);
      exit when (f(i) = 0);
    end loop;
    if f(n-1) > 1
     then f(n) := 1;
    end if;
    return f;
  end fvector;

  function fvector ( pts : Standard_Floating_VecVecs.VecVec ) return Vector is

  -- ALGORITHM : plain enumeration of vertices, edges, k-faces...

    use Floating_Face_Enumerators;

    n : constant integer32 := pts(pts'first).all'length;
    f : Vector(-1..n);
    tol : constant double_float := 10.0**(-8); --10.0**(-12);
 
    procedure Count_Vertex ( i : in integer32; cont : out boolean ) is
    begin
      f(0) := f(0) + 1;
      cont := true;
    end Count_Vertex;
    procedure Count_Vertices is new Enumerate_Vertices(Count_Vertex);

    procedure Count_Edge ( i,j : in integer32; cont : out boolean ) is
    begin
      f(1) := f(1) + 1;
      cont := true;
    end Count_Edge;
    procedure Count_Edges is new Enumerate_Edges(Count_Edge);

    procedure Count_Face ( face : in Vector; cont : out boolean ) is
    begin
      f(face'length-1) := f(face'length-1) + 1;
      cont := true;
    end Count_Face;
    procedure Count_Faces is new Enumerate_Faces(Count_Face);

  begin
    f(-1) := 1;
    f(0..n) := (0..n => 0);
    Count_Vertices(pts,tol);
    Count_Edges(pts,tol);
    for i in 2..(n-1) loop
      Count_Faces(i,pts,tol);
      exit when (f(i) = 0);
    end loop;
    if f(n-1) > 1
     then f(n) := 1;
    end if;
    return f;
  end fvector;

  function Face_Labels ( pts : Standard_Integer_VecVecs.VecVec;
                         k : natural32 ) return List is

    use Integer_Face_Enumerators;

    res,res_last : List;
 
    procedure Add_Vertex ( i : in integer32; cont : out boolean ) is
  
      face : Vector(1..1);

    begin
      face(1) := i;
      Append(res,res_last,face);
      cont := true;
    end Add_Vertex;
    procedure Add_Vertices is new Enumerate_Vertices(Add_Vertex);

    procedure Add_Edge ( i,j : in integer32; cont : out boolean ) is

      face : Vector(1..2);

    begin
      face(1) := i;
      face(2) := j;
      Append(res,res_last,face);
      cont := true;
    end Add_Edge;
    procedure Add_Edges is new Enumerate_Edges(Add_Edge);

    procedure Add_Face ( face : in Vector; cont : out boolean ) is
    begin
      Append(res,res_last,face);
      cont := true;
    end Add_Face;
    procedure Add_Faces is new Enumerate_Faces(Add_Face);

  begin
    if k = 0 then
      Add_Vertices(pts);
    elsif k = 1 then
      Add_Edges(pts);
    else
      Add_Faces(integer32(k),pts);
    end if;
    return res;
  end Face_Labels;

  function Face_Labels ( pts : Standard_Floating_VecVecs.VecVec;
                         k : natural32 ) return List is

    use Floating_Face_Enumerators;

    res,res_last : List;
    tol : constant double_float := 10.0**(-8); --10.0**(-12);
 
    procedure Add_Vertex ( i : in integer32; cont : out boolean ) is

      face : Vector(1..1);

    begin
      face(1) := i;
      Append(res,res_last,face);
      cont := true;
    end Add_Vertex;
    procedure Add_Vertices is new Enumerate_Vertices(Add_Vertex);

    procedure Add_Edge ( i,j : in integer32; cont : out boolean ) is

      face : Vector(1..2);

    begin
      face(1) := i;
      face(2) := j;
      Append(res,res_last,face);
      cont := true;
    end Add_Edge;
    procedure Add_Edges is new Enumerate_Edges(Add_Edge);

    procedure Add_Face ( face : in Vector; cont : out boolean ) is
    begin
      Append(res,res_last,face);
      cont := true;
    end Add_Face;
    procedure Add_Faces is new Enumerate_Faces(Add_Face);

  begin
    if k = 0 then
      Add_Vertices(pts,tol);
    elsif k = 1 then
      Add_Edges(pts,tol);
    else
      Add_Faces(integer32(k),pts,tol);
    end if;
    return res;
  end Face_Labels; 

  function Face_Labels ( pts : Standard_Integer_VecVecs.VecVec )
                       return Array_of_Lists is

    n : constant integer32 := pts(pts'first).all'length;
    res : Array_of_Lists(0..n);

  begin
    for i in res'range loop
      res(i) := Face_Labels(pts,natural32(i));
    end loop;
    return res;
  end Face_Labels;

  function Face_Labels ( pts : Standard_Floating_VecVecs.VecVec )
                       return Array_of_Lists is

    n : constant integer32 := pts(pts'first).all'length;
    res : Array_of_Lists(0..n);

  begin
    for i in res'range loop
      res(i) := Face_Labels(pts,natural32(i));
    end loop;
    return res;
  end Face_Labels;

end Face_Cardinalities;
