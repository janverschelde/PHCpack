with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Integer_Pruning_Methods;            use Integer_Pruning_Methods;

procedure ts_tropisms is

-- DESCRIPTION :
--   Tests on drivers to reduce computation of candidate tropisms
--   to mixed volume computation.

  function Append_Slack ( v : Standard_Integer_Vectors.Vector )
                        return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns v with one additional component at the end,
  --   equals to the value of some random number.

    res : Standard_Integer_Vectors.Vector(v'first..v'last+1);
    r : constant integer32 := Standard_Random_Numbers.Random(-99,99);

  begin
    res(v'range) := v;
    res(res'last) := r;
    return res;
  end Append_Slack;

  function Append_Slack ( v : List; s : integer32 := 0 ) return List is

  -- DESCRIPTION :
  --   Returns the vectors in v with one additional component 
  --   added at the end, equal to the value of s.

    res,res_last : List;
    tmp : List := v;
    ls : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Append_Slack(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Append_Slack;

  function Append_Slack ( v : Array_of_Lists ) return Array_of_Lists is

  -- DESCRIPTION :
  --   Adds slack variable equal to zero to all lists in v,
  --   except for the last list, for which slack one is given.

    res : Array_of_Lists(v'range);

  begin
    for i in v'first..v'last-1 loop
      res(i) := Append_Slack(v(i),i);
    end loop;
    res(res'last) := Append_Slack(v(v'last),0);
    return res;
  end Append_Slack;

  function Swap ( v : Standard_Integer_Vectors.Vector )
                return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Swaps the first coordinate to the end.

    res : Standard_Integer_Vectors.Vector(v'range);

  begin
    for i in v'first+1..v'last loop
       res(i-1) := v(i);
    end loop;
    res(res'last) := v(v'first);
    return res;
  end Swap;

  function Back ( v : Standard_Integer_Vectors.Vector )
                return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Opposite of the swap move.

    res : Standard_Integer_Vectors.Vector(v'range);

  begin
    for i in v'first..v'last-1 loop
       res(i+1) := v(i);
    end loop;
    res(res'first) := v(v'last);
    return res;
  end Back;

  function Move ( s : List; swp : boolean ) return List is

  -- DESCRIPTION :
  --   Moves the first coordinate to the last for every point in s.

    res,res_last : List;
    tmp : List := s;
    ls : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if swp
       then Append(res,res_last,Swap(ls.all));
       else Append(res,res_last,Back(ls.all));
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Move;

  function Move ( s : Array_of_Lists; swp : boolean )
                return Array_of_Lists is

  -- DESCRIPTION :
  --   Moves the first coordinate to the last for every point in s.

    res : Array_of_Lists(s'range);

  begin
    for i in res'range loop
       res(i) := Move(s(i),swp);
    end loop;
    return res;
  end Move;

  function Normals ( mcc : Mixed_Subdivision ) return List is

  -- DESCRIPTION :
  --   Returns the inner normals to the cells in mcc.

    res,res_last : List;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Append(res,res_last,mic.nor.all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Normals;

  function Candidate_Tropisms0 ( s : Array_of_Lists ) return List is

  -- DESCRIPTION :
  --   Returns a list of candidate tropisms for the supports s.

    n : constant integer32 := s'last-1;
    fs : Array_of_Lists(s'range);
    mix : Standard_Integer_Vectors.Vector(s'first..s'last-1);
    afa : Array_of_Faces(mix'range);
    nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range);
    mcc : Mixed_Subdivision;
    nrmls,res : List;

  begin
    fs := Move(s,true);
    put_line("The supports after the swap :"); put(fs);
    mix := (mix'range => 1);
    nbsucc := (nbsucc'range => 0.0);
    nbfail := (nbsucc'range => 0.0);
    put_line("Creating lower faces ...");
    for i in afa'range loop
      afa(i) := Create_Lower(mix(i),n+1,fs(i));
      put("  found "); put(Length_Of(afa(i)),1);
      put(" edges on lower hull of support ");
      put(i,1); new_line;
    end loop;
    put_line("pruning for mixed cells ...");
    Create_CS(n,mix,afa,fs(mix'range),nbsucc,nbfail,mcc);
    nrmls := Normals(mcc);
    put_line("The inner normals : "); put(nrmls);
    res := Move(nrmls,False);
    return res;
  end Candidate_Tropisms0;

  function Candidate_Tropisms ( s : Array_of_Lists ) return List is

  -- DESCRIPTION :
  --   Returns a list of candidate tropisms for the supports s.

    n : constant integer32 := s'last;
    s1 : constant Array_of_Lists(s'range) := Append_Slack(s);
    fs : constant Array_of_Lists(s'range) := Move(s1,true);
    mix : Standard_Integer_Vectors.Vector(s'range);
    afa : Array_of_Faces(mix'range);
    nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range);
    mcc : Mixed_Subdivision;
    nrmls,res : List;

  begin
    put_line("The supports after the swap :"); put(fs);
    mix := (mix'range => 1);
    nbsucc := (nbsucc'range => 0.0);
    nbfail := (nbsucc'range => 0.0);
    put_line("Creating lower faces ...");
    for i in afa'range loop
      afa(i) := Create_Lower(mix(i),n+1,fs(i));
      put("  found "); put(Length_Of(afa(i)),1);
      put(" edges on lower hull of support ");
      put(i,1); new_line;
    end loop;
    put_line("pruning for mixed cells ...");
    Create_CS(n,mix,afa,fs,nbsucc,nbfail,mcc);
    put("found "); put(Length_Of(mcc),1); put_line(" mixed cells");
    nrmls := Normals(mcc);
    put_line("The inner normals : "); put(nrmls);
    res := Move(nrmls,false);
    return res;
  end Candidate_Tropisms;

  procedure Candidate_Tropisms ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Computes the supports of p and computes candidate tropisms
  --   via a mixed volume computation.

    sp : constant Array_of_Lists(p'range)
       := Supports_of_Polynomial_Systems.Create(p);
    cdtrp : List;

  begin
    put_line("The supports :"); put(sp);
    cdtrp := Candidate_Tropisms(sp);
  --  cdtrp := Candidate_Tropisms0(sp);
    put_line("The candidate tropisms :"); put(cdtrp);
  end Candidate_Tropisms;

  procedure Main is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Computing candidate tropisms ...");
    new_line;
    get(lp);
    Candidate_Tropisms(lp.all);
  end Main;

begin
  Main;
end ts_tropisms;
