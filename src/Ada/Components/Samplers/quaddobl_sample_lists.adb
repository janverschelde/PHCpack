with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;

package body QuadDobl_Sample_Lists is

-- INTERNAL STATE :

  use_laurent : boolean := false;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean ) is
  begin
    use_laurent := laurent;
  end Set_Polynomial_Type;

-- CREATORS :

  function Create ( sols : QuadDobl_Complex_Solutions.Solution_List;
                    hyps : QuadDobl_Complex_VecVecs.VecVec )
                  return QuadDobl_Sample_List is

    res,res_last : QuadDobl_Sample_List;
    use QuadDobl_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        s : constant QuadDobl_Sample := Create(Head_Of(tmp).all,hyps);
      begin
        Append(res,res_last,s);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

-- SAMPLERS :

  procedure Sample ( s : in QuadDobl_Sample; m : in natural32;
                     samples,samples_last : in out QuadDobl_Sample_List ) is
  begin
    for i in 1..m loop
      declare
        s2 : QuadDobl_Sample;
      begin
        Sample(s,s2);
        Append(samples,samples_last,s2);
      end;
    end loop;
  end Sample;

  procedure Sample ( file : in file_type; full_output : in boolean;
                     s : in QuadDobl_Sample; m : in natural32;
                     samples,samples_last : in out QuadDobl_Sample_List ) is
  begin
    for i in 1..m loop
      declare
        s2 : QuadDobl_Sample;
      begin
        Sample(file,full_output,s,s2);
        Append(samples,samples_last,s2);
      end;
    end loop;
  end Sample;

  procedure Parallel_Sample
              ( s : in QuadDobl_Sample; m : in natural32;
                samples,samples_last : in out QuadDobl_Sample_List ) is
  begin
    for i in 1..m loop
      declare
        s2 : QuadDobl_Sample;
      begin
        Parallel_Sample(s,s2);
        Append(samples,samples_last,s2);
      end;
    end loop;
  end Parallel_Sample;

  procedure Parallel_Sample
              ( file : in file_type; full_output : in boolean;
                s : in QuadDobl_Sample; m : in natural32;
                samples,samples_last : in out QuadDobl_Sample_List ) is
  begin
    for i in 1..m loop
      declare
        s2 : QuadDobl_Sample;
      begin
        Parallel_Sample(file,full_output,s,s2);
        Append(samples,samples_last,s2);
      end;
    end loop;
  end Parallel_Sample;

  procedure Create_Samples
              ( s2sols : in QuadDobl_Complex_Solutions.Solution_List;
                hyps : in QuadDobl_Complex_VecVecs.VecVec;
                s2,s2_last : in out QuadDobl_Sample_List ) is

  -- DESCRIPTION :
  --   This routine is auxiliary, it creates samples from a list of
  --   solutions and a set of slices.

    use QuadDobl_Complex_Solutions;
    tmp : Solution_List;

  begin
    tmp := s2sols;
    while not Is_Null(tmp) loop
      declare
        newspt : constant QuadDobl_Sample := Create(Head_Of(tmp).all,hyps);
      begin
        Append(s2,s2_last,newspt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Create_Samples;

  procedure Sample ( s1 : in QuadDobl_Sample_List;
                     hyps : in QuadDobl_Complex_VecVecs.VecVec;
                     s2,s2_last : in out QuadDobl_Sample_List ) is

    s1sols : QuadDobl_Complex_Solutions.Solution_List := Sample_Points(s1);
    s2sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if use_laurent
     then QuadDobl_Sampling_Laurent_Machine.Sample(s1sols,hyps,s2sols);
     else QuadDobl_Sampling_Machine.Sample(s1sols,hyps,s2sols);
    end if;
    Create_Samples(s2sols,hyps,s2,s2_last);
    QuadDobl_Complex_Solutions.Deep_Clear(s1sols);
    QuadDobl_Complex_Solutions.Deep_Clear(s2sols);
  end Sample;

  procedure Sample ( file : in file_type;
                     s1 : in QuadDobl_Sample_List;
                     hyps : in QuadDobl_Complex_VecVecs.VecVec;
                     s2,s2_last : in out QuadDobl_Sample_List ) is

    s1sols : QuadDobl_Complex_Solutions.Solution_List := Sample_Points(s1);
    s2sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if use_laurent
     then QuadDobl_Sampling_Laurent_Machine.Sample(file,s1sols,hyps,s2sols);
     else QuadDobl_Sampling_Machine.Sample(file,s1sols,hyps,s2sols);
    end if;
    Create_Samples(s2sols,hyps,s2,s2_last);
    QuadDobl_Complex_Solutions.Deep_Clear(s1sols);
    QuadDobl_Complex_Solutions.Deep_Clear(s2sols);
  end Sample;

  procedure Sample_with_Stop
              ( s1 : in QuadDobl_Sample_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                tol : in double_float;
                hyps : in QuadDobl_Complex_VecVecs.VecVec;
                s2,s2_last : in out QuadDobl_Sample_List ) is

    s1sols : QuadDobl_Complex_Solutions.Solution_List := Sample_Points(s1);
    use QuadDobl_Complex_Solutions;
    s2sols,tmp : Solution_List;

  begin
    if use_laurent then
      QuadDobl_Sampling_Laurent_Machine.Sample_with_Stop
        (s1sols,x,tol,hyps,s2sols);
    else
      QuadDobl_Sampling_Machine.Sample_with_Stop(s1sols,x,tol,hyps,s2sols);
    end if;
    tmp := s2sols;
    while not Is_Null(tmp) loop
      declare
        newspt : constant QuadDobl_Sample := Create(Head_Of(tmp).all,hyps);
      begin
        Append(s2,s2_last,newspt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Deep_Clear(s1sols);
    Deep_Clear(s2sols);
  end Sample_with_Stop;

-- SAMPLING MEMBERSHIP TEST :

  function Equal ( x,y : QuadDobl_Complex_Vectors.Vector;
                   tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all components of x and y are as close
  --   to each other as the given tolerance.

    diff : Complex_Number;
    absdif : quad_double;

  begin
    for i in x'range loop
      diff := x(i) - y(i);
      absdif := AbsVal(diff);
      if absdif > tol
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  function Is_In ( s : QuadDobl_Sample_List; tol : double_float;
                   x : QuadDobl_Complex_Vectors.Vector ) return natural32 is

    tmp : QuadDobl_Sample_List := s;

  begin
    for i in 1..Length_Of(s) loop
      declare
        y : constant QuadDobl_Complex_Vectors.Vector
          := Sample_Point(Head_Of(tmp)).v;
      begin
        if Equal(x,y,tol)
         then return natural32(i);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return 0;
  end Is_In;

  function Is_In ( s : QuadDobl_Sample_List; tol : double_float;
                   spt : QuadDobl_Sample ) return natural32 is
  begin
    return Is_In(s,tol,Sample_Point(spt).v);
  end Is_In;

  procedure Membership_Test
                ( s1 : in QuadDobl_Sample_List;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  tol : in double_float; isin : out natural32;
                  s2,s2_last : in out QuadDobl_Sample_List ) is

    n : constant integer32 := x'last;
    m : constant integer32 := integer32(Number_of_Slices(Head_Of(s1)));
    nor : QuadDobl_Complex_Vectors.Vector(1..n);
    hyps : QuadDobl_Complex_VecVecs.VecVec(1..m);
    use QuadDobl_Complex_Vectors;
    eva : Complex_Number;

  begin
    for i in 1..m loop                             -- hyperplanes through x
      nor := Random_Vector(1,n);
      eva := x*nor;
      hyps(i) := new Vector(0..n);
      hyps(i)(nor'range) := nor;
      hyps(i)(0) := -eva;
    end loop;
   -- Sample(s1,hyps,s2,s2_last);                    -- conduct samples+test
    Sample_with_Stop(s1,x,tol,hyps,s2,s2_last);
    isin := Is_In(s2,tol,x);
    QuadDobl_Complex_VecVecs.Clear(hyps);
  end Membership_Test;

-- SELECTORS :

  function Sample_Points ( s : QuadDobl_Sample_List )
                         return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : QuadDobl_Sample_List := s;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Sample_Point(Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Sample_Points;

  function Select_Sub_List ( l : QuadDobl_Sample_List;
                             f : Standard_Natural_Vectors.Vector )
                           return QuadDobl_Sample_List is

    res,res_last : QuadDobl_Sample_List;
    tmp : QuadDobl_Sample_List := l;
    ind : integer32 := f'first;
  
  begin
    for i in 1..Length_Of(l) loop
      if i = f(ind) then
        Append(res,res_last,Head_Of(tmp));
        ind := ind + 1;
      end if;
      exit when (ind > f'last);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Select_Sub_List;

  function Select_Sub_Grid ( grid : Array_of_QuadDobl_Sample_Lists;
                             f : Standard_Natural_Vectors.Vector )
                           return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(grid'range);
 
  begin
    for i in grid'range loop
      res(i) := Select_Sub_List(grid(i),f);
    end loop;
    return res;
  end Select_Sub_Grid;

-- DESTRUCTORS :

  procedure Shallow_Clear ( samples : in out QuadDobl_Sample_List ) is
  begin
    Lists_of_QuadDobl_Samples.Clear(Lists_of_QuadDobl_Samples.List(samples));
  end Shallow_Clear;

  procedure Shallow_Clear
                 ( samples : in out Array_of_QuadDobl_Sample_Lists ) is
  begin
    for i in samples'range loop
      Shallow_Clear(samples(i));
    end loop;
  end Shallow_Clear;

  procedure Deep_Clear ( samples : in out QuadDobl_Sample_List ) is

    tmp : QuadDobl_Sample_List := samples;

  begin
    while not Is_Null(tmp) loop
      declare
        sample : QuadDobl_Sample := Head_Of(tmp);
      begin
        Deep_Clear(sample);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(samples);  
  end Deep_Clear;

  procedure Deep_Clear ( samples : in out Array_of_QuadDobl_Sample_Lists ) is
  begin
    for i in samples'range loop
      Deep_Clear(samples(i));
    end loop;
  end Deep_Clear;

begin
  use_laurent := false;
end QuadDobl_Sample_Lists;
