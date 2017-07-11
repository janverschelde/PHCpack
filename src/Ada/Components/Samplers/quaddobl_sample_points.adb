with unchecked_deallocation;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Numbers;            use QuadDobl_Random_Numbers;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with Continuation_Parameters;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;

package body QuadDobl_Sample_Points is

-- INTERNAL STATE :

  use_laurent : boolean := false;

-- DATA STRUCTURES :

  type QuadDobl_Sample_Rep ( n,k : integer32 ) is record
    point : QuadDobl_Complex_Solutions.Solution(n);
    slices : QuadDobl_Complex_VecVecs.VecVec(1..k);
  end record;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean ) is
  begin
    use_laurent := laurent;
  end Set_Polynomial_Type;

-- AUXILIARY :

  function New_Slices ( n,k : integer32 )
                      return QuadDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns k random hyperplanes in n-space.

    res : QuadDobl_Complex_VecVecs.VecVec(1..k);

  begin
    for i in 1..k loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'(Random_Vector(0,n));
    end loop;
    return res;
  end New_Slices;

  function Parallel_Slices ( n,k : integer32;
                             hyp : QuadDobl_Complex_VecVecs.VecVec )
                           return QuadDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns slices parallel to hyp.  Because of restrictions to the
  --   linear span of a component, the parameters n and k have to be
  --   provided as well.  Mainly n is important to be correct.

    res : QuadDobl_Complex_VecVecs.VecVec(1..k);

  begin
    for i in hyp'range loop
      res(i) := new QuadDobl_Complex_Vectors.Vector(0..n);
      res(i)(0) := Random1;
      for j in 1..n loop
        res(i)(j) := hyp(i)(j);
      end loop;
    end loop;
    return res;
  end Parallel_Slices;

-- CREATORS and COPY :

  function Create ( sol : in QuadDobl_Complex_Solutions.Solution;
                    hyp : in QuadDobl_Complex_VecVecs.VecVec )
                  return QuadDobl_Sample is

    res : QuadDobl_Sample;
    res_rep : QuadDobl_Sample_Rep(sol.n,hyp'last);

  begin
    res_rep.point := sol;
   -- res_rep.slices := hyp;
    for i in res_rep.slices'range loop
      res_rep.slices(i) := new QuadDobl_Complex_Vectors.Vector'(hyp(i).all);
    end loop;
    res := new QuadDobl_Sample_Rep'(res_rep);
    return res;
  end Create;

  procedure Copy ( s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample ) is
  begin
    if s1 /= null then
      declare
        s2sol : constant QuadDobl_Complex_Solutions.Solution(s1.n) := s1.point;
        s2hyp : QuadDobl_Complex_VecVecs.VecVec(1..s1.k);
      begin
        for i in s1.slices'range loop
          s2hyp(i) := new QuadDobl_Complex_Vectors.Vector'(s1.slices(i).all);
        end loop;
        s2 := Create(s2sol,s2hyp);
      end;
    end if;
  end Copy;

-- SAMPLERS :

  procedure Sample ( s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample ) is

    hyp : QuadDobl_Complex_VecVecs.VecVec(1..s1.k);
    sol : QuadDobl_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := New_Slices(s1.n,s1.k);
      if use_laurent then
        QuadDobl_Sampling_Laurent_Machine.Sample(s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Laurent_Machine.Satisfies(sol);
      else
        QuadDobl_Sampling_Machine.Sample(s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Sample;

  procedure Sample ( file : in file_type; full_output : in boolean;
                     s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample ) is

    hyp : QuadDobl_Complex_VecVecs.VecVec(1..s1.k);
    sol : QuadDobl_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := New_Slices(s1.n,s1.k);
      if use_laurent then
        QuadDobl_Sampling_Laurent_Machine.Sample
          (file,full_output,s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Laurent_Machine.Satisfies(sol);
      else
        QuadDobl_Sampling_Machine.Sample(file,full_output,s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Sample;

  procedure Parallel_Sample 
              ( s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample ) is

    hyp : QuadDobl_Complex_VecVecs.VecVec(1..s1.k);
    sol : QuadDobl_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := Parallel_Slices(s1.n,s1.k,s1.slices);
      if use_laurent then
        QuadDobl_Sampling_Laurent_Machine.Sample(s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Laurent_Machine.Satisfies(sol);
      else
        QuadDobl_Sampling_Machine.Sample(s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Parallel_Sample;

  procedure Parallel_Sample
              ( file : in file_type; full_output : in boolean;
                s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample ) is

    hyp : QuadDobl_Complex_VecVecs.VecVec(1..s1.k);
    sol : QuadDobl_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := Parallel_Slices(s1.n,s1.k,s1.slices);
      if use_laurent then
        QuadDobl_Sampling_Laurent_Machine.Sample
          (file,full_output,s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Laurent_Machine.Satisfies(sol);
      else
        QuadDobl_Sampling_Machine.Sample(file,full_output,s1.point,hyp,sol);
        exit when QuadDobl_Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Parallel_Sample;

  procedure Sample_on_Slices
              ( s1 : in QuadDobl_Sample;
                hyp : in QuadDobl_Complex_VecVecs.VecVec;
                s2 : out QuadDobl_Sample ) is

    sol : QuadDobl_Complex_Solutions.Solution(s1.n);

  begin
    if use_laurent
     then QuadDobl_Sampling_Laurent_Machine.Sample(s1.point,hyp,sol);
     else QuadDobl_Sampling_Machine.Sample(s1.point,hyp,sol);
    end if;
    s2 := Create(sol,hyp);
  end Sample_on_Slices;

  procedure Sample_on_Slices 
              ( file : in file_type; full_output : in boolean;
                s1 : in QuadDobl_Sample;
                hyp : in QuadDobl_Complex_VecVecs.VecVec;
                s2 : out QuadDobl_Sample ) is

    sol : QuadDobl_Complex_Solutions.Solution(s1.n);

  begin
    if use_laurent then
      QuadDobl_Sampling_Laurent_Machine.Sample
        (file,full_output,s1.point,hyp,sol);
    else
      QuadDobl_Sampling_Machine.Sample(file,full_output,s1.point,hyp,sol);
    end if;
    s2 := Create(sol,hyp);
  end Sample_on_Slices;

-- SELECTORS :

  function Number_of_Variables ( s : QuadDobl_Sample ) return integer32 is
  begin
    if s = null
     then return 0;
     else return s.n;
    end if;
  end Number_of_Variables;

  function Number_of_Slices ( s : QuadDobl_Sample ) return integer32 is
  begin
    if s = null
     then return 0;
     else return s.k;
    end if;
  end Number_of_Slices;

  function Sample_Point ( s : QuadDobl_Sample )
                        return QuadDobl_Complex_Solutions.Solution is
  begin
    return s.point;
  end Sample_Point;

  function Hyperplane_Sections ( s : QuadDobl_Sample )
                               return QuadDobl_Complex_VecVecs.VecVec is
  begin
    return s.slices;
  end Hyperplane_Sections;

-- DESTRUCTORS :

  procedure free is
    new unchecked_deallocation(QuadDobl_Sample_Rep,QuadDobl_Sample);

  procedure Shallow_Clear ( s : in out QuadDobl_Sample ) is
  begin
    if s /= null
     then free(s);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( s : in out QuadDobl_Sample ) is

    use QuadDobl_Complex_Vectors;

  begin
    if s /= null then
      for i in s.slices'range loop
        QuadDobl_Complex_Vectors.Clear(s.slices(i));
      end loop;
      free(s);
    end if;
  end Deep_Clear;

begin
  use_laurent := false;
end QuadDobl_Sample_Points;
