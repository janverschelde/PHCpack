with unchecked_deallocation;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Multprec_Complex_Vectors;
with Continuation_Parameters;
with Sampling_Machine;
with Sampling_Laurent_Machine;

package body Sample_Points is

-- INTERNAL STATE :

  use_laurent : boolean := false;

-- DATA STRUCTURES :

  type Standard_Sample_Rep ( n,k : integer32 ) is record
    point : Standard_Complex_Solutions.Solution(n);
    slices : Standard_Complex_VecVecs.VecVec(1..k);
    refined : Multprec_Sample;
  end record;

  type Multprec_Sample_Rep ( n,k : integer32 ) is record
    point : Multprec_Complex_Solutions.Solution(n);
    slices : Multprec_Complex_VecVecs.VecVec(1..k);
    original : Standard_Sample;
  end record;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean ) is
  begin
    use_laurent := laurent;
  end Set_Polynomial_Type;

-- AUXILIARY :

  function New_Slices ( n,k : integer32 )
                      return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns k random hyperplanes in n-space.

    res : Standard_Complex_VecVecs.VecVec(1..k);

  begin
    for i in 1..k loop
      res(i) := new Standard_Complex_Vectors.Vector'(Random_Vector(0,n));
    end loop;
    return res;
  end New_Slices;

  function Parallel_Slices ( n,k : integer32;
                             hyp : Standard_Complex_VecVecs.VecVec )
                           return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns slices parallel to hyp.  Because of restrictions to the
  --   linear span of a component, the parameters n and k have to be
  --   provided as well.  Mainly n is important to be correct.

    res : Standard_Complex_VecVecs.VecVec(1..k);

  begin
    for i in hyp'range loop
      res(i) := new Standard_Complex_Vectors.Vector(0..n);
      res(i)(0) := Random1;
      for j in 1..n loop
        res(i)(j) := hyp(i)(j);
      end loop;
    end loop;
    return res;
  end Parallel_Slices;

-- CREATORS and COPY :

  function Create ( sol : in Standard_Complex_Solutions.Solution;
                    hyp : in Standard_Complex_VecVecs.VecVec )
                  return Standard_Sample is

    res : Standard_Sample;
    res_rep : Standard_Sample_Rep(sol.n,hyp'last);

  begin
    res_rep.point := sol;
   -- res_rep.slices := hyp;
    for i in res_rep.slices'range loop
      res_rep.slices(i) := new Standard_Complex_Vectors.Vector'(hyp(i).all);
    end loop;
    res := new Standard_Sample_Rep'(res_rep);
    return res;
  end Create;

  function Create ( sol : in Multprec_Complex_Solutions.Solution;
                    hyp : in Multprec_Complex_VecVecs.VecVec )
                  return Multprec_Sample is

    res : Multprec_Sample;
    res_rep : Multprec_Sample_Rep(sol.n,hyp'last);

  begin
    res_rep.point := sol;
    res_rep.slices := hyp;
    res := new Multprec_Sample_Rep'(res_rep);
    return res;
  end Create;

  procedure Copy ( s1 : in Standard_Sample; s2 : out Standard_Sample ) is
  begin
    if s1 /= null then
      declare
        s2sol : constant Standard_Complex_Solutions.Solution(s1.n) := s1.point;
        s2hyp : Standard_Complex_VecVecs.VecVec(1..s1.k);
      begin
        for i in s1.slices'range loop
          s2hyp(i) := new Standard_Complex_Vectors.Vector'(s1.slices(i).all);
        end loop;
        s2 := Create(s2sol,s2hyp);
      end;
    end if;
  end Copy;

  procedure Copy ( s1 : in Multprec_Sample; s2 : out Multprec_Sample ) is
  begin
    if s1 /= null
     then
       declare
         s2sol : Multprec_Complex_Solutions.Solution(s1.n);
         s2hyp : Multprec_Complex_VecVecs.VecVec(1..s1.k);
       begin
         Multprec_Complex_Solutions.Copy(s1.point,s2sol);
         for i in s1.slices'range loop
           s2hyp(i) := new Multprec_Complex_Vectors.Vector(s1.slices(i)'range);
           Multprec_Complex_Vectors.Copy(s1.slices(i).all,s2hyp(i).all);
         end loop;
         s2 := Create(s2sol,s2hyp);
       end;
    end if;
  end Copy;

-- SAMPLERS and REFINERS :

  procedure Sample ( s1 : in Standard_Sample; s2 : out Standard_Sample ) is

    hyp : Standard_Complex_VecVecs.VecVec(1..s1.k);
    sol : Standard_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := New_Slices(s1.n,s1.k);
      if use_laurent then
        Sampling_Laurent_Machine.Sample(s1.point,hyp,sol);
        exit when Sampling_Laurent_Machine.Satisfies(sol);
      else
        Sampling_Machine.Sample(s1.point,hyp,sol);
        exit when Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Sample;

  procedure Sample ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Standard_Sample ) is

    hyp : Standard_Complex_VecVecs.VecVec(1..s1.k);
    sol : Standard_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := New_Slices(s1.n,s1.k);
      if use_laurent then
        Sampling_Laurent_Machine.Sample(file,full_output,s1.point,hyp,sol);
        exit when Sampling_Laurent_Machine.Satisfies(sol);
      else
        Sampling_Machine.Sample(file,full_output,s1.point,hyp,sol);
        exit when Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Sample;

  procedure Sample ( s1 : in Standard_Sample; s2 : out Multprec_Sample ) is

    s2st : Standard_Sample;

  begin
    Sample(s1,s2st);
    Refine(s2st,s2);       -- note: s2.orginal = s2st
  end Sample;

  procedure Sample ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Multprec_Sample ) is

    s2st : Standard_Sample;

  begin
    Sample(file,full_output,s1,s2st);
    Refine(file,full_output,s2st,s2);  -- note: s2.original = s2st
  end Sample;

  procedure Parallel_Sample 
                   ( s1 : in Standard_Sample; s2 : out Standard_Sample ) is

    hyp : Standard_Complex_VecVecs.VecVec(1..s1.k);
    sol : Standard_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := Parallel_Slices(s1.n,s1.k,s1.slices);
      if use_laurent then
        Sampling_Laurent_Machine.Sample(s1.point,hyp,sol);
        exit when Sampling_Laurent_Machine.Satisfies(sol);
      else
        Sampling_Machine.Sample(s1.point,hyp,sol);
        exit when Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Parallel_Sample;

  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Standard_Sample ) is

    hyp : Standard_Complex_VecVecs.VecVec(1..s1.k);
    sol : Standard_Complex_Solutions.Solution(s1.n);
    cnt_reruns : natural32 := 0;

  begin
    loop
      hyp := Parallel_Slices(s1.n,s1.k,s1.slices);
      if use_laurent then
        Sampling_Laurent_Machine.Sample(file,full_output,s1.point,hyp,sol);
        exit when Sampling_Laurent_Machine.Satisfies(sol);
      else
        Sampling_Machine.Sample(file,full_output,s1.point,hyp,sol);
        exit when Sampling_Machine.Satisfies(sol);
      end if;
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
    end loop;
    s2 := Create(sol,hyp);
  end Parallel_Sample;

  procedure Parallel_Sample
                   ( s1 : in Standard_Sample; s2 : out Multprec_Sample ) is

    s2st : Standard_Sample;

  begin
    Parallel_Sample(s1,s2st);
    Refine(s2st,s2);       -- note: s2.orginal = s2st
  end Parallel_Sample;

  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Multprec_Sample ) is

    s2st : Standard_Sample;

  begin
    Parallel_Sample(file,full_output,s1,s2st);
    Refine(file,full_output,s2st,s2);  -- note: s2.original = s2st
  end Parallel_Sample;

  procedure Sample_on_Slices
                   ( s1 : in Standard_Sample;
                     hyp : in Standard_Complex_VecVecs.VecVec;
                     s2 : out Standard_Sample ) is

    sol : Standard_Complex_Solutions.Solution(s1.n);

  begin
    if use_laurent
     then Sampling_Laurent_Machine.Sample(s1.point,hyp,sol);
     else Sampling_Machine.Sample(s1.point,hyp,sol);
    end if;
    s2 := Create(sol,hyp);
  end Sample_on_Slices;

  procedure Sample_on_Slices 
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample;
                     hyp : in Standard_Complex_VecVecs.VecVec;
                     s2 : out Standard_Sample ) is

    sol : Standard_Complex_Solutions.Solution(s1.n);

  begin
    if use_laurent
     then Sampling_Laurent_Machine.Sample(file,full_output,s1.point,hyp,sol);
     else Sampling_Machine.Sample(file,full_output,s1.point,hyp,sol);
    end if;
    s2 := Create(sol,hyp);
  end Sample_on_Slices;

  procedure Sample_on_Slices
                   ( s1 : in Standard_Sample;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample ) is

    s2st : Standard_Sample;

  begin
    Sample_on_Slices(s1,sthyp,s2st);
    Refine_on_Slices(s2st,mphyp,s2);
  end Sample_on_Slices;

  procedure Sample_on_Slices 
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample ) is

    s2st : Standard_Sample;

  begin
    Sample_on_Slices(file,full_output,s1,sthyp,s2st);
    Refine_on_Slices(file,full_output,s2st,mphyp,s2);
  end Sample_on_Slices;

  procedure Refine ( s1 : in out Standard_Sample; s2 : out Multprec_Sample ) is

    hyp : Multprec_Complex_VecVecs.VecVec(1..s1.k);
    sol : Multprec_Complex_Solutions.Solution(s1.n);

  begin
    if use_laurent
     then Sampling_Laurent_Machine.Refine(s1.point,s1.slices,sol,hyp);
     else Sampling_Machine.Refine(s1.point,s1.slices,sol,hyp);
    end if;
    s2 := Create(sol,hyp);
    s2.original := s1;
    s1.refined := s2;
  end Refine;

  procedure Refine ( file : in file_type; full_output : in boolean;
                     s1 : in out Standard_Sample; s2 : out Multprec_Sample ) is

    hyp : Multprec_Complex_VecVecs.VecVec(1..s1.k);
    sol : Multprec_Complex_Solutions.Solution(s1.n);

  begin
    if use_laurent then
      Sampling_Laurent_Machine.Refine
        (file,full_output,s1.point,s1.slices,sol,hyp);
    else
      Sampling_Machine.Refine
        (file,full_output,s1.point,s1.slices,sol,hyp);
    end if;
    s2 := Create(sol,hyp);
    s2.original := s1;
    s1.refined := s2;
  end Refine;

  procedure Refine ( s : in out Multprec_Sample ) is
  begin
    if use_laurent
     then Sampling_Laurent_Machine.Refine(s.point,s.slices);
     else Sampling_Machine.Refine(s.point,s.slices);
    end if;
  end Refine;

  procedure Refine ( file : in file_type; full_output : in boolean;
                     s : in out Multprec_Sample ) is
  begin
    if use_laurent
     then Sampling_Laurent_Machine.Refine(file,full_output,s.point,s.slices);
     else Sampling_Machine.Refine(file,full_output,s.point,s.slices);
    end if;
  end Refine;

  procedure Refine_on_Slices
                   ( s1 : in out Standard_Sample;
                     hyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample ) is

    sol : Multprec_Complex_Solutions.Solution(s1.n);
    s2hyp : Multprec_Complex_VecVecs.VecVec(hyp'range);

  begin
    if use_laurent
     then Sampling_Laurent_Machine.Refine_on_Slices(s1.point,s1.slices,hyp,sol);
     else Sampling_Machine.Refine_on_Slices(s1.point,s1.slices,hyp,sol);
    end if;
    for i in hyp'range loop
      s2hyp(i) := new Multprec_Complex_Vectors.Vector(hyp(i)'range);
      Multprec_Complex_Vectors.Copy(hyp(i).all,s2hyp(i).all);
    end loop;
    s2 := Create(sol,s2hyp);
    s2.original := s1;
    s1.refined := s2;
  end Refine_on_Slices;

  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s1 : in out Standard_Sample;
                     hyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample ) is

    sol : Multprec_Complex_Solutions.Solution(s1.n);
    s2hyp : Multprec_Complex_VecVecs.VecVec(hyp'range);

  begin
    if use_laurent then
      Sampling_Laurent_Machine.Refine_on_Slices
        (file,full_output,s1.point,s1.slices,hyp,sol);
    else
      Sampling_Machine.Refine_on_Slices
        (file,full_output,s1.point,s1.slices,hyp,sol);
    end if;
    for i in hyp'range loop
      s2hyp(i) := new Multprec_Complex_Vectors.Vector(hyp(i)'range);
      Multprec_Complex_Vectors.Copy(hyp(i).all,s2hyp(i).all);
    end loop;
    s2 := Create(sol,s2hyp);
    s2.original := s1;
    s1.refined := s2;
  end Refine_on_Slices;

  procedure Refine_on_Slices ( s : in out Multprec_Sample ) is
  begin
    if use_laurent
     then Sampling_Laurent_Machine.Refine_on_Slices(s.point,s.slices);
     else Sampling_Machine.Refine_on_Slices(s.point,s.slices);
    end if;
  end Refine_on_Slices;

  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s : in out Multprec_Sample ) is
  begin
    if use_laurent then
      Sampling_Laurent_Machine.Refine_on_Slices
        (file,full_output,s.point,s.slices);
    else
      Sampling_Machine.Refine_on_Slices(file,full_output,s.point,s.slices);
    end if;
  end Refine_on_Slices;

-- SELECTORS :

  function Number_of_Variables ( s : Standard_Sample ) return integer32 is
  begin
    if s = null
     then return 0;
     else return s.n;
    end if;
  end Number_of_Variables;

  function Number_of_Variables ( s : Multprec_Sample ) return integer32 is
  begin
    if s = null
     then return 0;
     else return s.n;
    end if;
  end Number_of_Variables;

  function Number_of_Slices ( s : Standard_Sample ) return integer32 is
  begin
    if s = null
     then return 0;
     else return s.k;
    end if;
  end Number_of_Slices;

  function Number_of_Slices ( s : Multprec_Sample ) return integer32 is
  begin
    if s = null
     then return 0;
     else return s.k;
    end if;
  end Number_of_Slices;

  function Sample_Point ( s : Standard_Sample )
                        return Standard_Complex_Solutions.Solution is
  begin
    return s.point;
  end Sample_Point;

  function Sample_Point ( s : Multprec_Sample )
                        return Multprec_Complex_Solutions.Solution is
  begin
    return s.point;
  end Sample_Point;

  function Hyperplane_Sections ( s : Standard_Sample )
                               return Standard_Complex_VecVecs.VecVec is
  begin
    return s.slices;
  end Hyperplane_Sections;

  function Hyperplane_Sections ( s : Multprec_Sample )
                               return Multprec_Complex_VecVecs.VecVec is
  begin
    return s.slices;
  end Hyperplane_Sections;

  function Refined ( s : Standard_Sample ) return Multprec_Sample is
  begin
    return s.refined;
  end Refined;

  function Original ( s : Multprec_Sample ) return Standard_Sample is
  begin
    return s.original;
  end Original;

-- DESTRUCTORS :

  procedure free is
    new unchecked_deallocation(Standard_Sample_Rep,Standard_Sample);
  procedure free is
    new unchecked_deallocation(Multprec_Sample_Rep,Multprec_Sample);

  procedure Shallow_Clear ( s : in out Standard_Sample ) is
  begin
    if s /= null
     then free(s);
    end if;
  end Shallow_Clear;

  procedure Shallow_Clear ( s : in out Multprec_Sample ) is
  begin
    if s /= null
     then free(s);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( s : in out Standard_Sample ) is

    use Standard_Complex_Vectors;

  begin
    if s /= null then
      for i in s.slices'range loop
        Standard_Complex_Vectors.Clear(s.slices(i));
      end loop;
      free(s);
    end if;
  end Deep_Clear;

  procedure Deep_Clear ( s : in out Multprec_Sample ) is
  begin
    if s /= null then
      Multprec_Complex_Solutions.Clear(s.point);
      for i in s.slices'range loop
        Multprec_Complex_Vectors.Clear(s.slices(i));
      end loop;
      free(s);
    end if;
  end Deep_Clear;

begin
  use_laurent := false;
end Sample_Points;
