with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Standard_Condition_Tables;

package body Multprec_Condition_Tables is

  function Create ( n : natural32 ) return Vector is

    res : constant Vector(0..integer32(n)) := (0..integer32(n) => 0);

  begin
    return res;
  end Create;

  function Truncate ( f : Floating_Number ) return integer32 is

  -- DESCRIPTION :
  --   Truncates a positive floating-point number to the largest
  --   natural number less than f.

    df : constant double_float := Round(f);
    res : integer32 := integer32(df);

  begin
    if double_float(res) > df
     then res := res - 1;
    end if;
    return res;
  end Truncate;

  procedure Update_Corrector ( t : in out Vector; e : in Floating_Number ) is

    tol : constant double_float := 10.0**(-integer(t'last+1));
    lco : Floating_Number;
    ind : integer32;

  begin
    if e > 1.0 then
      t(0) := t(0) + 1;
    else
      if e < tol then
        t(t'last) := t(t'last) + 1;
      else
        lco := Log10(e);
        Min(lco);
        ind := Truncate(lco);
        if ind < t'first then
          t(t'first) := t(t'first) + 1;
        elsif ind > t'last then
          t(t'last) := t(t'last) + 1;
        else
          t(ind) := t(ind) + 1;
        end if;
        Clear(lco);
      end if;
    end if;
  end Update_Corrector;

  procedure Update_Corrector ( t : in out Vector; s : in Solution ) is
  begin
    Update_Corrector(t,s.err);
  end Update_Corrector;

  procedure Update_Condition ( t : in out Vector; c : in Floating_Number ) is

    tol : constant double_float := 10.0**(-integer(t'last));
    lco : Floating_Number;
    ind : integer32;

  begin
    if c < tol then
      t(t'last) := t(t'last) + 1;
    else
      lco := Log10(c);
      Min(lco);
      ind := Truncate(lco);
      if ind < t'first then
        t(t'first) := t(t'first) + 1;
      elsif ind > t'last then
        t(t'last) := t(t'last) + 1;
      else
        t(ind) := t(ind) + 1;
      end if;
      Clear(lco);
    end if;
  end Update_Condition;

  procedure Update_Condition ( t : in out Vector; s : in Solution ) is
  begin
    Update_Condition(t,s.rco);
  end Update_Condition;

  procedure Update_Distance ( t : in out Vector; d : in Floating_Number ) is

    tol : constant double_float := 10.0**(-integer(t'last));
    lco : Floating_Number;
    ind : integer32;

  begin
    if d < tol then
      t(t'last) := t(t'last) + 1;
    else
      lco := Log10(d);
      Min(lco);
      if lco < 0.0
       then ind := 0;
       else ind := Truncate(lco);
      end if;
      Clear(lco);
      if ind < t'first then
        t(t'first) := t(t'first) + 1;
      elsif ind > t'last then
        t(t'last) := t(t'last) + 1;
      else
        t(ind) := t(ind) + 1;
      end if;
    end if;
  end Update_Distance;

  procedure Update_Residuals ( t : in out Vector; r : in Floating_Number ) is

    tol : constant double_float := 10.0**(-integer(t'last)+1);
    lco : Floating_Number;
    ind : integer32;

  begin
    if r > 1.0 then
      t(0) := t(0) + 1;
    else
      if r < tol then
        t(t'last) := t(t'last) + 1;
      else
        lco := Log10(r);
        Min(lco);
        ind := Truncate(lco);
        if ind < t'first then
          t(t'first) := t(t'first) + 1;
        elsif ind > t'last then
          t(t'last) := t(t'last) + 1;
        else
          t(ind) := t(ind) + 1;
        end if;
        Clear(lco);
      end if;
    end if;
  end Update_Residuals;

  procedure Update_Residuals ( t : in out Vector; s : in Solution ) is
  begin
    Update_Residuals(t,s.res);
  end Update_Residuals;

  procedure Corrector_Table ( t : in out Vector; s : in Solution_List ) is

    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Update_Corrector(t,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Corrector_Table;

  procedure Condition_Table ( t : in out Vector; s : in Solution_List ) is

    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Update_Condition(t,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Condition_Table;

  procedure Distances_Table ( t : in out Vector; s : in Solution_List ) is

    tmp1 : Solution_List := s;
    tmp2 : Solution_List;
    ls1,ls2 : Link_to_Solution;
    d,d12 : Floating_Number;

    use Multprec_Complex_Vectors;

  begin
    for i in 1..Length_Of(s) loop
      ls1 := Head_Of(tmp1);
      tmp2 := s;
      d := Create(1.0E+16);
      for j in 1..Length_Of(s) loop
        if j /= i then
          ls2 := Head_Of(tmp2);
          declare
            wrk : Multprec_Complex_Vectors.Vector(ls1.v'range)
                := ls2.v - ls1.v;
          begin
            d12 := Norm2(wrk);
            Clear(wrk);
          end;
           if d12 < d
            then Copy(d12,d);
           end if;
           Clear(d12);
        end if;
        tmp2 := Tail_Of(tmp2);
      end loop;
      Update_Distance(t,d);
      Clear(d);
      tmp1 := Tail_Of(tmp1);
    end loop;
  end Distances_Table;

  procedure Residuals_Table ( t : in out Vector; s : in Solution_List ) is

    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Update_Residuals(t,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Residuals_Table;

  procedure Write_Corrector_Table ( file : in file_type; t : in Vector ) is

    s : constant natural32 := Sum(t);
    i : constant integer32
      := Standard_Condition_Tables.First_Index_of_Nonzero(t);

  begin
    put_line(file,"Frequency table of logarithms of corrector :");
    put(file,"  FreqCorr : "); put(file,t);
    put(file," : "); put(file,s,1); new_line(file);
    put(file,"  ");
    put(file,"k-th number counts #solutions with k-1 <= -log10(corrector) < k");
    new_line(file);
    put(file,"  the largest corrector error is smaller than 1.0E-");
    put(file,i,1); put_line(file,".");
  end Write_Corrector_Table;

  procedure Write_Condition_Table ( file : in file_type; t : in Vector ) is

    s : constant natural32 := Sum(t);
    i : constant integer32
      := Standard_Condition_Tables.Last_Index_of_Nonzero(t);

  begin
    put_line(file,"Frequency table of logarithms of condition numbers :");
    put(file,"  FreqCond : "); put(file,t);
    put(file," : "); put(file,s,1); new_line(file);
    put(file,"  ");
    put(file,"k-th number counts #solutions with k-1 <= log10(cond#) < k");
    new_line(file);
    put(file,"  the largest condition number is about 1.0E+");
    put(file,i+1,1); put_line(file,".");
  end Write_Condition_Table;

  procedure Write_Distances_Table ( file : in file_type; t : in Vector ) is

    s : constant natural32 := Sum(t);

  begin
    put_line(file,"Frequency table of logarithms of distances to roots :");
    put(file,"FreqDist : "); put(file,t);
    put(file," : "); put(file,s,1); new_line(file);
    put(file,"k-th number counts #solutions with k-1 <= -log10(distance) < k");
    new_line(file);
  end Write_Distances_Table;

  procedure Write_Residuals_Table ( file : in file_type; t : in Vector ) is

    s : constant natural32 := Sum(t);
    i : constant integer32
      := Standard_Condition_Tables.First_Index_of_Nonzero(t);

  begin
    put_line(file,"Frequency table of absolute logarithms of residuals :");
    put(file,"  FreqResi : "); put(file,t);
    put(file," : "); put(file,s,1);  new_line(file);
    put(file,"  ");
    put(file,"k-th number counts #solutions with k-1 <= -log10(residual) < k");
    new_line(file);
    put(file,"  the largest residual is smaller than 1.0E-");
    put(file,i,1); put_line(file,".");
  end Write_Residuals_Table;

  procedure Write_Tables ( file : in file_type;
                           t_e,t_r,t_c : in Vector ) is
   
    s : natural32;

  begin
    put(file,"Frequency tables for correction, residual,");
    put_line(file," and condition numbers :");
    put(file,"FreqCorr : "); put(file,t_e); s := Sum(t_e);
    put(file," : "); put(file,s,1);  new_line(file);
    put(file,"FreqResi : "); put(file,t_r); s := Sum(t_r);
    put(file," : "); put(file,s,1);  new_line(file);
    put(file,"FreqCond : "); put(file,t_c); s := Sum(t_c);
    put(file," : "); put(file,s,1);  new_line(file);
    put_line(file,"Small correction terms and residuals"
                & " counted to the right.");
    put_line(file,"Well conditioned and distinct roots"
                & " counted to the left.");
  end Write_Tables;

  procedure Write_Tables ( file : in file_type;
                           t_e,t_r,t_c,t_d : in Vector ) is
   
    s : natural32;

  begin
    put(file,"Frequency tables for correction, residual, condition,");
    put_line(file," and distances :");
    put(file,"FreqCorr : "); put(file,t_e); s := Sum(t_e);
    put(file," : "); put(file,s,1);  new_line(file);
    put(file,"FreqResi : "); put(file,t_r); s := Sum(t_r);
    put(file," : "); put(file,s,1);  new_line(file);
    put(file,"FreqCond : "); put(file,t_c); s := Sum(t_c);
    put(file," : "); put(file,s,1);  new_line(file);
    put(file,"FreqDist : "); put(file,t_d); s := Sum(t_d);
    put(file," : "); put(file,s,1);  new_line(file);
    put_line(file,"Small correction terms and residuals"
                & " counted to the right.");
    put_line(file,"Well conditioned and distinct roots"
                & " counted to the left.");
  end Write_Tables;

end Multprec_Condition_Tables;
