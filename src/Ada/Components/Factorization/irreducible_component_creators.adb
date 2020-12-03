with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Multprec_Floating_Matrices;         use Multprec_Floating_Matrices;
with Standard_Complex_Poly_Systems;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_to_Standard_Convertors;    use Multprec_to_Standard_Convertors;
with Multprec_Complex_Poly_Systems;
with Standard_Polynomial_Interpolators;
with Multprec_Polynomial_Interpolators;
with Sampling_Machine;
with Standard_Subspace_Restrictions;     use Standard_Subspace_Restrictions;
with Multprec_Subspace_Restrictions;     use Multprec_Subspace_Restrictions;
with span_of_Component_io;               use Span_of_Component_io;
with Span_of_Component_Creators;         use Span_of_Component_Creators;
with Standard_Stacked_Sample_Grids;
with Multprec_Stacked_Sample_Grids;
with Make_Sample_Grids;                  use Make_Sample_Grids;

package body Irreducible_Component_Creators is

-- AUXILIARY RESTRICTS TO THE SPAN OF COMPONENT :

  procedure Standard_Restrict
              ( file : in file_type; s : in Standard_Span;
                dim : in natural32; tol : in double_float;
                sps,sps_last : in out Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the restriction of the embedded polynomial system and
  --   performs the Initialization in the Sampling Machine.
  --   Also the list of samples gets restricted to the span.

    ep : constant Standard_Complex_Poly_Systems.Poly_Sys
       := Sampling_Machine.Embedded_System;
    rp : constant Standard_Complex_Poly_Systems.Poly_Sys(ep'range)
       := Restrict(s,dim,tol,ep);
    cp : constant Standard_Complex_Poly_Systems.Poly_Sys
       := Collapse_Equations(rp,Dimension(s),dim);
    res,res_last : Standard_Sample_List;
    tmp : Standard_Sample_List := sps;

  begin
   -- put_line(file,"The restricted embedded polynomial system : "); 
   -- put_line(file,rp);
   -- put_line(file,"The collapsed restricted system : ");
   -- put_line(file,cp);
    Sampling_Machine.Initialize_Restricted(cp);
    while not Is_Null(tmp) loop
      declare
        res_spt : constant Standard_Sample := Restrict(s,dim,Head_Of(tmp)); 
       -- res_sol : Standard_Complex_Solutions.Solution(cp'last)
       --         := Sample_Point(res_spt);
       -- eva : Standard_Complex_Vectors.Vector := Eval(cp,res_sol.v);
      begin
        Append(res,res_last,res_spt);
       -- put_line(file,"Evaluation of restrictions : ");
       -- put_line(file,eva);
      end;
      tmp := Tail_Of(tmp);
    end loop;
   -- Deep_Clear(sps);
    sps := res;
    sps_last := res_last;
  end Standard_Restrict;

  procedure Multprec_Restrict
              ( file : in file_type; s : in Multprec_Span;
                dim,size : in natural32; tol : in double_float;
                sps,sps_last : in out Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the restriction of the embedded polynomial system and
  --   performs the Initialization in the Sampling Machine.
  --   Also the list of samples gets restricted to the span.

    ep : constant Standard_Complex_Poly_Systems.Poly_Sys
       := Sampling_Machine.Embedded_System;
    op : constant Multprec_Complex_Poly_Systems.Poly_Sys
       := Sampling_Machine.Original_System;
    mpep : Multprec_Complex_Poly_Systems.Poly_Sys(ep'range) := Convert(ep);
    mprp : Multprec_Complex_Poly_Systems.Poly_Sys(ep'range)
         := Restrict(s,dim,tol,mpep);
    rp : constant Standard_Complex_Poly_Systems.Poly_Sys(ep'range)
       := Convert(mprp);
    cp : constant Standard_Complex_Poly_Systems.Poly_Sys
       := Collapse_Equations(rp,Dimension(s),dim);
    rop : Multprec_Complex_Poly_Systems.Poly_Sys(op'range)
        := Restrict(s,0,tol,op);
    cop : constant Multprec_Complex_Poly_Systems.Poly_Sys
        := Collapse_Equations(rop,Dimension(s));
    res,res_last : Multprec_Sample_List;
    tmp : Multprec_Sample_List := sps;

  begin
   -- put_line(file,"The restricted embedded polynomial system : "); 
   -- put_line(file,rp);
   -- put_line(file,"The collapsed restricted system : ");
   -- put_line(file,cp);
   -- put_line(file,"The restricted original system : ");
   -- put(file,rop'last,1); new_line(file);
   -- put_line(file,rop);
   -- put_line(file,"The collapsed restricted original system :");
   -- put(file,cop'last,1); new_line(file);
   -- put_line(file,cop);
    Sampling_Machine.Initialize_Restricted(cp,cop,integer32(dim),size);
    while not Is_Null(tmp) loop
      declare
        res_spt : Multprec_Sample := Restrict(s,dim,Head_Of(tmp)); 
       -- res_sol : Multprec_Complex_Solutions.Solution(cp'last)
       --         := Sample_Point(res_spt);
       -- eva1 : Multprec_Complex_Vectors.Vector(rop'range)
       --      := Eval(rop,res_sol.v(cop'range));
       -- eva2 : Multprec_Complex_Vectors.Vector(cop'range)
       --      := Eval(cop,res_sol.v(cop'range));
      begin
        Append(res,res_last,res_spt);
       -- put_line(file,"Evaluation of restrictions :");
       -- put_line(file,eva1);
       -- put_line(file,"Evaluation at collapsed restrictions :");
       -- put_line(file,eva2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
   -- Deep_Clear(sps);
    Multprec_Complex_Poly_Systems.Clear(mpep);
    Multprec_Complex_Poly_Systems.Clear(mprp);
    Multprec_Complex_Poly_Systems.Clear(rop);
    sps := res;
    sps_last := res_last;
  end Multprec_Restrict;

  procedure Multprec_Restrict 
              ( file : in file_type; s : in Multprec_Span;
                dim,size : in natural32; tol : in double_float;
                sps,sps_last : in out Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the restriction of the embedded polynomial system and
  --   performs the Initialization in the Sampling Machine.
  --   Also the list of samples gets restricted to the span.

    ep : constant Standard_Complex_Poly_Systems.Poly_Sys
       := Sampling_Machine.Embedded_System;
    op : constant Multprec_Complex_Poly_Systems.Poly_Sys
       := Sampling_Machine.Original_System;
    mpep : Multprec_Complex_Poly_Systems.Poly_Sys(ep'range) := Convert(ep);
    mprp : Multprec_Complex_Poly_Systems.Poly_Sys(ep'range)
         := Restrict(s,dim,tol,mpep);
    rp : constant Standard_Complex_Poly_Systems.Poly_Sys(ep'range)
       := Convert(mprp);
    cp : constant Standard_Complex_Poly_Systems.Poly_Sys
       := Collapse_Equations(rp,Dimension(s),dim);
    rop : Multprec_Complex_Poly_Systems.Poly_Sys(op'range)
        := Restrict(s,0,tol,op);
    cop : constant Multprec_Complex_Poly_Systems.Poly_Sys
        := Collapse_Equations(rop,Dimension(s));
    res,res_last : Standard_Sample_List;
    tmp : Standard_Sample_List := sps;

  begin
   -- put_line(file,"The restricted embedded polynomial system : "); 
   -- put_line(file,rp);
   -- put_line(file,"The collapsed restricted system : ");
   -- put_line(file,cp);
   -- put_line(file,"The restricted original system : ");
   -- put(file,rop'last,1); new_line(file);
   -- put_line(file,rop);
   -- put_line(file,"The collapsed restricted original system :");
   -- put(file,cop'last,1); new_line(file);
   -- put_line(file,cop);
    Sampling_Machine.Initialize_Restricted(cp,cop,integer32(dim),size);
    while not Is_Null(tmp) loop
      declare
        res_spt : constant Standard_Sample := Restrict(s,dim,Head_Of(tmp)); 
       -- res_sol : Multprec_Complex_Solutions.Solution(cp'last)
       --         := Sample_Point(res_spt);
       -- eva1 : Multprec_Complex_Vectors.Vector(rop'range)
       --      := Eval(rop,res_sol.v(cop'range));
       -- eva2 : Multprec_Complex_Vectors.Vector(cop'range)
       --      := Eval(cop,res_sol.v(cop'range));
      begin
        Append(res,res_last,res_spt);
       -- put_line(file,"Evaluation of restrictions :");
       -- put_line(file,eva1);
       -- put_line(file,"Evaluation at collapsed restrictions :");
       -- put_line(file,eva2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
   -- Deep_Clear(sps);
    Multprec_Complex_Poly_Systems.Clear(mpep);
    Multprec_Complex_Poly_Systems.Clear(mprp);
    Multprec_Complex_Poly_Systems.Clear(rop);
    sps := res;
    sps_last := res_last;
  end Multprec_Restrict;

-- INTERPOLATION USING PLAIN LINEAR PROJECTION OPERATORS :

  procedure Standard_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample; 
                   maxdeg : in natural32; tol : in double_float;
                   p : in Standard_Projector; f : out Standard_Filter ) is

    use Standard_Polynomial_Interpolators;
    dim : constant integer32 := Number_of_Slices(spt);
    sps,sps_last,prev : Standard_Sample_List;
    nb,len : natural32;

  begin
    f := Create(p);
    Append(sps,sps_last,spt);
    for d in 1..maxdeg loop
      nb := Number_of_Terms(d,natural32(dim)+1)-1;
      len := Length_Of(sps);
      Sample(file,full_output,spt,nb-len,sps,sps_last);
      if d = 1
       then Sample_Update(file,f,sps,d);
       else Sample_Update(file,f,prev,d);
      end if;
      Sample(file,full_output,spt,1,sps,sps_last);
      exit when On_Component(file,f,tol,Sample_Point(Head_Of(sps_last)).v);
      prev := sps_last;
    end loop;
    Deep_Clear(sps);
  end Standard_Interpolate;

  procedure Standard_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample;
                   maxdeg : in natural32; tol : in double_float;
                   f : out Standard_Filter ) is

    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    p : constant Standard_Projector := Create(dim+1,n);

  begin
    Standard_Interpolate(file,full_output,spt,maxdeg,tol,p,f);
  end Standard_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in out Standard_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  p : in Multprec_Projector; f : out Multprec_Filter ) is

    mpspt : Multprec_Sample;

  begin
    Refine(file,full_output,spt,mpspt);
    Multprec_Interpolate(file,full_output,mpspt,maxdeg,size,tol,p,f);
  end Multprec_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Multprec_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  p : in Multprec_Projector; f : out Multprec_Filter ) is

    use Multprec_Polynomial_Interpolators;
    dim : constant natural32 := natural32(Number_of_Slices(spt));
    sps,sps_last,prev : Multprec_Sample_List;
    nb,len : natural32;

  begin
    Append(sps,sps_last,spt);
    f := Create(p);
    for d in 1..maxdeg loop
      nb := Number_of_Terms(d,dim+1)-1;
      len := Length_Of(sps);
      Sample(file,full_output,Original(spt),nb-len,sps,sps_last);
      if d = 1
       then Sample_Update(file,f,sps,d);
       else Sample_Update(file,f,prev,d);
      end if;
      Sample(file,full_output,Original(spt),1,sps,sps_last);
      exit when On_Component(file,f,tol,Sample_Point(Head_Of(sps_last)).v);
      prev := sps_last;
    end loop;
    Deep_Clear(sps);
  end Multprec_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in out Standard_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  f : out Multprec_Filter ) is

    mpspt : Multprec_Sample;

  begin
    Refine(file,full_output,spt,mpspt);
    Multprec_Interpolate(file,full_output,mpspt,maxdeg,size,tol,f);
  end Multprec_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Multprec_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  f : out Multprec_Filter ) is

    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    p : constant Multprec_Projector := Create(dim+1,n,size);

  begin
    Multprec_Interpolate(file,full_output,spt,maxdeg,size,tol,p,f);
  end Multprec_Interpolate;

-- INTERPOLATION WITH LINEAR PROJECTIONS AND EXPLOITATION OF SPAN :

  procedure Standard_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Standard_Sample;
                  maxdeg : in natural32; tol : in double_float;
                  p : in Standard_Projector;
                  f : out Standard_Filter; s : out Standard_Span ) is

    use Standard_Polynomial_Interpolators;
    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    sps,sps_last,prev : Standard_Sample_List;
    nb,len : natural32;

  begin
    f := Create(p);
    Append(sps,sps_last,spt);
    s := Null_Standard_Span;
    for d in 1..maxdeg loop
      nb := Number_of_Terms(d,dim+1)-1;
      len := Length_Of(sps);
      Sample(file,full_output,Head_Of(sps),nb-len,sps,sps_last);
      if Empty(s) and (nb <= n+1) then
        s := Create(sps,tol);
        if not Empty(s) then
          Standard_Restrict(file,s,dim,tol,sps,sps_last);
          Shallow_Clear(f);  -- deep clear also erases p
          f := Create(p);
          prev := sps;
        end if;
      end if;
      if d = 1
       then Sample_Update(file,f,sps,d);
       else Sample_Update(file,f,prev,d);
      end if;
      Sample(file,full_output,Head_Of(sps),1,sps,sps_last);
      exit when On_Component(file,f,tol,Sample_Point(Head_Of(sps_last)).v);
      prev := sps_last;
    end loop;
    if not Empty(s)
     then Sampling_Machine.Clear_Restricted;
    -- elsif Degree(f) = 1
    --     then s := Create(sps,tol);
    end if;
    Deep_Clear(sps);
  end Standard_Interpolate;

  procedure Standard_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Standard_Sample; 
                  maxdeg : in natural32; tol : in double_float;
                  f : out Standard_Filter; s : out Standard_Span ) is

    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    p : constant Standard_Projector := Create(dim+1,n);

  begin
    Standard_Interpolate(file,full_output,spt,maxdeg,tol,p,f,s);
  end Standard_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in out Standard_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  p : in Multprec_Projector;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    mpspt : Multprec_Sample;

  begin
    Refine(file,full_output,spt,mpspt);
    Multprec_Interpolate(file,full_output,mpspt,maxdeg,size,tol,p,f,s);
  end Multprec_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in out Standard_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    mpspt : Multprec_Sample;

  begin
    Refine(file,full_output,spt,mpspt);
    Multprec_Interpolate(file,full_output,mpspt,maxdeg,size,tol,f,s);
  end Multprec_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Multprec_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  p : in Multprec_Projector;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    use Multprec_Polynomial_Interpolators;
    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    sps,sps_last,prev : Multprec_Sample_List;
    nb,len : natural32;
    start_spt : Standard_Sample := Original(spt);

  begin
    f := Create(p);
    Append(sps,sps_last,spt);
    s := Null_Multprec_Span;
    for d in 1..maxdeg loop
      nb := Number_of_Terms(d,dim+1)-1;
      len := Length_Of(sps);
      Sample(file,full_output,start_spt,nb-len,sps,sps_last);
      if Empty(s) and (nb <= n+1) then
        s := Create(sps,size,tol);
        if not Empty(s) then
          Multprec_Restrict(file,s,dim,size,tol,sps,sps_last);
          Shallow_Clear(f);  -- deep clear also erases p
          f := Create(p);
          prev := sps;
          start_spt := Restrict(s,dim,Original(spt));
        end if;
      end if;
      if d = 1
       then Sample_Update(file,f,sps,d);
       else Sample_Update(file,f,prev,d);
      end if;
      Sample(file,full_output,start_spt,1,sps,sps_last);
      exit when On_Component(file,f,tol,Sample_Point(Head_Of(sps_last)).v);
      prev := sps_last;
    end loop;
    if not Empty(s)
     then Sampling_Machine.Clear_Restricted;
    -- elsif Degree(f) = 1
    --     then s := Create(sps,size,tol);
    end if;
    Deep_Clear(sps);
  end Multprec_Interpolate;

  procedure Multprec_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Multprec_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    p : constant Multprec_Projector := Create(dim+1,n,size);

  begin
    Multprec_Interpolate(file,full_output,spt,maxdeg,size,tol,p,f,s);
  end Multprec_Interpolate;

-- INTERPOLATION USING CENTRAL PROJECTIONS :

  procedure Standard_Central_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Standard_Sample; 
                  maxdeg : in natural32; tol : in double_float;
                  p : in Standard_Projector;
                  f : out Standard_Filter; s : out Standard_Span ) is

    use Standard_Polynomial_Interpolators;
    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant integer32 := Number_of_Variables(spt);
    sps,sps_last,prev : Standard_Sample_List;
    nb,len : natural32;
    todo : boolean := true;
    testpt : Standard_Sample;

  begin
    f := Create(p);
    Append(sps,sps_last,spt);
    s := Null_Standard_Span;
    for d in 1..maxdeg loop
      nb := Number_of_Terms(d,dim+1)-1;
      len := Length_Of(sps);
      Sample(file,full_output,Head_Of(sps),nb-len,sps,sps_last);
      if Empty(s) and (integer32(nb) <= n+1) then
        s := Create(sps,tol);
        if not Empty(s) then
          Standard_Restrict(file,s,dim,tol,sps,sps_last);
          Shallow_Clear(f);  -- deep clear also erases p
          f := Create(p);
          prev := sps;
          for i in 1..Dimension(s)-dim-1 loop
            declare
              basept : Standard_Sample;
              basehyp : constant Standard_Complex_Vectors.Vector(0..n)
                      := Random_Vector(0,n);
            begin
              put(file,"Sampling central point "); put(file,i,1);
              put_line(file,".");
              Sample(file,full_output,Head_Of(sps),basept);
              Central_Update(f,basept,basehyp);
              Sample_Update(file,f,sps,d);
              if i = 1 then
                Sample(file,full_output,Head_Of(sps),testpt);
                todo := false;
              end if;
              if On_Component(file,f,tol,
                              Sample_Point(Head_Of(sps_last)).v) then
                Sampling_Machine.Clear_Restricted;
                Deep_Clear(sps); return;
              elsif i < Dimension(s)-dim-1 then
                Shallow_Clear(f);
                f := Create(p);
              else
                Append(sps,sps_last,testpt);
              end if;
            end;
          end loop;
        end if;
      end if;
      if todo then
        if d = 1
         then Sample_Update(file,f,sps,d);
         else Sample_Update(file,f,prev,d);
        end if;
        Sample(file,full_output,Head_Of(sps),1,sps,sps_last);
        exit when On_Component(file,f,tol,Sample_Point(Head_Of(sps_last)).v);
      else
        todo := true;
      end if;
      prev := sps_last;
    end loop;
    if not Empty(s)
     then Sampling_Machine.Clear_Restricted;
    -- elsif Degree(f) = 1
    --     then s := Create(sps,tol);
    end if;
    Deep_Clear(sps);
  end Standard_Central_Interpolate;

  procedure Standard_Central_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Standard_Sample;
                  maxdeg : in natural32; tol : in double_float;
                  f : out Standard_Filter; s : out Standard_Span ) is

    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    p : constant Standard_Projector := Create(dim+1,n);

  begin
    Standard_Central_Interpolate(file,full_output,spt,maxdeg,tol,p,f,s);
  end Standard_Central_Interpolate;

  procedure Multprec_Central_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in out Standard_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  p : in Multprec_Projector;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    mpspt : Multprec_Sample;

  begin
    Refine(file,full_output,spt,mpspt);
    Multprec_Central_Interpolate(file,full_output,mpspt,maxdeg,size,tol,p,f,s);
  end Multprec_Central_Interpolate;

  procedure Multprec_Central_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in out Standard_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    mpspt : Multprec_Sample;

  begin
    Refine(file,full_output,spt,mpspt);
    Multprec_Central_Interpolate(file,full_output,mpspt,maxdeg,size,tol,f,s);
  end Multprec_Central_Interpolate;

  procedure Multprec_Central_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Multprec_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  p : in Multprec_Projector;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    use Multprec_Polynomial_Interpolators;
    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant integer32 := Number_of_Variables(spt);
    sps,sps_last,prev : Multprec_Sample_List;
    nb,len : natural32;
    start_spt : Standard_Sample := Original(spt);
    todo : boolean := true;
    testpt : Multprec_Sample;

  begin
    f := Create(p);
    s := Null_Multprec_Span;
    for d in 1..maxdeg loop
      nb := Number_of_Terms(d,dim+1)-1;
      len := Length_Of(sps);
      Sample(file,full_output,start_spt,nb-len,sps,sps_last);
     -- Parallel_Sample(file,full_output,start_spt,nb-len,sps,sps_last);
      if Empty(s) and (integer32(nb) <= n+1) then
        s := Create(sps,size,tol);
        if not Empty(s) then
          Multprec_Restrict(file,s,dim,size,tol,sps,sps_last);
          Shallow_Clear(f);  -- deep clear also erases p
          f := Create(p);
          prev := sps;
          start_spt := Restrict(s,dim,Original(spt));
          for i in 1..Dimension(s)-dim-1 loop
            declare
              basept : Multprec_Sample;
              basehyp : constant Multprec_Complex_Vectors.Vector(0..n)
                      := Create(Random_Vector(0,n));
                      -- := Random_Vector(0,n,size);
            begin
              put(file,"Sampling central point "); put(file,i,1);
              put_line(file,".");
              Sample(file,full_output,start_spt,basept);
             -- Parallel_Sample(file,full_output,start_spt,basept);
              Central_Update(f,basept,basehyp);
              Sample_Update(file,f,sps,d);
              if i = 1
               then Sample(file,full_output,start_spt,testpt);
                   -- Parallel_Sample(file,full_output,start_spt,testpt);
                    todo := false;
              end if;
              if On_Component(file,f,tol,Sample_Point(testpt).v) then
                Sampling_Machine.Clear_Restricted;
                Deep_Clear(sps); return;
              elsif i < Dimension(s)-dim-1 then
                Shallow_Clear(f);
                f := Create(p);
              else
                Append(sps,sps_last,testpt);
              end if;
            end;
          end loop;
        end if;
      end if;
      if todo
       then if d = 1
             then Sample_Update(file,f,sps,d);
             else Sample_Update(file,f,prev,d);
            end if;
            Sample(file,full_output,start_spt,1,sps,sps_last);
           -- Parallel_Sample(file,full_output,start_spt,1,sps,sps_last);
            exit when On_Component
                           (file,f,tol,Sample_Point(Head_Of(sps_last)).v);
       else todo := true;
      end if;
      prev := sps_last;
    end loop;
    if not Empty(s)
     then Sampling_Machine.Clear_Restricted;
    -- elsif Degree(f) = 1
    --     then s := Create(sps,size,tol);
    end if;
    Deep_Clear(sps);
  end Multprec_Central_Interpolate;

  procedure Multprec_Central_Interpolate
                ( file : in file_type; full_output : in boolean;
                  spt : in Multprec_Sample;
                  maxdeg,size : in natural32; tol : in double_float;
                  f : out Multprec_Filter; s : out Multprec_Span ) is

    dim : constant natural32 := natural32(Number_of_Slices(spt));
    n : constant natural32 := natural32(Number_of_Variables(spt));
    p : constant Multprec_Projector := Create(dim+1,n,size);

  begin
    Multprec_Central_Interpolate(file,full_output,spt,maxdeg,size,tol,p,f,s);
  end Multprec_Central_Interpolate;

-- NEWTON INTERPOLATION USING DIVIDED DIFFERENCES :

  procedure Standard_Create_Span
                ( file : in file_type; sps : in Standard_Sample_List;
                  dim : out natural32; sp : out Standard_Span;
                  spansps,rst_spansps : out Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the linear span of the space spanned by the samples,
  --   starting the sampling at the first element of the list.

  -- ON ENTRY :
  --   file       file to write diagnostics on;
  --   sps        list of samples.

  -- ON RETURN :
  --   dim        number of slices;
  --   sp         linear span of the component;
  --   spansps    extra samples computed to determine the span;
  --   rst_spansps is the list of samples restricted to the span.

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;
    spansps_last,rst_spansps_last : Standard_Sample_List;
    spt : constant Standard_Sample := Head_Of(sps);

  begin
    tstart(timer);
    dim := natural32(Number_of_Slices(spt));
    Create_Span(spt,tol,spansps,spansps_last,sp);
    rst_spansps := sps;
    if Empty(sp) then
      put_line(file,"The span of the component is empty.");
    else
      put_line(file,"The span of the component : "); put(file,sp);
      Standard_Restrict(file,sp,dim,tol,rst_spansps,rst_spansps_last);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creating Span and Subspace Restriction");
    new_line(file);
  end Standard_Create_Span;

  procedure Multprec_Create_Span
              ( file : in file_type; sps : in Standard_Sample_List;
                size : in natural32; dim : out natural32;
                sp : out Multprec_Span; spansps : out Multprec_Sample_List;
                rst_spansps : out Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the linear span of the space spanned by the samples,
  --   starting the sampling at the first element of the list.

  -- ON ENTRY :
  --   file       file to write diagnostics on;
  --   sps        list of samples;
  --   size       size of the multi-precision numbers.

  -- ON RETURN :
  --   dim        number of slices;
  --   sp         linear span of the component;
  --   spansps    extra samples computed to determine the span;
  --   rst_spansps is the list of samples restricted to the span.

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;
    spansps_last : Multprec_Sample_List;
    rst_spansps_last : Standard_Sample_List;
    spt : constant Standard_Sample := Head_Of(sps);

  begin
    tstart(timer);
    dim := natural32(Number_of_Slices(spt));
    Create_Span(spt,tol,size,spansps,spansps_last,sp);
    rst_spansps := sps;
    if Empty(sp) then
      put_line(file,"The span of the component is empty.");
    else
      put_line(file,"The span of the component : "); put(file,sp);
      Multprec_Restrict(file,sp,dim,size,tol,rst_spansps,rst_spansps_last);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creating Span and Subspace Restriction");
    new_line(file);
  end Multprec_Create_Span;

  procedure Standard_Create_Grid
              ( file : in file_type; sps : in Standard_Sample_List;
                eps,dist : out double_float;
                grid
                  : out Standard_Stacked_Sample_Grids.Stacked_Sample_Grid ) is

  -- DESCRIPTION :
  --   Creates a rectangular grid to interpolate a curve.

  -- ON ENTRY :
  --   file       to write intermediate diagnostics on;
  --   sps        list of samples.

  -- ON RETURN :
  --   eps        accuracy of the samples in the grid;
  --   dist       minimal distance between the points in the grid;
  --   grid       the grid of samples points.

    timer : Timing_Widget;
    maxerr : double_float;

    use Standard_Stacked_Sample_Grids;

  begin
    tstart(timer);
    grid := Create(file,sps);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creating the grid of samples");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,maxerr,3); new_line(file);
    eps := maxerr;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid.");
    new_line(file);
    dist := 1.0;
  end Standard_Create_Grid;

  procedure Standard_Create_Newton_Interpolator1
                ( file : in file_type;
                  grid : in Array_of_Standard_Sample_Lists;
                  q : out Standard_Divided_Differences.Newton_Interpolator1;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Newton interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          Newton form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;
    maxerr : double_float;

    use Standard_Divided_Differences;

  begin
    tstart(timer);
    q := Create(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Newton Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := maxerr;
    new_line(file);
    print_times(file,timer,"Testing quality of Newton Interpolator");
    new_line(file);
  end Standard_Create_Newton_Interpolator1;

  procedure Standard_Create_Newton_Interpolator
                ( file : in file_type; len : in natural32;
                  grid : in Standard_Stacked_Sample_Grids.Stacked_Sample_Grid;
                  q : out Standard_Divided_Differences.Newton_Taylor_Form;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Newton interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          Newton form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    use Standard_Stacked_Sample_Grids;
    use Standard_Divided_Differences;

    timer : Timing_Widget;
    maxerr : double_float;
    y : constant Standard_Complex_Vectors.Vector(0..integer32(len))
      := Random_Vector(0,integer32(len));
    c : constant Standard_Complex_Numbers.Complex_Number := y(0);

  begin
    tstart(timer);
    q := Create(grid,c,y);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Newton Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := maxerr;
    new_line(file);
    print_times(file,timer,"Testing quality of Newton Interpolator");
    new_line(file);
  end Standard_Create_Newton_Interpolator;

  procedure Multprec_Create_Newton_Interpolator1
                ( file : in file_type;
                  grid : in Array_of_Multprec_Sample_Lists;
                  q : out Multprec_Divided_Differences.Newton_Interpolator1;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Newton interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          Newton form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;
    maxerr : Floating_Number;

    use Multprec_Divided_Differences;

  begin
    tstart(timer);
    q := Create(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Newton Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := Round(maxerr);
    Clear(maxerr);
    new_line(file);
    print_times(file,timer,"Testing the quality of the Newton Interpolator");
    new_line(file);
  end Multprec_Create_Newton_Interpolator1;

  procedure Multprec_Create_Newton_Interpolator
                ( file : in file_type; len,size : in natural32;
                  grid : in Multprec_Stacked_Sample_Grids.Stacked_Sample_Grid;
                  q : out Multprec_Divided_Differences.Newton_Taylor_Form;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Newton interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          Newton form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    use Multprec_Stacked_Sample_Grids;
    use Multprec_Divided_Differences;

    timer : Timing_Widget;
    maxerr : Floating_Number;
    y : constant Multprec_Complex_Vectors.Vector(0..integer32(len))
      := Random_Vector(0,integer32(len),size);
    c : constant Multprec_Complex_Numbers.Complex_Number := y(0);

  begin
    tstart(timer);
    q := Create(grid,c,y);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Newton Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := Round(maxerr);
    new_line(file);
    print_times(file,timer,"Testing quality of Newton Interpolator");
    new_line(file);
  end Multprec_Create_Newton_Interpolator;

  procedure Standard_Test_Evaluation
                ( file : in file_type; sps : Standard_Sample_List;
                  dim : in natural32; sp : in Standard_Span;
                  q : in Standard_Divided_Differences.Newton_Interpolator1;
                  testres : out double_float ) is

  -- DESCRIPTION :
  --   Evaluates the interpolator at the sample points, eventually,
  --   in case the span is not empty, after restriction to this span.

  -- ON ENTRY :
  --   file       to write intermedate diagnostics on;
  --   dim        number of slices;
  --   sps        test samples;
  --   sp         linear span of the component;
  --   q          Newton form of interpolating polynomial;

  -- ON RETURN :
  --   testres    maximal residual in the test samples.

    use Standard_Divided_Differences;
    eq : Newton_Form_Evaluator1;
    timer : Timing_Widget;
    tmp : Standard_Sample_List := sps;
    maxerr : double_float := 0.0;
    abseva : double_float;

  begin
    tstart(timer);
    eq := Create(q);
    for i in 1..Length_Of(sps) loop
      put(file,"Residual of evaluation at test sample ");
      put(file,i,2); put(file," : ");
      declare
        test_spt : constant Standard_Sample := Head_Of(tmp);
        rest_spt : Standard_Sample;
      begin
        if Empty(sp) then
          abseva := AbsVal(Eval(eq,Sample_Point(test_spt).v));
          put(file,abseva,3);
        else
          rest_spt := Restrict(sp,dim,test_spt);
          abseva := AbsVal(Eval(eq,Sample_Point(rest_spt).v));
          put(file,abseva,3);
          Deep_Clear(rest_spt);
          if abseva > maxerr
           then maxerr := abseva;
          end if;
        end if;
      end;
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(eq);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Evaluation at test samples");
    new_line(file);
    testres := maxerr;
  end Standard_Test_Evaluation;

  procedure Multprec_Test_Evaluation
                ( file : in file_type; sps : Multprec_Sample_List;
                  dim : in natural32; sp : in Multprec_Span;
                  q : in Multprec_Divided_Differences.Newton_Interpolator1;
                  testres : out double_float ) is

  -- DESCRIPTION :
  --   Evaluates the interpolator at the sample points, eventually,
  --   in case the span is not empty, after restriction to this span.

  -- ON ENTRY :
  --   file       to write intermedate diagnostics on;
  --   dim        number of slices;
  --   sps        test samples;
  --   sp         linear span of the component;
  --   q          Newton form of interpolating polynomial;

  -- ON RETURN :
  --   testres    maximal residual in the test samples.

    use Multprec_Divided_Differences;
    timer : Timing_Widget;
    eq : Newton_Form_Evaluator1;
    maxerr : Floating_Number := Create(integer(0));
    abseva : Floating_Number;
    eva : Multprec_Complex_Numbers.Complex_Number;
    tmp : Multprec_Sample_List := sps;

  begin
    tstart(timer);
    eq := Create(q);
    for i in 1..Length_Of(sps) loop
      declare
        test_spt : constant Multprec_Sample := Head_Of(tmp);
        rest_spt : Multprec_Sample;
      begin
        if Empty(sp)
         then eva := Eval(eq,Sample_Point(test_spt).v);
         else rest_spt := Restrict(sp,dim,test_spt);
              eva := Eval(eq,Sample_Point(rest_spt).v);
             -- Deep_Clear(rest_spt);
        end if;
        abseva := AbsVal(eva);
        put(file,"Residual of evaluation at test sample ");
        put(file,i,2); put(file," : ");
        put(file,abseva,3); new_line(file);
        if abseva > maxerr
         then Copy(abseva,maxerr);
        end if;
        Clear(eva); Clear(abseva);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Evaluation at test samples");
    new_line(file);
    testres := Round(maxerr);
  end Multprec_Test_Evaluation;

  procedure Standard_Newton_Interpolate1
                ( file : in file_type; sps : in Standard_Sample_List;
                  q : out Standard_Divided_Differences.Newton_Interpolator1;
                  eps,dist,gridres,testres : out double_float ) is

    use Standard_Divided_Differences;

    dim : natural32;
    spsps,sps1 : Standard_Sample_List;
    sp : Standard_Span;
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Standard_Sample_Lists(0..integer32(len));

  begin
    Standard_Create_Span(file,sps,dim,sp,spsps,sps1);
    Standard_Rectangular_Grid_Creator(file,sps1,len,grid,eps,dist);
    Standard_Create_Newton_Interpolator1(file,grid,q,gridres);
    Standard_Test_Evaluation(file,spsps,dim,sp,q,testres);
    Deep_Clear(grid); Deep_Clear(spsps);
    if not Empty(sp)
     then Sampling_Machine.Clear_Restricted;
    end if;
  end Standard_Newton_Interpolate1;

  procedure Standard_Newton_Interpolate
                ( file : in file_type; sps : in Standard_Sample_List;
                  q : out Standard_Divided_Differences.Newton_Taylor_Form;
                  eps,dist,gridres,testres : out double_float ) is

    use Standard_Stacked_Sample_Grids;
    use Standard_Divided_Differences;

    timer : Timing_Widget;
   -- dim : natural32;
    tmp,testsps,testsps_last : Standard_Sample_List;
    len : constant natural32 := Length_Of(sps);
    dimsps : constant integer32 := Number_of_Slices(Head_Of(sps));
    grid : Stacked_Sample_Grid(dimsps,integer32(len));
    eva : Standard_Complex_Numbers.Complex_Number;
    abseva : double_float;
    maxerr : double_float := -1.0;

  begin
   -- Standard_Create_Span(file,sps,dim,sp,spsps,sps1);
   -- Standard_Create_Grid(file,sps1,eps,dist,grid);
    tstart(timer);
    grid := Create(file,sps);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid");
    tstart(timer);
    eps := Maximal_Error(grid);
    dist := Minimal_Distance(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid");
    Standard_Create_Newton_Interpolator(file,len,grid,q,gridres);
    tstart(timer);
    Sample(Head_Of(sps),5,testsps,testsps_last);
    tmp := testsps;
    new_line(file);
    for i in 1..integer32(5) loop
      eva := Eval(q,Sample_Point(Head_Of(tmp)).v);
      abseva := AbsVal(eva);
      put(file,"Residual of evaluation at test sample ");
      put(file,i,2); put(file," : ");
      put(file,abseva,3); new_line(file);
      if abseva > maxerr
       then maxerr := abseva;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    testres := maxerr;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Generating and evaluating at test samples");
  end Standard_Newton_Interpolate;

  procedure Multprec_Newton_Interpolate1
                ( file : in file_type; sps : in Standard_Sample_List;
                  size : in natural32;
                  q : out Multprec_Divided_Differences.Newton_Interpolator1;
                  eps,dist,gridres,testres : out double_float ) is

    use Multprec_Divided_Differences;

    dim : natural32;
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Multprec_Sample_Lists(0..integer32(len));
    spsps : Multprec_Sample_List;
    sps1 : Standard_Sample_List;
    sp : Multprec_Span;

  begin
    Multprec_Create_Span(file,sps,size,dim,sp,spsps,sps1);
    Multprec_Rectangular_Grid_Creator(file,sps1,len,size,grid,eps,dist);
    Multprec_Create_Newton_Interpolator1(file,grid,q,gridres);
    Multprec_Test_Evaluation(file,spsps,dim,sp,q,testres);
    Deep_Clear(grid); Deep_Clear(spsps);
    if not Empty(sp)
     then Sampling_Machine.Clear_Restricted;
    end if;
  end Multprec_Newton_Interpolate1;

  procedure Multprec_Newton_Interpolate
                ( file : in file_type; sps : in Standard_Sample_List;
                  size : in natural32;
                  q : out Multprec_Divided_Differences.Newton_Taylor_Form;
                  eps,dist,gridres,testres : out double_float ) is

    use Multprec_Stacked_Sample_Grids;
    use Multprec_Divided_Differences;

    timer : Timing_Widget;
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    len : constant natural32 := Length_Of(sps);
    grid : Stacked_Sample_Grid(dim,integer32(len));
    maxerr,mindist : Floating_Number;
    testsps,testsps_last,tmp : Multprec_Sample_List;
    eva : Multprec_Complex_Numbers.Complex_Number;
    abseva : Floating_Number;
    fltabseva : double_float;

  begin
    tstart(timer);
    grid := Create(file,sps,size);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid");
    tstart(timer);
    maxerr := Maximal_Error(grid);
    eps := Round(maxerr); Clear(maxerr);
    mindist := Minimal_Distance(grid);
    dist := Round(mindist); Clear(mindist);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid");
    Multprec_Create_Newton_Interpolator(file,len,size,grid,q,gridres);
    Sample(Head_Of(sps),5,testsps,testsps_last);
    testres := -1.0;
    tmp := testsps;
    new_line(file);
    for i in 1..integer32(5) loop
      eva := Eval(q,Sample_Point(Head_Of(tmp)).v);
      abseva := AbsVal(eva);
      put(file,"Residual of evaluation at test sample ");
      put(file,i,2); put(file," : ");
      put(file,abseva,3); new_line(file);
      fltabseva := Round(maxerr);
      Clear(maxerr);
      if fltabseva > testres
       then testres := fltabseva;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Generating and Evaluating at Test Samples");
  end Multprec_Newton_Interpolate;

-- INTERPOLATION WITH TRACES :

  procedure Standard_Create_Trace_Interpolator1
                ( file : in file_type;
                  grid : in Array_of_Standard_Sample_Lists;
                  q : out Standard_Trace_Interpolators.Trace_Interpolator1;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Traces interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          trace form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;
    maxerr : double_float;

    use Standard_Trace_Interpolators;

  begin
    tstart(timer);
    q := Create(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Traces Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := maxerr;
    new_line(file);
    print_times(file,timer,"Testing Quality of Traces Interpolator");
    new_line(file);
  end Standard_Create_Trace_Interpolator1;

  procedure Standard_Create_Power_Trace_Interpolator1
                ( file : in file_type;
                  grid : in Array_of_Standard_Sample_Lists;
                  q : out Standard_Trace_Interpolators.Trace_Interpolator1;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Traces interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          trace form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;
    maxerr : double_float;

    use Standard_Trace_Interpolators;

  begin
    tstart(timer);
    q := Create_on_Triangle(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Traces Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := maxerr;
    new_line(file);
    print_times(file,timer,"Testing Quality of Traces Interpolator");
    new_line(file);
  end Standard_Create_Power_Trace_Interpolator1;

  procedure Standard_Create_Trace_Interpolator
                ( file : in file_type; d : in natural32;
                  grid : in Standard_Stacked_Sample_Grids.Stacked_Sample_Grid;
                  q : out Standard_Trace_Interpolators.Trace_Interpolator;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Traces interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          trace form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;
    maxerr : double_float;

    use Standard_Trace_Interpolators;

  begin
    tstart(timer);
    q := Create(grid,integer32(d));
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Traces Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := maxerr;
    new_line(file);
    print_times(file,timer,"Testing Quality of Traces Interpolator");
    new_line(file);
  end Standard_Create_Trace_Interpolator;

  procedure Standard_Create_Linear_Trace1
                ( file : in file_type;
                  grid : in Array_of_Standard_Sample_Lists;
                  q : out Standard_Complex_Vectors.Vector;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the linear trace through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          linear trace part of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;

    use Standard_Trace_Interpolators;

  begin
    tstart(timer);
    q := Create(grid,1);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Linear Trace");
    new_line(file);
    tstart(timer);
    gridres := Maximal_Error(q,1,grid(0..1));
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,gridres,3); new_line(file);
    new_line(file);
    print_times(file,timer,"Testing Quality of the Linear Trace");
    new_line(file);
  end Standard_Create_Linear_Trace1;

  procedure Multprec_Create_Trace_Interpolator1
                ( file : in file_type;
                  grid : in Array_of_Multprec_Sample_Lists;
                  q : out Multprec_Trace_Interpolators.Trace_Interpolator1;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Traces interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          trace form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;
    maxerr : Floating_Number;

    use Multprec_Trace_Interpolators;

  begin
    tstart(timer);
    q := Create(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Traces Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := Round(maxerr);
    new_line(file);
    print_times(file,timer,"Testing quality of Traces Interpolator");
    new_line(file);
  end Multprec_Create_Trace_Interpolator1;

  procedure Multprec_Create_Power_Trace_Interpolator1
                ( file : in file_type;
                  grid : in Array_of_Multprec_Sample_Lists;
                  size : in natural32;
                  q : out Multprec_Trace_Interpolators.Trace_Interpolator1;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Traces interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          trace form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    timer : Timing_Widget;
    maxerr : Floating_Number;

    use Multprec_Trace_Interpolators;

  begin
    tstart(timer);
    q := Create_on_Triangle(file,grid,size);     -- "file" for diagnostics
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Traces Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := Round(maxerr);
    new_line(file);
    print_times(file,timer,"Testing quality of Traces Interpolator");
    new_line(file);
  end Multprec_Create_Power_Trace_Interpolator1;

  procedure Multprec_Create_Trace_Interpolator
                ( file : in file_type; d : in natural32;
                  grid : in Multprec_Stacked_Sample_Grids.Stacked_Sample_Grid;
                  q : out Multprec_Trace_Interpolators.Trace_Interpolator;
                  gridres : out double_float ) is

  -- DESCRIPTION :
  --   Computes the Traces interpolator through the samples in the grid.

  -- ON ENTRY :
  --   file       for intermediate diagnostics;
  --   grid       grid of sample points.

  -- ON RETURN :
  --   q          trace form of interpolating polynomial;
  --   gridres    maximal residual of grid samples in q.

    use Multprec_Stacked_Sample_Grids;
    use Multprec_Trace_Interpolators;

    timer : Timing_Widget;
    maxerr : Floating_Number;

  begin
    tstart(timer);
    q := Create(grid,integer32(d));
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Constructing the Traces Interpolator");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(q,grid);
    tstop(timer);
    put(file,"Maximal residual of evaluation at grid   : ");
    put(file,maxerr,3); new_line(file);
    gridres := Round(maxerr);
    new_line(file);
    print_times(file,timer,"Testing quality of Traces Interpolator");
  end Multprec_Create_Trace_Interpolator;

  procedure Standard_Test_Evaluation
                ( file : in file_type; sps : Standard_Sample_List;
                  dim : in natural32; sp : in Standard_Span;
                  q : in Standard_Trace_Interpolators.Trace_Interpolator1;
                  testres : out double_float ) is

  -- DESCRIPTION :
  --   Evaluates the interpolator at the sample points, eventually,
  --   in case the span is not empty, after restriction to this span.

  -- ON ENTRY :
  --   file       to write intermedate diagnostics on;
  --   dim        number of slices;
  --   sps        test samples;
  --   sp         linear span of the component;
  --   q          trace form of interpolating polynomial;

  -- ON RETURN :
  --   testres    maximal residual in the test samples.

    use Standard_Trace_Interpolators;

    timer : Timing_Widget;
    tmp : Standard_Sample_List := sps;
    maxerr : double_float := 0.0;
    abseva : double_float;

  begin
    tstart(timer);
    for i in 1..Length_Of(sps) loop
      put(file,"Residual of evaluation at test sample ");
      put(file,i,2); put(file," : ");
      declare
        test_spt : Standard_Sample := Head_Of(tmp);
        rest_spt : Standard_Sample;
      begin
        if Empty(sp) then
          abseva := AbsVal(Eval(q,Sample_Point(test_spt).v));
          put(file,abseva,3);
        else
          rest_spt := Restrict(sp,dim,test_spt);
          abseva := AbsVal(Eval(q,Sample_Point(rest_spt).v));
          put(file,abseva,3);
          Deep_Clear(rest_spt);
          if abseva > maxerr
           then maxerr := abseva;
          end if;
        end if;
      end;
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Evaluation of Trace Form at Test Samples");
    new_line(file);
    testres := maxerr;
  end Standard_Test_Evaluation;

  procedure Standard_Test_Linear_Evaluation
                ( file : in file_type; sps : Standard_Sample_List;
                  dim : in natural32; sp : in Standard_Span;
                  q : in Standard_Complex_Vectors.Vector;
                  testres : out double_float ) is

  -- DESCRIPTION :
  --   Evaluates the linear trace at the sample points, eventually,
  --   in case the span is not empty, after restriction to this span.

  -- ON ENTRY :
  --   file       to write intermedate diagnostics on;
  --   dim        number of slices;
  --   sps        test samples;
  --   sp         linear span of the component;
  --   q          trace form of interpolating polynomial;

  -- ON RETURN :
  --   testres    maximal residual in the test samples.

    use Standard_Trace_Interpolators;

    timer : Timing_Widget;
   -- restsps : Standard_Sample_List;
    val,eva : Standard_Complex_Numbers.Complex_Number;

  begin
    tstart(timer);
    --if Empty(sp)
    -- then
    Eval_Trace(q,1,sps,val,eva);
    -- else restsps := Restrict(sp,dim,sps);
    --      Eval_Trace(q,1,restsps,val,eva);
    --end if;
    tstop(timer);
    put(file,"Value at trace : "); put(file,val); new_line(file);
    put(file,"Computed sum   : "); put(file,eva); new_line(file);
    testres := AbsVal(val-eva);
    put(file,"Absolute Value of Difference : "); put(file,testres);
    new_line(file);
    new_line(file);
    print_times
      (file,timer,"Evaluation of Linear Trace at Test Samples");
    new_line(file);
  end Standard_Test_Linear_Evaluation;

  procedure Multprec_Test_Evaluation
                ( file : in file_type; sps : Multprec_Sample_List;
                  dim : in natural32; sp : in Multprec_Span;
                  q : in Multprec_Trace_Interpolators.Trace_Interpolator1;
                  testres : out double_float ) is

  -- DESCRIPTION :
  --   Evaluates the interpolator at the sample points, eventually,
  --   in case the span is not empty, after restriction to this span.

  -- ON ENTRY :
  --   file       to write intermedate diagnostics on;
  --   dim        number of slices;
  --   sps        test samples;
  --   sp         linear span of the component;
  --   q          trace form of interpolating polynomial;

  -- ON RETURN :
  --   testres    maximal residual in the test samples.

    use Multprec_Trace_Interpolators;

    timer : Timing_Widget;
    tmp : Multprec_Sample_List := sps;
    maxerr : Floating_Number := Create(integer(0));
    abseva : Floating_Number;

  begin
    tstart(timer);
    for i in 1..Length_Of(sps) loop
      put(file,"Residual of evaluation at test sample ");
      put(file,i,2); put(file," : ");
      declare
        test_spt : constant Multprec_Sample := Head_Of(tmp);
        rest_spt : Multprec_Sample;
      begin
        if Empty(sp) then
          put(file,AbsVal(Eval(q,Sample_Point(test_spt).v)),3);
        else
          rest_spt := Restrict(sp,dim,test_spt);
          abseva := AbsVal(Eval(q,Sample_Point(rest_spt).v));
          put(file,abseva,3);
          Deep_Clear(rest_spt);
          if abseva > maxerr
           then Copy(abseva,maxerr);
          end if;
        end if;
      end;
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Evaluation at test samples");
    new_line(file);
    testres := Round(maxerr);
  end Multprec_Test_Evaluation;

  procedure Standard_Trace_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Trace_Interpolators.Trace_Interpolator1;
                   eps,dist,gridres,testres : out double_float ) is

    dim : natural32;
    spsps,sps1 : Standard_Sample_List;
    sp : Standard_Span;
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Standard_Sample_Lists(0..integer32(len));

  begin
    Standard_Create_Span(file,sps,dim,sp,spsps,sps1);
    Standard_Rectangular_Grid_Creator(file,sps1,len,grid,eps,dist);
    Standard_Create_Trace_Interpolator1(file,grid,t,gridres);
    Standard_Test_Evaluation(file,spsps,dim,sp,t,testres);
    Deep_Clear(grid); Deep_Clear(spsps);
    if not Empty(sp)
     then Sampling_Machine.Clear_Restricted;
    end if;
  end Standard_Trace_Interpolate1;

  procedure Standard_Power_Trace_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Trace_Interpolators.Trace_Interpolator1;
                   eps,dist,gridres,testres : out double_float ) is

    dim : natural32;
    spsps,sps1 : Standard_Sample_List;
    sp : Standard_Span;
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Standard_Sample_Lists(0..integer32(len));

  begin
    Standard_Create_Span(file,sps,dim,sp,spsps,sps1);
    Standard_Triangular_Grid_Creator(file,sps1,len,grid,eps,dist);
    Standard_Create_Power_Trace_Interpolator1(file,grid,t,gridres);
    Standard_Test_Evaluation(file,spsps,dim,sp,t,testres);
    Deep_Clear(grid); Deep_Clear(spsps);
    if not Empty(sp)
     then Sampling_Machine.Clear_Restricted;
    end if;
  end Standard_Power_Trace_Interpolate1;

  procedure Standard_Trace_Interpolate
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Trace_Interpolators.Trace_Interpolator;
                   eps,dist,gridres,testres : out double_float ) is

    use Standard_Stacked_Sample_Grids;
    use Standard_Trace_Interpolators;

    timer : Timing_Widget;
    len : constant natural32 := Length_Of(sps);
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    grid : Stacked_Sample_Grid(dim,integer32(len));
    testsps,testsps_last,tmp : Standard_Sample_List;
    eva : Standard_Complex_Numbers.Complex_Number;
    abseva,maxerr : double_float;

  begin
    tstart(timer);
    grid := Create_Full(file,sps);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid");
    tstart(timer);
    eps := Maximal_Error(grid);
    dist := Minimal_Distance(grid);
    tstop(timer);
    new_line(file);
    put(file,"The maximal error on samples in grid : ");
    put(file,eps,3); new_line(file);
    put(file,"The minimal distance between samples : ");
    put(file,dist,3); new_line(file);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid");
    Standard_Create_Trace_Interpolator(file,len,grid,t,gridres);
    tstart(timer);
    Sample(Head_Of(sps),5,testsps,testsps_last);
    tmp := testsps;
    maxerr := -1.0;
    new_line(file);
    for i in 1..integer32(5) loop
      eva := Eval(t,Sample_Point(Head_Of(tmp)).v);
      abseva := AbsVal(eva);
      put(file,"Residual of evaluation at test sample ");
      put(file,i,2); put(file," : ");
      put(file,abseva,3); new_line(file);
      if abseva > maxerr
       then maxerr := abseva;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    testres := maxerr;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Generating and evaluating at test samples");
  end Standard_Trace_Interpolate;

  procedure Multprec_Trace_Interpolate1
                ( file : in file_type; sps : in Standard_Sample_List;
                  size : in natural32;
                  t : out Multprec_Trace_Interpolators.Trace_Interpolator1;
                  eps,dist,gridres,testres : out double_float ) is

    dim : natural32;
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Multprec_Sample_Lists(0..integer32(len));
    spsps : Multprec_Sample_List;
    sps1 : Standard_Sample_List;
    sp : Multprec_Span;

  begin
    Multprec_Create_Span(file,sps,size,dim,sp,spsps,sps1);
    Multprec_Rectangular_Grid_Creator(file,sps1,len,size,grid,eps,dist);
    Multprec_Create_Trace_Interpolator1(file,grid,t,gridres);
    Multprec_Test_Evaluation(file,spsps,dim,sp,t,testres);
    Deep_Clear(grid); --Deep_Clear(spsps); --> CAUSES SEGMENTATION FAULT
    if not Empty(sp)
     then Sampling_Machine.Clear_Restricted;
    end if;
  end Multprec_Trace_Interpolate1;

  procedure Multprec_Power_Trace_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   size : in natural32;
                   t : out Multprec_Trace_Interpolators.Trace_Interpolator1;
                   eps,dist,gridres,testres : out double_float ) is

    dim : natural32;
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Multprec_Sample_Lists(0..integer32(len));
    spsps : Multprec_Sample_List;
    sps1 : Standard_Sample_List;
    sp : Multprec_Span;

  begin
    Multprec_Create_Span(file,sps,size,dim,sp,spsps,sps1);
    Multprec_Triangular_Grid_Creator(file,sps1,len,size,grid,eps,dist);
    Multprec_Create_Power_Trace_Interpolator1(file,grid,size,t,gridres);
    Multprec_Test_Evaluation(file,spsps,dim,sp,t,testres);
    Deep_Clear(grid); --Deep_Clear(spsps); --> CAUSES SEGMENTATION FAULT
    if not Empty(sp)
     then Sampling_Machine.Clear_Restricted;
    end if;
  end Multprec_Power_Trace_Interpolate1;

  procedure Multprec_Trace_Interpolate
                ( file : in file_type; sps : in Standard_Sample_List;
                  size : in natural32;
                  t : out Multprec_Trace_Interpolators.Trace_Interpolator;
                  eps,dist,gridres,testres : out double_float ) is

    use Multprec_Stacked_Sample_Grids;
    use Multprec_Trace_Interpolators;

    timer : Timing_Widget;
    len : constant natural32 := Length_Of(sps);
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    grid : Stacked_Sample_Grid(dim,integer32(len));
    f : Floating_Number;
    testsps,testsps_last,tmp : Multprec_Sample_List;
    eva : Multprec_Complex_Numbers.Complex_Number;
    abseva : Floating_Number;
    fltabseva : double_float;

  begin
    tstart(timer);
    grid := Create_Full(file,sps,size);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid");
    tstart(timer);
    f := Maximal_Error(grid);
    eps := Round(f);
    Clear(f);
    f := Minimal_Distance(grid);
    dist := Round(f);
    Clear(f);
    tstop(timer);
    new_line(file);
    put(file,"The maximal error on samples in grid : ");
    put(file,eps,3); new_line(file);
    put(file,"The minimal distance between samples : ");
    put(file,dist,3); new_line(file);
    new_line(file);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid");
    Multprec_Create_Trace_Interpolator(file,len,grid,t,gridres);
    Sample(Head_Of(sps),5,testsps,testsps_last);
    testres := -1.0;
    tmp := testsps;
    new_line(file);
    for i in 1..integer32(5) loop
      eva := Eval(t,Sample_Point(Head_Of(tmp)).v);
      abseva := AbsVal(eva);
      put(file,"Residual of evaluation at test sample ");
      put(file,i,2); put(file," : ");
      put(file,abseva,3); new_line(file);
      fltabseva := Round(abseva);
      Clear(abseva);
      if fltabseva > testres
       then testres := fltabseva;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Generating and Evaluating at Test Samples");
  end Multprec_Trace_Interpolate;

  procedure Standard_Linear_Trace_Interpolate
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Complex_Vectors.Vector;
                   eps,dist,gridres,testres : out double_float ) is

  -- NOTE :
  --   Since a linear form requires relatively few examples, it may
  --   not be worth all the trouble going through the determination
  --   of the linear span of the component.

    dim : constant natural32 := 1;
   -- spsps,sps1 : Standard_Sample_List;
    sp : Standard_Span;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
   -- Standard_Create_Span(file,sps,dim,sp,spsps,sps1);
   -- Standard_Rectangular_Grid_Creator(file,sps1,2,grid,eps,dist);
    Standard_Rectangular_Grid_Creator(file,sps,2,grid,eps,dist);
    Standard_Create_Linear_Trace1(file,grid,t,gridres);
    Standard_Test_Linear_Evaluation(file,grid(2),dim,sp,t,testres);
    Deep_Clear(grid);
   -- Deep_Clear(spsps);
   -- if not Empty(sp)
   --  then Sampling_Machine.Clear_Restricted;
   -- end if;
  end Standard_Linear_Trace_Interpolate; 

end Irreducible_Component_Creators;
