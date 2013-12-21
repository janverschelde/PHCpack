with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Multprec_Complex_VecVecs;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Sampling_Machine;
with Sample_Points,Sample_Point_Lists;   use Sample_Points,Sample_Point_Lists;
with Standard_Polynomial_Interpolators;  use Standard_Polynomial_Interpolators;
with Multprec_Polynomial_Interpolators;  use Multprec_Polynomial_Interpolators;
with Projection_Operators;               use Projection_Operators;
with Interpolation_Filters;              use Interpolation_Filters;

procedure ts_filter is

  procedure Test_Standard_Membership
               ( f : in Standard_Filter; spt : in Standard_Sample;
                 ind : in natural32; tol : in double_float ) is

  -- DESCRIPTION :
  --   Evaluates, computes the residual and membership test.

  begin
    put("Residual at sample "); put(ind,1); put(" : ");
    put(AbsVal(Evaluate(f,Sample_Point(spt).v)));
    if On_Component(f,tol,Sample_Point(spt).v)
     then put_line(" on component.");
     else put_line(" not on component.");
    end if;
  end Test_Standard_Membership;

  procedure Test_Multprec_Membership
               ( f : in Multprec_Filter; spt : in Multprec_Sample;
                 ind : in natural32; tol : in double_float ) is

  -- DESCRIPTION :
  --   Evaluates, computes the residual and membership test.

    eva : Multprec_Complex_Numbers.Complex_Number;
    abseva : Floating_Number;
    flteva : double_float;

  begin
    put("Residual at sample "); put(ind,1); put(" : ");
    eva := Evaluate(f,Sample_Point(spt).v);
    abseva := AbsVal(eva); Clear(eva);
    flteva := Trunc(abseva); Clear(abseva);
    put(flteva);
    if On_Component(f,tol,Sample_Point(spt).v)
     then put_line(" on component.");
     else put_line(" not on component.");
    end if;
  end Test_Multprec_Membership;

  procedure Test_Standard_Membership 
               ( f : in Standard_Filter; sps : in Standard_Sample_List;
                 tol : in double_float ) is

  -- DESCRIPTION :
  --   Test whether the given samples belong to the component.

    tmp : Standard_Sample_List := sps;
    spt : Standard_Sample;

  begin
    for i in 1..Length_Of(sps) loop
      spt := Head_Of(tmp);
      Test_Standard_Membership(f,spt,i,tol);
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Standard_Membership;

  procedure Test_Multprec_Membership 
               ( f : in Multprec_Filter; sps : in Multprec_Sample_List;
                 tol : in double_float ) is

  -- DESCRIPTION :
  --   Test whether the given samples belong to the component.

    tmp : Multprec_Sample_List := sps;
    spt : Multprec_Sample;

  begin
    for i in 1..Length_Of(sps) loop
      spt := Head_Of(tmp);
      Test_Multprec_Membership(f,spt,i,tol);
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Multprec_Membership;

  procedure Test_Standard_Sample
               ( f : in Standard_Filter; spt : in Standard_Sample;
                 tol : in double_float; ind : in natural32;
                 newspt : out Standard_Sample ) is

  -- DESCRIPTION :
  --   Computes a new sample and tests whether it satisfies the filter.

  begin
    put_line("Computing one new sample...");
    Sample(spt,newspt);
    Test_Standard_Membership(f,newspt,ind,tol);
  end Test_Standard_Sample;

  procedure Test_Multprec_Sample
               ( f : in Multprec_Filter; spt : in Standard_Sample;
                 tol : in double_float; ind : in natural32;
                 newspt : out Multprec_Sample ) is

  -- DESCRIPTION :
  --   Computes a new sample and tests whether it satisfies the filter.

  begin
    put_line("Computing one new sample...");
    Sample(spt,newspt);
    Test_Multprec_Membership(f,newspt,ind,tol);
  end Test_Multprec_Sample;

  procedure Test_Standard_Update
               ( f : in out Standard_Filter; sps : in Standard_Sample_List;
                 d : in natural32 ) is

  -- DESCRIPTION :
  --   Updates the filter with new samples to arrive at degree d.

  begin
    Sample_Update(Standard_Output,f,sps,d);
    put("Dimension : "); put(Dimension(f),1);
    put("  Degree : "); put(Degree(f),1);
    put("  Centrality : "); put(Centrality(f),1);
    put_line("  Interpolator : "); put_line(Interpolator(f));
  end Test_Standard_Update;

  procedure Test_Multprec_Update
               ( f : in out Multprec_Filter; sps : in Multprec_Sample_List;
                 d : in natural32 ) is

  -- DESCRIPTION :
  --   Updates the filter with new samples to arrive at degree d.

  begin
    Sample_Update(Standard_Output,f,sps,d);
    put("Dimension : "); put(Dimension(f),1);
    put("  Degree : "); put(Degree(f),1);
    put("  Centrality : "); put(Centrality(f),1);
    put_line("  Interpolator : "); put_line(Interpolator(f));
  end Test_Multprec_Update;

  procedure Standard_Sample_Central_Points
               ( f : in out Standard_Filter; spt : in Standard_Sample;
                 ncp : in natural32 ) is

    n : constant integer32 := Sample_Point(spt).n;

  begin
    put("Sampling "); put(ncp,1); put(" central points...");
    for i in 1..ncp loop
      declare
        newspt : Standard_Sample;
        ranhyp : Standard_Complex_Vectors.Vector(0..n);
      begin
        Sample(spt,newspt);
        ranhyp := Random_Vector(0,n);
        Central_Update(f,newspt,ranhyp);
      end;
    end loop;                        
    put("  Centrality of the filter : "); put(Centrality(f),1); new_line;
  end Standard_Sample_Central_Points;

  procedure Multprec_Sample_Central_Points
               ( f : in out Multprec_Filter; spt : in Standard_Sample;
                 ncp : in natural32 ) is

    n : constant integer32 := Sample_Point(spt).n;

  begin
    put("Sampling "); put(ncp,1); put(" central points...");
    for i in 1..ncp loop
      declare
        newspt : Multprec_Sample;
        stranhyp : Standard_Complex_Vectors.Vector(0..n);
        mpranhyp : Multprec_Complex_Vectors.Vector(0..n);
      begin
        Sample(spt,newspt);
        stranhyp := Random_Vector(0,n);
        mpranhyp := Create(stranhyp);
        Central_Update(f,newspt,mpranhyp);
      end;
    end loop;                        
    put("  Centrality of the filter : "); put(Centrality(f),1); new_line;
  end Multprec_Sample_Central_Points;

  procedure Test_Standard_Construction
               ( sol : in Standard_Complex_Solutions.Solution;
                 slihyp,prohyp : in Standard_Complex_VecVecs.VecVec;
                 tol : in double_float ) is

  -- DESCRIPTION :
  --   Allows to build the interpolating filter for increasing degrees.

    p : constant Standard_Projector := Create(prohyp);
    f : Standard_Filter := Create(p);
    dim : constant natural32 := natural32(slihyp'last);
    spt : constant Standard_Sample := Create(sol,slihyp);
    newspt : Standard_Sample;
    d,nb,cnt,ncp : natural32 := 0;
    sps,sps_last,prev : Standard_Sample_List;
    ans : character;

  begin
    put("Give the degree : "); get(d);
    put("Give the number of central points : "); get(ncp);
    Standard_Sample_Central_Points(f,spt,ncp);
    d := d-ncp;
    nb := Standard_Polynomial_Interpolators.Number_of_Terms(d,dim+1) - 1;
    put("Sampling "); put(nb,1); put_line(" points...");
    Sample(spt,nb,sps,sps_last);
    Test_Standard_Update(f,sps,d);
    loop
      prev := sps_last;
      Test_Standard_Membership(f,sps,tol);
      cnt := Length_Of(sps);
      loop
        cnt := cnt+1;
        Test_Standard_Sample(f,spt,tol,cnt,newspt);
        Append(sps,sps_last,newspt);
        put("Do you want to test more samples ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end loop;
      put("Give higher degree or 0 to exit : "); get(d);
      exit when (d = 0);
      put("Give the number of central points : "); get(ncp);
      Standard_Sample_Central_Points(f,spt,ncp);
      d := d-Centrality(f);
      nb := Standard_Polynomial_Interpolators.Number_of_Terms(d,dim+1) - 1;
      if Length_Of(sps) < nb then
        nb := nb - Length_Of(sps);
        put("Sampling "); put(nb,1); put_line(" additional points...");
        Sample(spt,nb,sps,sps_last);
      else
        put_line("No new samples needed.");
      end if;
      prev := Tail_Of(prev);
      Test_Standard_Update(f,prev,d);
    end loop;
  end Test_Standard_Construction;

  procedure Test_Multprec_Construction
               ( sol : in Standard_Complex_Solutions.Solution;
                 slihyp : in Standard_Complex_VecVecs.VecVec;
                 prohyp : in Multprec_Complex_VecVecs.VecVec;
                 tol : in double_float ) is

  -- DESCRIPTION :
  --   Allows to build the interpolating filter for increasing degrees.

    p : constant Multprec_Projector := Create(prohyp);
    f : Multprec_Filter := Create(p);
    dim : constant natural32 := natural32(slihyp'last);
    spt : constant Standard_Sample := Create(sol,slihyp);
    newspt : Multprec_Sample;
    d,nb,cnt,ncp : natural32 := 0;
    sps,sps_last,prev : Multprec_Sample_List;
    ans : character;

  begin
    put("Give the degree : "); get(d);
    put("Give the number of central points : "); get(ncp);
    Multprec_Sample_Central_Points(f,spt,ncp);
    d := d - ncp;
    nb := Multprec_Polynomial_Interpolators.Number_of_Terms(d,dim+1) - 1;
    put("Sampling "); put(nb,1); put_line(" points...");
    Sample(spt,nb,sps,sps_last);
    Test_Multprec_Update(f,sps,d);
    loop
      prev := sps_last;
      Test_Multprec_Membership(f,sps,tol);
      cnt := Length_Of(sps);
      loop
        cnt := cnt+1;
        Test_Multprec_Sample(f,spt,tol,cnt,newspt);
        Append(sps,sps_last,newspt);
        put("Do you want to test more samples ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end loop;
      put("Give higher degree or 0 to exit : "); get(d);
      exit when (d = 0);
      put("Give the number of central points : "); get(ncp);
      Multprec_Sample_Central_Points(f,spt,ncp);
      d := d-Centrality(f);
      nb := Multprec_Polynomial_Interpolators.Number_of_Terms(d,dim+1) - 1;
      if Length_Of(sps) < nb then
        nb := nb - Length_Of(sps);
        put("Sampling "); put(nb,1); put_line(" additional points...");
        Sample(spt,nb,sps,sps_last);
      else
        put_line("No new samples needed.");
      end if;
      prev := Tail_Of(prev);
      Test_Multprec_Update(f,prev,d);
    end loop;
  end Test_Multprec_Construction;

  procedure Test_Standard_Filter
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim : in natural32 ) is

  -- DESCRIPTION :
  --   Runs through the list of generic points and tests the creation
  --   of the interpolation filters.  

    tol : constant double_float := 1.0E-8;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    prohyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim)+1)
           := Random_Hyperplanes(dim+1,natural32(p'last));
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Length_Of(sols) loop
      new_line;
      put("For sample no "); put(i,1); put(", ");
      Test_Standard_Construction(Head_Of(tmp).all,hyp,prohyp,tol);
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Test_Standard_Filter;
  
  procedure Test_Multprec_Filter
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim : in natural32 ) is

   -- f : Multprec_Filter;
    tol : double_float;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    size,deci : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    prohyp : Multprec_Complex_VecVecs.VecVec(1..integer32(dim)+1);
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci);
    tol := 10.0**integer(-deci/2);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,dim);
    prohyp := Random_Hyperplanes(dim+1,natural32(p'last),size);
    Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner(size);
    for i in 1..Length_Of(sols) loop
      new_line;
      put("For sample no "); put(i,1); put(", ");
      Test_Multprec_Construction(Head_Of(tmp).all,hyp,prohyp,tol);
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Test_Multprec_Filter;

  procedure Main is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    new_line;
    put_line("Testing the construction of filters for solution components");
    Standard_Read_Embedding(lp,sols,dim);
    loop
      new_line;
      put_line("MENU for the construction of component filters : ");
      put_line("  0. exit this program;");
      put_line("  1. test construction with standard arithmetic;");
      put_line("  2. test construction with multi-precision arithmetic");
      put("Type 0,1, or 2 to select : "); Ask_Alternative(ans,"012");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Standard_Filter(lp.all,sols,dim);
        when '2' => Test_Multprec_Filter(lp.all,sols,dim);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_filter;
