with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Multprec_Random_Vectors;           use Multprec_Random_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vector_Tools;     use Multprec_Complex_Vector_Tools;
with Multprec_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Sample_Point_Lists_io;             use Sample_Point_Lists_io;
with Projection_Operators;              use Projection_Operators;
with Multprec_Polynomial_Interpolators; use Multprec_Polynomial_Interpolators;
with Interpolation_Filters;             use Interpolation_Filters;
with Witness_Sets_io;                   use Witness_Sets_io;
with Multiplicity_Homotopies;           use Multiplicity_Homotopies;

procedure ts_mulsam is

-- DESCRIPTION :
--   Interactive test on the multiplicity > 1 sampling machine.

  procedure Test_Multprec_Membership
               ( file : in file_type;
                 f : in Multprec_Filter; spt : in Multprec_Sample;
                 ind : in integer32; tol : in double_float ) is

  -- DESCRIPTION :
  --   Evaluates, computes the residual and membership test.

    use Multprec_Complex_Numbers;

    eva : Complex_Number;
    abseva : Floating_Number;
    flteva : double_float;

  begin
    put(file,"Residual at sample "); put(file,ind,1); put(file," : ");
    eva := Evaluate(f,Sample_Point(spt).v);
    abseva := AbsVal(eva); Clear(eva);
    flteva := Trunc(abseva); Clear(abseva);
    put(file,flteva);
    if On_Component(f,tol,Sample_Point(spt).v)
     then put_line(file," on component.");
     else put_line(file," not on component.");
    end if;
  end Test_Multprec_Membership;

  procedure Test_Multprec_Membership 
               ( file : in file_type; f : in Multprec_Filter;
                 sps : in Multprec_Sample_List; tol : in double_float ) is

  -- DESCRIPTION :
  --   Test whether the given samples belong to the component.

    tmp : Multprec_Sample_List := sps;
    spt : Multprec_Sample;

  begin
    for i in 1..Length_Of(sps) loop
      spt := Head_Of(tmp);
      Test_Multprec_Membership(file,f,spt,i,tol);
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Multprec_Membership;

  function Cluster_Radius
                   ( sols : in Multprec_Complex_Solutions.Solution_List;
                     center : in Multprec_Complex_Vectors.Vector )
                   return Floating_Number is

  -- DESCRIPTION :
  --   Returns an approximation for the radius of the cluster.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    res : Floating_Number := Create(0);
    dif : Complex_Number;
    absdif : Floating_Number;
    ls : Link_to_Solution := Head_Of(sols);

  begin
    dif := ls.v(1) - center(1);
    res := AbsVal(dif);
    Clear(dif);
    for i in 2..center'last loop
      dif := ls.v(i) - center(i);
      absdif := AbsVal(dif);
      if absdif > res
       then Copy(absdif,res);
      end if;
      Clear(dif); Clear(absdif);
    end loop;
    return res;
  end Cluster_Radius;

  procedure Center ( file : in file_type; n : in natural32;
                     hyp : in Multprec_Complex_VecVecs.VecVec;
                     sols : in Multprec_Complex_Solutions.Solution_List;
                     first,last : in out Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Adds the center of the cluster as new sample to the list first.

  -- ON ENTRY :
  --   file          for intermediate output and diagnostics;
  --   n             dimension of the ambient space;
  --   hyp           hyperplanes that cut out the generic points in sols;
  --   sols          generic points in cluster on hyperplanes;
  --   first         pointer to first element in list of sample points;
  --   last          pointer to last element in the list of sample points.

  -- ON RETURN :
  --   first         list of sample points updated with new sample;
  --   last          pointer to center point last added to the the list.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_Solutions;

    ctr : Vector(1..integer32(n));
    tmp : Solution_List := Tail_Of(sols);
    ls : Link_to_Solution;
    nb : constant natural32 := Length_Of(sols);
    cnb : Complex_Number := Create(nb);
    samsol : Solution(n);
    spt : Multprec_Sample;
    rad : Floating_Number;

  begin
    Copy(Head_Of(sols).v,ctr);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Add(ctr,ls.v);
      tmp := Tail_Of(tmp);
    end loop;
    for i in 1..n loop
      Div(ctr(i),cnb);
    end loop;
    Clear(cnb);
    put_line(file,"The center of the cluster :");
    put_line(file,ctr);
    rad := Cluster_Radius(sols,ctr);
    put(file,"The cluster radius : ");
    put(file,rad,3); new_line(file);
    Clear(rad);
    samsol.v := ctr;
    samsol.t := ls.t;
    samsol.m := ls.m;
    spt := Create(samsol,hyp);
    Append(first,last,spt);
  end Center;

  procedure Write_Summary
                   ( file : in file_type;
                     sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Writes the (err,rco,res) results for every solution in the list.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put(file,"solution "); put(file,i,1); put(file," : ");
      put(file,"  err : "); put(file,ls.err,3);
      put(file,"  rco : "); put(file,ls.rco,3);
      put(file,"  res : "); put(file,ls.res,3);
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Summary;

  procedure Sample ( file : in file_type;
                     p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in out Standard_Complex_Solutions.Solution_List;
                     eps : in double_float; nbstep,ind : in natural;
                     first,last : in out Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Does one step with tube sampling and one step with incoming homotopy.

    use Standard_Complex_Numbers;

    zeps : Complex_Number
         := Standard_Complex_Solutions.Head_Of(sols).v(p'last);
    mpsols : Multprec_Complex_Solutions.Solution_List;
    sthyp : Standard_Complex_VecVecs.VecVec(1..1);
    mphyp : Multprec_Complex_VecVecs.VecVec(1..1);

  begin
    for i in 1..nbstep loop
      Sampling_Homotopy(file,eps,p,sthyp,sols);
    end loop;
    mphyp := Create(sthyp);
    zeps := Create(REAL_PART(zeps));
    Incoming_Homotopy(file,p,sols,ind,zeps,mpsols);
    Center(file,p'last,mphyp,mpsols,first,last);
  end Sample;

  procedure Mover_Index
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 index : out integer32 ) is

  -- DESCRIPTION :
  --   Interactive determination of the index of the equation that
  --   controls the distance to the multiple component.
  --   The default is good for a special embedding to investigate curves.

    ans : integer32;

  begin
    index := p'last-1;
    loop
      put("Current moving equation is equation ");
      put(index,1); put(" with "); 
      put(Number_of_Terms(p(index)),1);
      put_line(" terms :");
      put(p(index)); new_line;
      put("Type zero to accept or choose moving equation : ");
      get(ans);
      exit when (ans = 0);
      index := ans;
    end loop;
  end Mover_Index;

  procedure Linear_Interpolate
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 dim : in natural;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Sampling and incremental interpolation after linear projections.

    ind,len,deg,nbt : natural32;
    ans : integer32;
    s,s_last : Multprec_Sample_List;
    mpi : Multprec_Projector;
    mfl : Multprec_Filter;
    tol : constant double_float := 1.0E-8;
    eps : constant double_float := 0.001;
    nbsteps : constant natural := 500;

  begin
    Mover_Index(p,ind);
    len := Standard_Complex_Solutions.Length_Of(sols);
    mpi := Create(2,p'last,len*2);
    deg := 1;
    loop
      nbt := Number_of_Terms(deg,2);
      put("Number of terms for degree ");
      put(deg,1); put(" in two variables : ");
      put(nbt,1); new_line;
      new_line;
      put("We have "); put(Length_Of(s),1); put_line(" samples.");
      put("How many new samples do you want ? (<0 to exit) ");
      get(ans);
      exit when (ans < 0);
      put("Computing "); put(ans,1);
      put_line(" new samples.  Be patient...");
      for i in 1..ans loop
        Sample(file,p,sols,eps,nbsteps,ind,s,s_last);
      end loop;
      if Length_Of(s) >= nbt then
        put("Interpolating at degree "); put(deg,1); put_line(".");
        mfl := Create(mpi);
        Sample_Update(file,mfl,s,deg);
        Test_Multprec_Membership(file,mfl,s,tol);
        Shallow_Clear(mfl);
        deg := deg+1;
      else
        put("Too few samples to interpolate at degree ");
        put(deg,1); put_line(".");
      end if;
    end loop;
    put_line(file,"The list of sample points : ");
    put(file,s);
  end Linear_Interpolate;

  procedure Central_Interpolate
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 dim,nbc : in natural32;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Sampling and incremental interpolation after central projections.

  -- ON ENTRY :
  --   file      to write intermediate output and diagnostics;
  --   p         embedded system with slack variables;
  --   dim       dimension of the solution component;
  --   nbc       number of central points used in the projection;
  --   sols      one cluster of generic points.

  -- ON RETURN :
  --   p         embedded system with modified moving equations;
  --   sols      modified cluster at last sample.

    ind,len,deg,nbt,cnt_nbc : natural32;
    ans : integer32;
    s,s_last,cp,cp_last,ptr_cp : Multprec_Sample_List;
    mpi : Multprec_Projector;
    mfl : Multprec_Filter;
    tol : constant double_float := 1.0E-8;
    eps : constant double_float := 0.001;
    nbsteps : constant natural := 500;

  begin
    Mover_Index(p,ind);
    len := Standard_Complex_Solutions.Length_Of(sols);
    mpi := Create(2,p'last,len*2);
    put("Computing "); put(nbc,1); put_line(" central points ...");
    for i in 1..nbc loop
      Sample(file,p,sols,eps,nbsteps,ind,cp,cp_last);
    end loop;
    cnt_nbc := 0;    -- number of central points used
    deg := 1;
    loop
      if (deg > 1) and (cnt_nbc < nbc) then
        put("How many central points you wish to use ? ");
        get(ans);
        if ans > integer32(nbc)
         then cnt_nbc := nbc;
         else cnt_nbc := ans;
        end if;
      end if; 
      nbt := Number_of_Terms(deg-cnt_nbc,2);
      put("Number of terms for degree ");
      put(deg,1); put(" in two variables : ");
      put(nbt,1); put(" with #central points : ");
      put(cnt_nbc,1);
      new_line;
      new_line;
      put("We have "); put(Length_Of(s),1); put_line(" samples.");
      put("How many new samples do you want ? (<0 to exit) ");
      get(ans);
      exit when (ans < 0);
      put("Computing "); put(ans,1);
      put_line(" new samples.  Be patient...");
      new_line;
      put(file,"Computing "); put(file,ans,1);
      put_line(file," new samples :");
      for i in 1..ans loop
        Sample(file,p,sols,eps,nbsteps,ind,s,s_last);
      end loop;
      if Length_Of(s) >= nbt then
        put("Interpolating at degree "); put(deg,1); put_line(".");
        mfl := Create(mpi);
        ptr_cp := cp;
        for i in 1..cnt_nbc loop
          declare
            rv : Multprec_Complex_Vectors.Vector(0..p'last)
               := Random_Vector(0,p'last,len*2);
          begin
            Central_Update(mfl,Head_Of(ptr_cp),rv);
          end;
          ptr_cp := Tail_Of(ptr_cp);
        end loop;
        Sample_Update(file,mfl,s,deg);
        Test_Multprec_Membership(file,mfl,s,tol);
        Shallow_Clear(mfl);
        deg := deg+1;
      else
        put("Too few samples to interpolate at degree ");
        put(deg,1); put_line(".");
      end if;
    end loop;
    put_line(file,"The list of sample points : ");
    put(file,s);
  end Central_Interpolate;

  procedure Main is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim,nbc : natural32;

  begin
    new_line;
    put_line("Sampling generic points from multiplicity > 1 components.");
    new_line;
    put("Reading a reconditioned system with a cluster.");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give number of central points : "); get(nbc);
    new_line;
    if nbc = 0
     then Linear_Interpolate(file,lp.all,dim,sols);
     else Central_Interpolate(file,lp.all,dim,nbc,sols);
    end if;
  end Main;

begin
  Main;
end ts_mulsam;
