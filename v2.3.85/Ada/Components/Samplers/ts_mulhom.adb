with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Solutions;
with Witness_Sets_io;                    use Witness_Sets_io;
with Standard_Solution_Clusters;         use Standard_Solution_Clusters;
with Multiplicity_Homotopies;            use Multiplicity_Homotopies;

procedure ts_mulhom is

-- DESCRIPTION :
--   Interactive development of homotopies to sample components with
--   multiplicities higher than one.

-- STEP 0 : IDENTIFICATION OF THE CLUSTERS :

  procedure Start_Cluster
               ( sols : in Standard_Complex_Solutions.Solution_List;
                 clusnb : in integer32;
                 startsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                 startsols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Creates a file with start system and corresponding start solutions
  --   for the cluster with number clusnb in sols.

    file : file_type;
    use Standard_Complex_Solutions;
    ptr_sols,ptr_start,res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    put("Reading the name of file for start of cluster ");
    put(clusnb,1); put_line(".");
    Read_Name_and_Create_File(file);
    put_line(file,startsys);
    ptr_sols := sols;
    ptr_start := startsols;
    while not Is_Null(ptr_sols) loop
      ls := Head_Of(ptr_sols);
      if ls.m = clusnb
       then Append(res,res_last,Head_Of(ptr_start).all);
      end if;
      ptr_sols := Tail_Of(ptr_sols);
      ptr_start := Tail_Of(ptr_start);
    end loop;
    new_line(file);
    put(file,"TITLE : start system for cluster ");
    put(file,clusnb,1); new_line(file);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    new_line(file);
    put(file,Length_Of(res),natural32(Head_Of(res).n),res);
    close(file);
  end Start_Cluster;

  procedure Cluster_Analysis
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   This routine compares the solutions in the list sols with each other,
  --   w.r.t. a given tolerance set by the user.

  -- ON ENTRY :
  --   file      to write diagnostics and results;
  --   p         embedded system;
  --   sols      list of solutions, possibly clustered.

  -- ON RETURN :
  --   sols      multiplicity flag identifies cluster number.

    tol : double_float := 1.0E-4;
    newtol : double_float := 0.0;
    nbclus,dim : integer32;
    startsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    startsols : Standard_Complex_Solutions.Solution_List;

  begin
    loop
      put("The current cluster tolerance is ");
      put(tol,3); new_line;
      put("Give another value, or zero to accept this tolerance : ");
      get(newtol);
      exit when (newtol = 0.0);
      tol := newtol;
    end loop;
    Separate_Clusters(file,tol,sols,nbclus);
    skip_line;
    new_line;
    put("Reading the start system and start solutions.");
    Standard_Read_Embedding(startsys,startsols,natural32(dim));
    put_line(file,"The start system");
    put_line(file,startsys.all);
    for i in 1..nbclus loop
      declare
        use Standard_Complex_Solutions;
        clussols : Solution_List := Select_Cluster(sols,i);
        clusfile : file_type;
      begin
        put(file,"SOLUTIONS in cluster "); put(file,i,1);
        put_line(file," :");
        put(file,Length_Of(clussols),natural32(Head_Of(clussols).n),clussols);
        new_line;
        put("Reading file name to write embedding + cluster ");
        put(i,1); put_line(".");
        Read_Name_and_Create_File(clusfile);
        put_line(clusfile,p);
        new_line(clusfile);
        put(clusfile,"TITLE : embedding and solutions in cluster ");
        put(clusfile,i,1); new_line(clusfile);
        new_line(clusfile);
        put_line(clusfile,"THE SOLUTIONS :");
        put(clusfile,Length_Of(clussols),
            natural32(Head_Of(clussols).n),clussols);
        Close(clusfile);
        Start_Cluster(sols,i,startsys.all,startsols);
      end;
    end loop;
  end Cluster_Analysis;

-- STEP 1 : RECONDITIONING THE CLUSTER :

  procedure Call_Reconditioning_Homotopy
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    ind : constant integer32 := p'last-1;
    zeps,a : Complex_Number;
    ans,dim,k : integer32 := 0;
    startsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    startsols : Solution_List;

  begin
    put("Reading the corresponding start solutions for the cluster.");
    Standard_Read_Embedding(startsys,startsols,natural32(dim));
    put("Give relaxation parameter in homotopy : "); get(k);
    put("Give accessibility constant in homotopy : "); get(a);
    loop
      put("Current moving equation is equation ");
      put(ind,1); put(" with "); put(Number_of_Terms(p(ind)),1);
      put_line(" terms.");
      put("Type zero to accept or choose moving equation : ");
      get(ans);
      exit when (ans = 0);
    end loop;
    Reconditioning_Homotopy(file,p,startsys.all,startsols,k,a,ind,zeps);
  end Call_Reconditioning_Homotopy;

-- STEP 2 : INCOMING HOMOTOPY TO SHARPEN THE CLUSTER :

  function Distance ( v1,v2 : Multprec_Complex_Vectors.Vector )
                    return Floating_Number is

  -- DESCRIPTION :
  --   Returns the max norm of the difference of v1 with v2.

    res : Floating_Number;
    use Multprec_Complex_Vectors;
    diff : Multprec_Complex_Vectors.Vector(v1'range) := v1 - v2;

  begin
    res := Max_Norm(diff);
    Clear(diff);
    return res;
  end Distance;

  procedure Write_Distance
              ( file : in file_type;
                sols : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Writes the distance between the solutions in the list.

    use Multprec_Complex_Solutions;
    d : Floating_Number;

  begin
    if Length_Of(sols) > 1 then
      d := Distance(Head_Of(sols).v,Head_Of(Tail_Of(sols)).v);
      put(file,"The distance is "); put(file,d,3); new_line(file);
    end if;
  end Write_Distance;

  procedure Call_Incoming_Homotopy
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Calls the homotopy that lets the solutions in the cluster move in
  --   to the singularity.

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    ind : integer32 := 0;
    zeps : constant Complex_Number := Head_Of(sols).v(p'last);
    mpsols : Multprec_Complex_Solutions.Solution_List;

  begin
    put("There are "); put(p'last,1);
    put_line(" equations in the system.");
    put("Give the index of the moving equation : "); get(ind);
    new_line;
    put_line("See the output file for results...");
    new_line;
    new_line(file);
    put(file,"The moving equation p(");
    put(file,ind,1); put(file,") = ");
    put(file,p(ind)); new_line(file);
    put(file,"Starting at zz = ");
    put(file,zeps); new_line(file);
    Incoming_Homotopy(file,p,sols,ind,zeps,mpsols);
    Write_Distance(file,mpsols);
  end Call_Incoming_Homotopy;

-- STEP 3 : THE SAMPLING HOMOTOPY :

  procedure Call_Sampling_Homotopy
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    eps : double_float := 0.0;

  begin
    put("Give the magnitude of the perturbation : "); get(eps);
    new_line;
    put_line("See the output file for results...");
    new_line;
   -- Sampling_Homotopy(file,eps,p,sols);
  end Call_Sampling_Homotopy;

  procedure Main is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    put_line(file,lp.all);
    new_line;
    put_line("MENU for homotopies to sample multiple components :");
    put_line("  0. Separate solution list into different clusters;");
    put_line("  1. Test reconditioning of the given embedding;");
    put_line("  2. Let the cluster come into the singularity;");
   -- put_line("  3. Sample from a reconditioned embedded system.");
   -- put("Type 0, 1, 2, or 3 to select : ");
    put("Type 0, 1, or 2 to select : ");
   -- Ask_Alternative(ans,"0123");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Cluster_Analysis(file,lp.all,sols);
      when '1' => Call_Reconditioning_Homotopy(file,lp.all,sols);
      when '2' => Call_Incoming_Homotopy(file,lp.all,sols);
     -- when '3' => Call_Sampling_Homotopy(file,lp.all,sols);
      when others => null;
    end case;
  end Main;

begin
  new_line;
  put_line("Homotopies to Sample high Multiplicity Components.");
  Main;
end ts_mulhom;
