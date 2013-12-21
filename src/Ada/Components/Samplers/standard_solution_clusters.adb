with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;

package body Standard_Solution_Clusters is

  function Distance ( v1,v2 : Standard_Complex_Vectors.Vector )
                    return double_float is

  -- DESCRIPTION :
  --   Returns the max norm of the difference of v1 and v2.

    use Standard_Complex_Vectors;
    diff : constant vector := v1 - v2;

  begin
    return Max_Norm(diff);
  end Distance;

  procedure Separate_Clusters
               ( file : in file_type; tol : in double_float;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 nbclus : out integer32 ) is

  -- ASSUMPTION :
  --   The radius of a cluster is about 10^(-16/m), m = multiplicity.
  --   As we may have only one cluster, we take m = Length_Of(sols) in
  --   the tolerance to separate the clusters, and allow two extra
  --   decimal places for roundoff.

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(sols));
    distsols : double_float;
    tmp1 : Solution_List := Tail_Of(sols);
    ls1,ls2 : Link_to_Solution;
    tmp2 : Solution_List;
    found : boolean;

  begin
    new_line(file);
    put(file,"Separating clusters with tolerance");
    put(file,tol,3); put_line(file," ...");
    nbclus := 1;
    for i in 2..len loop                 -- classify the i-th solution
      ls1 := Head_Of(tmp1);
      tmp2 := sols;
      found := false;
      put(file,"  Checking solution "); put(file,i,1); put_line(file," :");
      for j in 1..i-1 loop
        ls2 := Head_Of(tmp2);
        distsols := Distance(ls1.v,ls2.v);
        put(file,"    Distance with solution "); put(file,j,1); 
        put(file," is"); put(file,distsols,3);
        if distsols < tol then               -- found matching cluster
          ls1.m := ls2.m;
          found := true;
          put_line(file,", is in cluster.");
        else
          tmp2 := Tail_Of(tmp2);
          put_line(file,", not in cluster.");
        end if;
        exit when found;
      end loop;
      if not found then                      -- found new cluster
        nbclus := nbclus + 1;
        ls1.m := nbclus;
        put(file,"  Found new cluster : ");
        put(file,nbclus,1); put_line(file,".");
      end if;
      Set_Head(tmp1,ls1);
      tmp1 := Tail_Of(tmp1);
    end loop;
    put(file,"Number of clusters : "); put(file,nbclus,1); new_line(file);
    new_line(file);
  end Separate_Clusters;

  function Select_Cluster
               ( sols : Standard_Complex_Solutions.Solution_List;
                 k : integer32 )
               return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if ls.m = k
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Select_Cluster;

  function Center_of_Cluster
               ( sols : Standard_Complex_Solutions.Solution_List )
               return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Vectors,Standard_Complex_Solutions;

    n : constant integer32 := Head_Of(sols).n;
    res : Vector(1..n) := Head_Of(sols).v;
    tmp : Solution_List := Tail_Of(sols);
    cnt : integer32 := 1;
    fac : double_float;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt+1;
      Add(res,Head_Of(tmp).v);
      tmp := Tail_Of(tmp);
    end loop;
    fac := 1.0/double_float(cnt);
    Mul(res,Create(fac));
    return res;
  end Center_of_Cluster;

  procedure Distances_to_Center
               ( sols : in Standard_Complex_Solutions.Solution_List;
                 center : in Standard_Complex_Vectors.Vector;
                 max,min,avg : out double_float ) is

  begin
    null;
  end Distances_to_Center;

end Standard_Solution_Clusters;
