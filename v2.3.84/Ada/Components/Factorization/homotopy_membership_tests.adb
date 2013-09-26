with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Witness_Sets;                       use Witness_Sets;
with Sampling_Machine;

package body Homotopy_Membership_Tests is

-- AUXILIARY OPERATIONS :

  function Embed ( s : Solution; m : integer32 ) return Solution is

  -- DESCRIPTION :
  --   Returns the solution with m additional coordinates added to it.

    res : Solution(s.n+m);

  begin
    res.t := s.t;
    res.m := s.m;
    res.v(s.v'range) := s.v;
    for i in s.v'last+1..res.n loop
      res.v(i) := Create(0.0);
    end loop;
    return res;
  end Embed;

  function Adjusted_Slices 
                 ( sli : Standard_Complex_VecVecs.VecVec;
                   sol : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Adjusts the constant coefficient of each slice such that
  --   the solution vector satisfies the equations.

    res : Standard_Complex_VecVecs.VecVec(sli'range);

  begin
    for i in res'range loop
      res(i) := new Standard_Complex_Vectors.Vector'(sli(i).all);
      res(i)(0) := -res(i)(1)*sol(1);
      for j in 2..sol'last loop
        res(i)(0) := res(i)(0) - res(i)(j)*sol(j);
      end loop;
    end loop;
    return res;
  end Adjusted_Slices;

-- TARGET ROUTINES :

  procedure Homotopy_Membership_Test
                 ( file : in file_type;
                   genpts : in Solution_List; x : in Vector;
                   adjsli : in Standard_Complex_VecVecs.VecVec;
                   homtol : in double_float; found : out boolean ) is

    newsols,tmp : Solution_List;
    test_sol : Link_to_Solution;
    difference : Vector(x'range);
    diff : double_float;

  begin
    Sampling_Machine.Sample_with_Stop(file,genpts,x,homtol,adjsli,newsols);
    found := false;
    put_line(file,"  Scanning the new list of generic points :");
    tmp := newsols;
    for i in 1..Length_Of(newsols) loop
      test_sol := Head_Of(tmp);
      difference := test_sol.v - x;
      diff := Max_Norm(difference);
      if diff <= homtol then
        put(file,"    match with generic point "); put(file,i,1);
        put(file,", as difference is");
        put(file,diff,3); put(file," <="); put(file,homtol,3);
        new_line(file);
        found := true;
      else
        put(file,"    difference with generic point "); put(file,i,1);
        put(file," is ");
        put(file,diff,3); put(file," >"); put(file,homtol,3);
        new_line(file);
      end if;
      exit when found;
      tmp := Tail_Of(tmp);
    end loop;
    if found
     then put_line(file,"  Point lies on a solution component.");
     else put_line(file,"  Point does not lie on a solution component.");
    end if;
    Clear(newsols);
  end Homotopy_Membership_Test;

  function Embed_Solution ( s : Solution; n : integer32 ) return Solution is

  -- DESCRIPTION :
  --   Adds as many slack variables to the solution s so that its
  --   number of variables equals the ambient dimension n.

  begin
    if s.n >= n
     then return s;
     else return Embed(s,n-s.n);
    end if;
  end Embed_Solution;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Poly_Sys; dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Solution_List; s : in Solution;
                restol,homtol : in double_float;
                success,found : out boolean ) is

   -- es : Solution(s.n + dim) := Embed(s,dim);
    es : constant Solution(ep'last) := Embed_Solution(s,ep'last);
    y : constant Vector(ep'range) := Eval(ep,es.v);
    res : constant double_float := Max_Norm(y(y'first..y'last-integer32(dim)));
    newsli : Standard_Complex_VecVecs.VecVec(sli'range);
   -- eva : Complex_Number;

  begin
   -- put_line(file,"The solution vector : "); put_line(file,es.v);
   -- put_line(file,"The evaluation : "); put_line(file,y);
    put(file,"residual is "); put(file,res); new_line(file);
    if res <= restol then
      put(file,"  point satisfies the equations, as residual <=");
      put(file,restol,3); new_line(file);
      success := true;
      newsli := Adjusted_Slices(sli,es.v);
     -- for i in newsli'range loop
     --   put(file,"eval at adjusted slice : ");
     --   eva := newsli(i)(0) + newsli(i)(es.v'range)*es.v;
     --   put(file,eva); new_line(file);
     -- end loop;
      Homotopy_Membership_Test(file,genpts,es.v,newsli,homtol,found);
    else
      put(file,"  point does not lie on the component, as residual >");
      put(file,restol,3); new_line(file);
      success := false;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Poly_Sys; dim : in natural32;
                genpts,sols : in Solution_List;
                restol,homtol : in double_float ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    isols,isols_last : Solution_List;
    sli : constant Standard_Complex_VecVecs.VecVec := Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Machine.Initialize(ep);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
    put(file,"Tested "); put(file,Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not Is_Null(isols) then
      put(file,"  Found "); put(file,Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,Length_Of(isols),natural32(Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Poly_Sys; dim : in natural32;
                genpts,sols : in Solution_List;
                restol,homtol : in double_float;
                isols,isols_last : in out Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec := Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Machine.Initialize(ep);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
    put(file,"Tested "); put(file,Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not Is_Null(isols) then
      put(file,"  Found "); put(file,Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,Length_Of(isols),natural32(Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

end Homotopy_Membership_Tests;
