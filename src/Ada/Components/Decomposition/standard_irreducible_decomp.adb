with unchecked_deallocation;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Linear_Projections;        use Standard_Linear_Projections;
with Standard_Polynomial_Interpolators;  use Standard_Polynomial_Interpolators;
with Standard_Breakup_Components;        use Standard_Breakup_Components;
with Filtered_Points;

package body Standard_Irreducible_Decomp is

  procedure Count_Components
                ( file : in file_type; p : in out Link_to_Poly_Sys;
                  tol : in double_float ) is

    dis : double_float;
    found : array(p'range) of boolean := (p'range => false);
    cnt : integer32 := 0;

  begin
    for i in p'range loop                        -- compute mutual distances
      if not found(i) then
        for j in i+1..p'last loop
          if not found(j) then
            if Degree(p(i)) = Degree(p(j)) then
              dis := Distance(p(i),p(j));
              put(file,"Distance between p("); put(file,i,1);
              put(file,") and p("); put(file,j,1);
              put(file,") is "); put(file,dis,3); new_line(file);
              if dis < tol
               then found(j) := true;
              end if;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    cnt := 0;                               -- count the number of components
    for i in found'range loop
      if not found(i)
       then cnt := cnt+1;
      end if;
    end loop;
    put(file,"Found "); put(file,cnt,1);
    put(file," component");
    if cnt = 1
     then put_line(file,".");
     else put_line(file,"s.");
    end if;
    declare                                         -- collect the equations
      equ : Poly_Sys(1..cnt);
      ind : integer32 := 0;
    begin
      for i in found'range loop
        if not found(i)
         then ind := ind+1; Copy(p(i),equ(ind));
        end if;
      end loop;
      Clear(p);
      p := new Poly_Sys'(equ);
    end;
  end Count_Components;

-- CREATOR :

  procedure Create_Hypersurfaces
              ( file : in file_type; n,d,itp : in natural32;
                skewproj : in boolean;
                embsys : in Poly_Sys; genpts : in Solution_List;
                surf : out Link_to_Hypersurfaces; timing : out Duration ) is

    res : Hypersurfaces(integer32(n),integer32(d));
    tol : double_float := 1.0E-8;
    timer : Timing_Widget;

  begin
    Copy(genpts,res.pts);
    for i in res.hyp'range loop
      res.hyp(i) := new Vector'(Random_Vector(0,integer32(n)));
    end loop;
    tstart(timer);
    case itp is
      when 1 =>
        res.equ := new Poly_Sys'
            (Massive_Interpolate(file,embsys,genpts,res.hyp,d));
        Count_Components(file,res.equ,tol);
      when 2 =>
        res.equ := new Poly_Sys'
            (Incremental_Interpolate(file,embsys,genpts,res.hyp,d));
      when 3 =>
        declare
          len : constant integer32 := integer32(Length_Of(genpts));
          subspaces : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(1..len);
          pivots : Standard_Integer_VecVecs.VecVec(1..len);
          basepts : Standard_Complex_VecVecs.Array_of_VecVecs(1..len);
          basecard : Standard_Natural_Vectors.Vector(1..len) := (1..len => 0);
        begin
          Dynamic_Interpolate
            (file,embsys,d,genpts,res.hyp,skewproj,subspaces,pivots,
             basepts,basecard,res.equ);
        end;
      when others => null;
    end case;
    tstop(timer);
    surf := new Hypersurfaces'(res);
    timing := Elapsed_User_Time(timer);
  end Create_Hypersurfaces;

-- AUXILIARY TO SELECTORS :

  function Zero_Component
             ( v : Vector; tol : double_float ) return integer32 is

  -- DESCRIPTION :
  --   If |v(i)| < tol, then returns i, otherwise returns v'first-1.
 
  begin
    for i in v'range loop
      if AbsVal(v(i)) < tol
       then return i;
      end if;
    end loop;
    return v'first-1;
  end Zero_Component;

-- SELECTORS :

  procedure Component_Test
              ( file : in file_type; level : in natural32;
                surf : in Hypersurfaces; tol : in double_float;
                sols : in Solution_List; solslvl : in natural32;
                num : out natural32; filpts : in out List;
                rem_first,rem_last : in out Solution_List ) is

    n : constant natural32 := Number_of_Unknowns(surf.equ(surf.equ'first));
    projsols : constant VecVec := Evaluate(surf.hyp,sols,integer32(n));
    eva : Vector(surf.equ'range);
    cnt : integer32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    ind : integer32;

  begin
    put(file,"COMPONENT TEST AT LEVEL "); put(file,level,1);
    put_line(file," :");
    for i in projsols'range loop
      ls := Head_Of(tmp);
      eva := Eval(surf.equ.all,projsols(i).all);
      put(file," eva("); put(file,i,1); put(file,") : ");
      put(file,eva,3); new_line(file);
      put(file,"Solution "); put(file,i,1);
      ind := Zero_Component(eva,tol);
      if ind >= eva'first then
        put(file," lies on component ");
        put(file,ind,1); put_line(file,".");
        cnt := cnt+1;
        Filtered_Points.Update(filpts,integer32(level),ind,integer32(solslvl));
      else
        put_line(file," fails the component test.");
        Append(rem_first,rem_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    num := natural32(cnt);
  end Component_Test;

  procedure Component_Test
              ( file : in file_type; soco : in Solution_Components;
                tol : in double_float; sols : in Solution_List;
                level : in natural32; num : out natural32; 
                filpts : in out List; rem_sols : out Solution_List ) is

    cnt,icnt : natural32 := 0;
    wrk,rem_first,rem_last : Solution_List;
    first_time : boolean := true;

  begin
    Copy(sols,wrk);
    for i in reverse 1..soco'last loop
      if soco(i) /= null then
        if first_time then
          first_time := false;
        else
          Copy(rem_first,wrk);
          Clear(rem_first); rem_last := rem_first;
        end if;
        Component_Test(file,natural32(i),soco(i).all,tol,wrk,level,
                       icnt,filpts,rem_first,rem_last);
        cnt := cnt + icnt;
      end if;
    end loop;
    num := cnt;
    rem_sols := rem_first;
  end Component_Test;

-- DESTRUCTORS :

  procedure Clear ( h : in out Hypersurfaces ) is
  begin
    Clear(h.pts);
    Clear(h.equ);
  end Clear;

  procedure Clear ( h : in out Link_to_Hypersurfaces ) is

    procedure free is
      new unchecked_deallocation(Hypersurfaces,Link_to_Hypersurfaces);

  begin
    if h /= null
     then Clear(h.all); free(h);
    end if;
  end Clear;

  procedure Clear ( s : in out Solution_Components ) is
  begin
    for i in s'range loop
      Clear(s(i));
    end loop;
  end Clear;

end Standard_Irreducible_Decomp;
