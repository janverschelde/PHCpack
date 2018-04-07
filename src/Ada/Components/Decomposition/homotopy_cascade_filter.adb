with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Multprec_Complex_Polynomials;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;          use Standard_Solution_Filters;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with Witness_Sets;                       use Witness_Sets;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with Filtered_Points;

package body Homotopy_Cascade_Filter is

-- AUXILIARY :

  procedure Down_Continuation
              ( file : in file_type; 
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in integer32; sols : in out Solution_List;
                pocotime : out duration ) is

  -- DESCRIPTION :
  --   Performs a continuation to remove the slice from the embedded system.
  --   On entry, sols contains the start solutions, on return, the
  --   computed solutions are in the list sols.

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(0.0));
    Black_Box_Polynomial_Continuation(file,true,target,embsys,sols,pocotime);
  end Down_Continuation;

-- TARGET ROUTINES :

  procedure Standard_Initialize
              ( soco : in out Standard_Irreducible_Decomp.Solution_Components;
                n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    surf : Standard_Irreducible_Decomp.Hypersurfaces(n,0);

  begin
    surf.equ := new Standard_Complex_Poly_Systems.Poly_Sys'(p);
    soco(0) := new Standard_Irreducible_Decomp.Hypersurfaces'(surf);
  end Standard_Initialize;

  procedure Multprec_Initialize
              ( soco : in out Multprec_Irreducible_Decomp.Solution_Components;
                n : in integer32;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    surf : Multprec_Irreducible_Decomp.Hypersurfaces(n,0);

  begin
    surf.equ := new Multprec_Complex_Poly_Systems.Poly_Sys'(p);
    soco(0) := new Multprec_Irreducible_Decomp.Hypersurfaces'(surf);
  end Multprec_Initialize;

  procedure Standard_Update_Filter
              ( fp,fp_last : in out List; n,k : in integer32;
                soco : in Standard_Irreducible_Decomp.Solution_Components ) is

    deg : integer32;
    use Standard_Complex_Polynomials;

  begin
    for i in soco(k).equ'range loop
      if soco(k).equ(i) /= Null_Poly then
        deg := Standard_Complex_Polynomials.Degree(soco(k).equ(i));
        Append(fp,fp_last,Filtered_Points.Create(n,k,deg));
      end if;
    end loop;
  end Standard_Update_Filter;

  procedure Multprec_Update_Filter
              ( fp,fp_last : in out List; n,k : in integer32;
                soco : in Multprec_Irreducible_Decomp.Solution_Components ) is

    deg : integer32;
    use Multprec_Complex_Polynomials;

  begin
    for i in soco(k).equ'range loop
      if soco(k).equ(i) /= Null_Poly then
        deg := Multprec_Complex_Polynomials.Degree(soco(k).equ(i));
        Append(fp,fp_last,Filtered_Points.Create(n,k,deg));
      end if;
    end loop;
  end Multprec_Update_Filter;

  procedure Standard_Update_Hypersurfaces
              ( file : in file_type;
                soco : in out Standard_Irreducible_Decomp.Solution_Components;
                n,top,k,itp : in natural32; skewproj : in boolean;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Solution_List; tol_sing : in double_float;
                timings : out Duration; fp,fp_last : in out List ) is

    reg_sols : constant Solution_List := Regular_Filter(sols,tol_sing);

  begin
    if not Is_Null(reg_sols) then
      Standard_Irreducible_Decomp.Create_Hypersurfaces
        (file,n+k,k,itp,skewproj,embsys,reg_sols,soco(integer32(k)),timings);
      Standard_Update_Filter(fp,fp_last,integer32(top)+1,integer32(k),soco);
    else
      timings := 0.0;
    end if;
  end Standard_Update_Hypersurfaces;

  procedure Multprec_Update_Hypersurfaces
              ( file : in file_type;
                soco : in out Multprec_Irreducible_Decomp.Solution_Components;
                n,top,k,size,itp : in natural32; skewproj : in boolean;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Solution_List; tol_sing : in double_float;
                timings : out Duration; fp,fp_last : in out List ) is

    reg_sols : constant Solution_List := Regular_Filter(sols,tol_sing);

  begin
    if not Is_Null(reg_sols) then
      Multprec_Irreducible_Decomp.Create_Hypersurfaces
        (file,n+k,k,size,itp,skewproj,embsys,orgsys,reg_sols,
         soco(integer32(k)),timings);
      Multprec_Update_Filter(fp,fp_last,integer32(top)+1,integer32(k),soco);
    else
      timings := 0.0;
    end if;
  end Multprec_Update_Hypersurfaces;

  procedure Standard_Cascade_Loop
              ( file : in file_type; n,k,itp : in natural32;
                skewproj : in boolean;
                embp : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                sols : in out Solution_List;
                soco : in out Standard_Irreducible_Decomp.Solution_Components;
                tol,tol_sing,tol_eval : in double_float;
                npa,ns0,ns1,div : in out Standard_Natural_Vectors.Vector;
                fp,fp_last : in out List;
                gentimes,clatimes : in out Array_of_Duration ) is

    firstzero : boolean := (ns0(integer32(k)) /= 0);
    sols0,sols1,rem_sols : Solution_List;
    cnt : natural32;

    use Standard_Irreducible_Decomp;

  -- ALGORITHM :
  --   The classification is organized in three cases :

    procedure Top_Component ( i : in integer32 ) is

    -- DESCRIPTION :
    --   Here were are for the first time a component of solutions.

    begin
      firstzero := (ns0(i-1) /= 0);
      if (firstzero and (i > 1)) then
        Remove_Component(sols0);
        Standard_Update_Hypersurfaces
          (file,soco,n,k,natural32(i-1),itp,skewproj,embp(i-1).all,sols0,
           tol_sing,clatimes(integer(i)-1),fp,fp_last);
      end if;
    end Top_Component;

    procedure General_Component ( i : in integer32 ) is

    -- DESCRIPTION :
    --   In this case we have already one component and found
    --   also solution with zero slack variables.

      comptesttimer : Timing_Widget;

    begin
      tstart(comptesttimer);
      Component_Test(file,soco,tol_eval,sols0,natural32(i-1),cnt,fp,rem_sols);
      tstop(comptesttimer);
      if cnt /= Length_Of(sols0) then
        put(file,"WARNING : only "); put(file,cnt,1);
        put(file," out of "); put(file,ns0(i-1),1);
        put_line(file," solutions could be classified.");
        put_line(file,"The unclassified solutions : ");
        put(file,Length_Of(rem_sols),natural32(Head_Of(rem_sols).n),rem_sols);
        Remove_Component(rem_sols);
        Standard_Update_Hypersurfaces
          (file,soco,n,k,natural32(i-1),itp,skewproj,embp(i-1).all,rem_sols,
           tol_sing,clatimes(integer(i)-1),fp,fp_last);
      end if;
      new_line(file);
      print_times(file,comptesttimer,"Component Membership Test");
      new_line(file);
      clatimes(integer(i)-1)
        := clatimes(integer(i)-1) + Elapsed_User_Time(comptesttimer);
    end General_Component;

    procedure Bottom_Component ( i : in integer32 ) is

    -- DESCRIPTION :
    --   This is the final classification stage at the end of the cascade,
    --   where we may discover actual isolated solutions.

      comptesttimer : Timing_Widget;

    begin
      tstart(comptesttimer);
      Component_Test(file,soco,tol_eval,sols1,natural32(i-1),cnt,fp,rem_sols);
      tstop(comptesttimer);
      ns0(i-1) := ns0(i-1) + cnt;
      ns1(i-1) := ns1(i-1) - cnt;
      if Length_Of(rem_sols) = 0 then
        put_line(file,"NO ISOLATED SOLUTIONS");
      else
        put_line(file,"ISOLATED SOLUTIONS : ");
        put(file,Length_Of(rem_sols),natural32(Head_Of(rem_sols).n),rem_sols);
        Remove_Component(rem_sols);
        Copy(rem_sols,soco(0).pts);
        Append(fp,fp_last,
               Filtered_Points.Create(integer32(k+1),0,
                                      integer32(Length_Of(rem_sols))));
      end if;
      new_line(file);
      print_times(file,comptesttimer,"Component Membership Test");
      new_line(file);
      clatimes(0) := clatimes(0) + Elapsed_User_Time(comptesttimer);
    end Bottom_Component;

  begin
    for i in reverse 1..integer32(k) loop
      npa(i-1) := Length_Of(sols); 
      Down_Continuation(file,embp(i).all,i,sols,gentimes(integer(i)-1));
      Clear(sols0); Clear(sols1);
      Filter_and_Split_Solutions
        (file,sols,integer32(n),i-1,tol,sols0,sols1);
      ns0(i-1) := Length_Of(sols0); 
      ns1(i-1) := Length_Of(sols1); 
      div(i-1) := npa(i-1) - ns0(i-1) - ns1(i-1);
      if not firstzero then
        Top_Component(i);
      elsif ns0(i-1) > 0 then
        General_Component(i);
      elsif i = 1 and ns1(i-1) > 0 then
        Bottom_Component(i);
      end if;
      Clear(sols);
      exit when Is_Null(sols1);
      sols := Remove_Component(sols1);
    end loop;
  end Standard_Cascade_Loop;

  procedure Multprec_Cascade_Loop
              ( file : in file_type; n,k,size,itp : in natural32;
                skewproj : in boolean;
                embp : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in out Solution_List;
                soco : in out Multprec_Irreducible_Decomp.Solution_Components;
                tol,tol_sing,tol_eval : in double_float;
                npa,ns0,ns1,div : in out Standard_Natural_Vectors.Vector;
                fp,fp_last : in out List;
                gentimes,clatimes : in out Array_of_Duration ) is

    firstzero : boolean := (ns0(integer32(k)) /= 0);
    sols0,sols1,rem_sols : Solution_List;
    cnt : natural32;

    use Multprec_Irreducible_Decomp;

  -- ALGORITHM :
  --   The classification is organized in three cases :

    procedure Top_Component ( i : in integer32) is

    -- DESCRIPTION :
    --   Here were are for the first time a component of solutions.

    begin
      firstzero := (ns0(i-1) /= 0);
      if (firstzero and (i > 1)) then
        Remove_Component(sols0);
        Multprec_Update_Hypersurfaces
          (file,soco,n,k,natural32(i-1),size,itp,skewproj,embp(i-1).all,orgsys,
           sols0,tol_sing,clatimes(integer(i)-1),fp,fp_last);
      end if;
    end Top_Component;

    procedure General_Component ( i : in integer32 ) is

    -- DESCRIPTION :
    --   In this case we have already one component and found
    --   also solution with zero slack variables.

      comptesttimer : Timing_Widget;

    begin
      tstart(comptesttimer);
      Component_Test(file,soco,tol_eval,sols0,natural32(i-1),cnt,fp,rem_sols);
      tstop(comptesttimer);
      if cnt /= Length_Of(sols0) then
        put(file,"WARNING : only "); put(file,cnt,1);
        put(file," out of "); put(file,ns0(i-1),1);
        put_line(file," solutions could be classified.");
        put_line(file,"The unclassified solutions : ");
        put(file,Length_Of(rem_sols),natural32(Head_Of(rem_sols).n),rem_sols);
        Remove_Component(rem_sols);
        Multprec_Update_Hypersurfaces
          (file,soco,n,k,natural32(i-1),size,itp,skewproj,embp(i-1).all,orgsys,
           rem_sols,tol_sing,clatimes(integer(i)-1),fp,fp_last);
      end if;
      new_line(file);
      print_times(file,comptesttimer,"Component Membership Test");
      new_line(file);
      clatimes(integer(i)-1)
        := clatimes(integer(i)-1) + Elapsed_User_Time(comptesttimer);
    end General_Component;

    procedure Bottom_Component ( i : in integer32 ) is

    -- DESCRIPTION :
    --   This is the final classification stage at the end of the cascade,
    --   where we may discover actual isolated solutions.

      comptesttimer : Timing_Widget;

    begin
      tstart(comptesttimer);
      Component_Test(file,soco,tol_eval,sols1,natural32(i-1),cnt,fp,rem_sols);
      tstop(comptesttimer);
      ns0(i-1) := ns0(i-1) + cnt;
      ns1(i-1) := ns1(i-1) - cnt;
      if Length_Of(rem_sols) = 0 then
        put_line(file,"NO ISOLATED SOLUTIONS");
      else
        put_line(file,"ISOLATED SOLUTIONS : ");
        put(file,Length_Of(rem_sols),natural32(Head_Of(rem_sols).n),rem_sols);
        Remove_Component(rem_sols);
        Copy(rem_sols,soco(0).pts);
        Append(fp,fp_last,
               Filtered_Points.Create(integer32(k)+1,0,
                                      integer32(Length_Of(rem_sols))));
      end if;
      new_line(file);
      print_times(file,comptesttimer,"Component Membership Test");
      new_line(file);
      clatimes(0) := clatimes(0) + Elapsed_User_Time(comptesttimer);
    end Bottom_Component;

  begin
    for i in reverse 1..integer32(k) loop
      npa(i-1) := Length_Of(sols); 
      Down_Continuation(file,embp(i).all,i,sols,gentimes(integer(i)-1));
      Clear(sols0); Clear(sols1);
      Filter_and_Split_Solutions
        (file,sols,integer32(n),i-1,tol,sols0,sols1);
      ns0(i-1) := Length_Of(sols0); 
      ns1(i-1) := Length_Of(sols1); 
      div(i-1) := npa(i-1) - ns0(i-1) - ns1(i-1);
      if not firstzero then
        Top_Component(i);
      elsif ns0(i-1) > 0 then
        General_Component(i);
      elsif i = 1 and ns1(i-1) > 0 then
        Bottom_Component(i);
      end if;
      Clear(sols);
      exit when Is_Null(sols1);
      sols := Remove_Component(sols1);
    end loop;
  end Multprec_Cascade_Loop;

end Homotopy_Cascade_Filter;
