with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Sampling_Machine;
with Interpolation_Filters;              use Interpolation_Filters;
with Span_of_Component;                  use Span_of_Component;
with Irreducible_Components;             use Irreducible_Components;
with Main_Poly_Continuation;             use Main_Poly_Continuation;

package body Drivers_to_Component_Creators is

  procedure Display_Filter_Settings
               ( file : in file_type; full_output : in boolean;
                 stoptol,membtol : in double_float ) is
  begin
    new_line(file);
    put_line(file,"Current settings in the filter creator :");
    put(file,"  1. tolerance for stop test  : "); put(file,stoptol,3);
    new_line(file);
    put(file,"  2. tolerance for membership : "); put(file,membtol,3);
    new_line(file);
    put(file,"  3. full intermediate output : ");
    if full_output
     then put_line(file,"  yes");
     else put_line(file,"  no, only summary");
    end if;
  end Display_Filter_Settings;

  procedure Standard_Tuner
              ( file : in file_type; full_output : in out boolean;
                stoptol,membtol : in out double_float ) is

    ans : character;
    oc : natural32;

  begin
    loop
      Display_Filter_Settings(Standard_Output,full_output,stoptol,membtol);
      put("Type 0 to exit, 1,2 or 3 to change : ");
      Ask_Alternative(ans,"0123");
      exit when (ans = '0');
      case ans is
        when '1' => put("Give new tolerance for stop test : ");
                    get(stoptol);
        when '2' => put("Give new tolerance for membership : ");
                    get(membtol);
        when '3' => full_output := not full_output;
        when others => null;
      end case;
    end loop;
    Display_Filter_Settings(file,full_output,stoptol,membtol);
    new_line;
    Sampling_Machine.Interactive_Tune_Sampler(file);
   -- Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Refiner(file);
    if full_output
     then new_line;
          Driver_for_Process_io(file,oc);
    end if;
  end Standard_Tuner;

  procedure Multprec_Tuner
              ( file : in file_type; full_output : in out boolean;
                size : in natural32;
                stoptol,membtol : in out double_float ) is

    ans : character;
    oc : natural32;

  begin
    loop
      Display_Filter_Settings(Standard_Output,full_output,stoptol,membtol);
      put("Type 0 to exit, 1,2 or 3 to change : ");
      Ask_Alternative(ans,"0123");
      exit when (ans = '0');
      case ans is
        when '1' => put("Give new tolerance for stop test : ");
                    get(stoptol);
        when '2' => put("Give new tolerance for membership : ");
                    get(membtol);
        when '3' => full_output := not full_output;
        when others => null;
      end case;
    end loop;
    Display_Filter_Settings(file,full_output,stoptol,membtol);
    new_line;
    Sampling_Machine.Interactive_Tune_Sampler(file);
   -- Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(file,size);
    if full_output
     then new_line;
          Driver_for_Process_io(file,oc);
    end if;
  end Multprec_Tuner;

  function Cardinality ( v : Standard_Natural_Vectors.Vector )
                       return natural32 is

  -- DESCRIPTION :
  --   Returns the number of nonzero entries in the vector v.

    res : natural32 := 0;

  begin
    for i in v'range loop
      if v(i) = 0
       then return res;
       else res := res+1;
      end if;
    end loop;
    return res;
  end Cardinality;

  procedure Write_Labels ( file : in file_type;
                           g : in Standard_Natural_Vectors.Link_to_Vector ) is

    use Standard_Natural_Vectors;
    card : natural32;

  begin
    if g /= null then
      card := Cardinality(g.all);
      put(file,"    contains ");
      put(file,card,1);
      if card = 1
       then put(file," generic point");
       else put(file," generic points");
      end if;
      for i in 1..card loop
        put(file," ");
        put(file,g(integer32(i)),1);
      end loop;
      put_line(file,".");
    end if;
  end Write_Labels;

  procedure Write_Summary
               ( file : in file_type; k : in natural32;
                 deco : in Standard_Irreducible_Component_List ) is

    len : constant natural32 := Length_Of(deco);
    tmp : Standard_Irreducible_Component_List := deco;
    sum : natural32 := 0;

  begin
    put(file,"Found "); put(file,len,1);
    if len = 1
     then put(file," irreducible component");
     else put(file," irreducible components");
    end if; 
    put(file," of dimension ");
    put(file,k,1);
    put_line(file," :");
    for i in 1..len loop
      put(file,"  Component "); put(file,i,1);
      declare
        c : constant Standard_Irreducible_Component := Head_Of(tmp);
        dc : constant natural32 := Degree(c);
        f : constant Standard_Filter := Filter(c);
        df : constant natural32 := Degree(f);
        s : constant Standard_Span := Span(c);
      begin
        put(file," has degree "); put(file,dc,1);
        if df > 0 and df < dc then
          put(file,"(="); put(file,df,1);
          put(file,"+"); put(file,dc-df,1);
          put(file,")");
        end if;
        if Empty(s) then
          put_line(file,";");
        else
          put(file," and spans space of dimension ");
          put(file,Dimension(s),1); put_line(file,";");
        end if;
        Write_Labels(file,Labels(c));
        sum := sum + dc;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"The sum of the degrees is "); put(file,sum,1);
    put_line(file,".");
  end Write_Summary;

  procedure Write_Summary
               ( file : in file_type; k : in natural32;
                 deco : in Multprec_Irreducible_Component_List ) is

    len : constant natural32 := Length_Of(deco);
    tmp : Multprec_Irreducible_Component_List := deco;
    sum : natural32 := 0;

  begin
    put(file,"Found "); put(file,len,1);
    if len = 1
     then put(file," irreducible component");
     else put(file," irreducible components");
    end if;
    put(file," of dimension ");
    put(file,k,1);
    put_line(file," :");
    for i in 1..len loop
      put(file,"  Component "); put(file,i,1);
      declare
        c : constant Multprec_Irreducible_Component := Head_Of(tmp);
        dc : constant natural32 := Degree(c);
        f : constant Multprec_Filter := Filter(c);
        df : constant natural32 := Degree(f);
        s : constant Multprec_Span := Span(c);
      begin
        put(file," has degree "); put(file,dc,1);
        if df > 0 and df < dc then
          put(file,"(="); put(file,df,1);
          put(file,"+"); put(file,dc-df,1);
          put(file,")");
        end if;
        if Empty(s) then
          put_line(file,";");
        else
          put(file," and spans space of dimension ");
          put(file,Dimension(s),1); put_line(file,";");
        end if;
        Write_Labels(file,Labels(c));
        sum := sum + dc;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"The sum of the degrees is "); put(file,sum,1);
    put_line(file,".");
  end Write_Summary;

end Drivers_to_Component_Creators;
