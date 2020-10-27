with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Integer_Lifting_Functions;          use Integer_Lifting_Functions;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Floating_Lifting_Functions;         use Floating_Lifting_Functions;
with Floating_Lifting_Utilities;         use Floating_Lifting_Utilities;
with Floating_Integer_Convertors;
with Floating_Mixed_Subdivisions_io;

package body Main_Lifting_Functions is

-- AUXILIARIES :

  procedure put ( file : in file_type;
                  v : in Standard_Floating_Vectors.Vector;
                  fore,aft,exp : in natural32 ) is
  begin
    for i in v'range loop
      put(file,v(i),fore,aft,exp);
    end loop;
  end put;

  procedure put ( v : in Standard_Floating_Vectors.Vector;
                  fore,aft,exp : in natural32 ) is
  begin
    put(Standard_Output,v,fore,aft,exp);
  end put;

  procedure Read_Integer_Vector
              ( v : in out Standard_Integer_Vectors.Vector ) is
  begin
    for i in v'range loop
      put("Give integer for component "); put(i,1); put(" : ");
      Read_Integer(v(i));
    end loop;
  end Read_Integer_Vector;

  procedure Read_Float_Vector
              ( v : in out Standard_Floating_Vectors.Vector ) is
  begin
    for i in v'range loop
      put("Give float for component "); put(i,1); put(" : ");
      Read_Double_Float(v(i));
    end loop;
  end Read_Float_Vector;

  procedure Set_Integer_Bounds
              ( file : in file_type;
                low,upp : in out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows the default values for lower and upper bounds and allows
  --   the user to modify these.

    ans : character;
    m : constant integer32 := low'length;

  begin
    loop
      new_line;
      put_line("Current lower and upper bounds on lifting values");
      put("  lower : "); put(low); new_line;
      put("  upper : "); put(upp); new_line;
      put("Do you want to change these values ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      put("Reading "); put(m,1); put(" lower bounds ");
      put_line("for the lifting."); Read_Integer_Vector(low);
      put("Reading "); put(m,1); put(" upper bounds ");
      put_line("for the lifting."); Read_Integer_Vector(upp);
    end loop;
    put(file,"  Lower bounds : "); put(file,low); new_line(file);
    put(file,"  Upper bounds : "); put(file,upp); new_line(file);
  end Set_Integer_Bounds;

  procedure Set_Float_Bounds
              ( file : in file_type;
                low,upp : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows the default values for lower and upper bounds and allows
  --   the user to modify these.

    ans : character;
    m : constant integer32 := low'length;

  begin
    loop
      new_line;
      put_line("Current lower and upper bounds on lifting values");
      put("  lower : "); put(low,2,3,3); new_line;
      put("  upper : "); put(upp,2,3,3); new_line;
      put("Do you want to change these values ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      put("Reading "); put(m,1); put(" lower bounds ");
      put_line("for the lifting."); Read_Float_Vector(low);
      put("Reading "); put(m,1); put(" upper bounds ");
      put_line("for the lifting."); Read_Float_Vector(upp);
    end loop;
    put(file,"  Lower bounds : "); put(file,low,2,3,3); new_line(file);
    put(file,"  Upper bounds : "); put(file,upp,2,3,3); new_line(file);
  end Set_Float_Bounds;

  function Read_Integer_Lifting
              ( L : Lists_of_Integer_Vectors.List )
              return Lists_of_Integer_Vectors.List is

    use Lists_of_Integer_Vectors;
    use Standard_Integer_Vectors;

    tmp : List := L;
    res,res_last : List;

  begin
    put_line("Give a lifting value for the following points :");
    while not Is_Null(tmp) loop
      declare
        v : constant Vector := Head_Of(tmp).all;
        lv : Vector(v'first..v'last+1);
      begin
        lv(v'range) := v;
        put(v); put(" : "); Read_Integer(lv(lv'last));
        Append(res,res_last,lv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Read_Integer_Lifting;

  function Read_Float_Lifting
              ( L : Lists_of_Floating_Vectors.List ) 
              return Lists_of_Floating_Vectors.List is

    use Lists_of_Floating_Vectors;
    use Standard_Floating_Vectors;

    tmp : List := L;
    res,res_last : List;

  begin
    put_line("Give a lifting value for the following points :");
    while not Is_Null(tmp) loop
      declare
        v : constant Vector := Head_Of(tmp).all;
        lv : Vector(v'first..v'last+1);
      begin
        lv(v'range) := v;
        put(v,0,0,0); put(" : "); Read_Double_Float(lv(lv'last));
        Append(res,res_last,lv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Read_Float_Lifting;

  procedure Integer_Random_Linear_Lifting
             ( file : in file_type; n : in integer32;
               points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lilifu : in out Standard_Integer_VecVecs.Link_to_VecVec ) is

    low : Standard_Integer_Vectors.Vector(points'range) := (points'range => 0);
    upp : Standard_Integer_Vectors.Vector(points'range)
        := Adaptive_Lifting(points);
    ranvec : Standard_Integer_Vectors.Vector(1..n);

  begin
    Set_Integer_Bounds(file,low,upp);
    lilifu := new Standard_Integer_VecVecs.VecVec(points'range);
    for i in lilifu'range loop
      for j in ranvec'range loop
        ranvec(j) := Random(low(i),upp(i));
      end loop;
      lilifu(i) := new Standard_Integer_Vectors.Vector'(ranvec);
    end loop;
    lifted := Linear_Lift(lilifu.all,points);
  end Integer_Random_Linear_Lifting;

  procedure Float_Random_Linear_Lifting
             ( file : in file_type; n : in integer32;
               points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               lilifu : in out Standard_Floating_VecVecs.Link_to_VecVec ) is

    m : constant integer32 := points'last;
    low : Standard_Floating_Vectors.Vector(points'range)
        := (points'range => 0.0);
    upp : Standard_Floating_Vectors.Vector(points'range)
        := Adaptive_Lifting(points);
    ranvec : Standard_Floating_Vectors.Vector(1..n);

  begin
    Set_Float_Bounds(file,low,upp);
    lilifu := new Standard_Floating_VecVecs.VecVec(points'range);
    for i in lilifu'range loop
      ranvec := Random(m,low(i),upp(i));
      lifted(i) := Linear_Lift(points(i),ranvec);
      lilifu(i) := new Standard_Floating_Vectors.Vector'(ranvec);
    end loop;
  end Float_Random_Linear_Lifting;

  procedure Integer_User_Linear_Lifting
             ( file : in file_type; n : in integer32;
               points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lilifu : in out Standard_Integer_VecVecs.Link_to_VecVec ) is

    lift : Standard_Integer_VecVecs.VecVec(points'range);

  begin
    for k in lift'range loop
      put("Give "); put(n,1); put(" integers for ");
      put("vector "); put(k,1); put(" : "); get(natural32(n),lift(k));
      put(file," vector "); put(file,k,1);
      put(file," : "); put(file,lift(k)); new_line(file);
    end loop;
    lifted := Linear_Lift(lift,points);
    lilifu := new Standard_Integer_VecVecs.VecVec'(lift);
  end Integer_User_Linear_Lifting;

  procedure Float_User_Linear_Lifting
             ( file : in file_type; n : in integer32;
               points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               lilifu : in out Standard_Floating_VecVecs.Link_to_VecVec ) is

    lift : Standard_Floating_VecVecs.VecVec(points'range);

  begin
    for k in lift'range loop
      put("Give "); put(n,1); put(" floats for ");
      put("vector "); put(k,1); put(" : "); get(natural32(n),lift(k));
      put(file," vector "); put(file,k,1);
      put(file," : "); put(file,lift(k)); new_line(file);
      lifted(k) := Linear_Lift(points(k),lift(k).all);
    end loop;
    lilifu := new Standard_Floating_VecVecs.VecVec'(lift);
  end Float_User_Linear_Lifting;

  procedure Integer_Polynomial_Lifting 
            ( file : in file_type;
              points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lftd : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is

    lift : Link_to_Poly_Sys;

  begin
    get(lift);
    put(file,lift.all);
    lftd := Polynomial_Lift(lift.all,points);
  end Integer_Polynomial_Lifting;

  procedure Float_Polynomial_Lifting
          ( file : in file_type;
            points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
            lftd : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists ) is

    lift : Link_to_Poly_Sys;

  begin
    get(lift);
    put(file,lift.all);
    lftd := Polynomial_Lift(lift.all,points);
  end Float_Polynomial_Lifting;

  procedure Integer_User_Point_Wise_Lifting
            ( points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lftd : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is
  begin
    for k in points'range loop
      lftd(k) := Read_Integer_Lifting(points(k));
    end loop;
  end Integer_User_Point_Wise_Lifting;

  procedure Float_User_Point_Wise_Lifting
          ( points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
            lftd : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists ) is
  begin
    for k in points'range loop
      lftd(k) := Read_Float_Lifting(points(k));
    end loop;
  end Float_User_Point_Wise_Lifting;

  procedure Integer_Random_Point_Wise_Lifting
            ( file : in file_type;
              points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is

    low : Standard_Integer_Vectors.Vector(points'range) := (points'range => 0);
    upp : Standard_Integer_Vectors.Vector(points'range)
        := Adaptive_Lifting(points);

  begin
    Set_Integer_Bounds(file,low,upp);
    lifted := Random_Lift(low,upp,points);
  end Integer_Random_Point_Wise_Lifting;

  procedure Float_Random_Point_Wise_Lifting
          ( file : in file_type;
            points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
            lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists ) is

    low : Standard_Floating_Vectors.Vector(points'range)
        := (points'range => 0.0);
    upp : Standard_Floating_Vectors.Vector(points'range)
        := Adaptive_Lifting(points); -- := (points'range => 1.0); 

  begin
    Set_Float_Bounds(file,low,upp);
    lifted := Random_Lift(points,low,upp);
  end Float_Random_Point_Wise_Lifting;

  procedure Float_Stable_Lifting
          ( file : in file_type; b : in double_float;
            points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
            lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists ) is

  begin
    put(file,"  Lifting bound for stable mixed volume : "); 
    put(file,b); new_line(file);
    lifted := Stable_Lift(points,b);
  end Float_Stable_Lifting;

  function Menu_for_Lifting_Functions ( gl : in boolean ) return character is

  -- DESCRIPTION :
  --   Displays the menu and returns the selected lifting function.

  -- ON ENTRY :
  --   gl       true when the system is "genuine Laurent" and computing
  --            stable mixed volumes is not an option.

    ans : character;

  begin
    put_line("MENU for Lifting Functions :");
    put_line("  Point-wise : 0. For each point a random lifting value.");
    put_line("             : 1. You can choose all lifting values yourself.");
    put_line("  Linear     : 2. Random linear vectors will be chosen.");
    put_line("             : 3. You can choose the vectors yourself.");
    put_line("  Polynomial : 4. The system to be solved will be used.");
    put_line("             : 5. You can choose the polynomials yourself.");
    if not gl then
      put_line("  6. One single lifting to compute stable mixed volumes.");
      put("Type a number between 0 and 6 to select lifting : ");
      Ask_Alternative(ans,"0123456"); return ans;
    else
      put("Type a number between 0 and 5 to select lifting : ");
      Ask_Alternative(ans,"012345"); return ans;
    end if;
  end Menu_for_Lifting_Functions;

  procedure Dispatch_Integer_Lifting
              ( file : in file_type; p : in Poly_Sys; choice : in character;
                points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lilifu : in out Standard_Integer_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Calls the appropriate lifting function to perform the lifting.

    n : constant integer32 := p'length;

  begin
    new_line(file); put(file,"INTEGER LIFTING FUNCTION : ");
    case choice is
      when '0' => put_line(file,"random point-wise.");
                  Integer_Random_Point_Wise_Lifting(file,points,lifted);
      when '1' => put_line(file,"point-wise provided by user.");
                  Integer_User_Point_Wise_Lifting(points,lifted);
      when '2' => put_line(file,"random linear.");
                  Integer_Random_Linear_Lifting(file,n,points,lifted,lilifu);
      when '3' => put_line(file,"linear provided by user.");
                  Integer_User_Linear_Lifting(file,n,points,lifted,lilifu);
      when '4' => put_line(file,"polynomial system.");
                  lifted := Polynomial_Lift(p,points);
      when '5' => put_line(file,"polynomials provided by user.");
                  Integer_Polynomial_Lifting(file,points,lifted);
      when others => null;
    end case;
    new_line(file); put_line(file,"THE LIFTED SUPPORTS :"); new_line(file);
    put(file,lifted);
  end Dispatch_Integer_Lifting;

  procedure Dispatch_Integer_Lifting
              ( file : in file_type; p : in Laur_Sys; choice : in character;
                points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lilifu : in out Standard_Integer_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Calls the appropriate lifting function to perform the lifting.
  --   This version allows Laurent polynomial systems on input.

    n : constant integer32 := p'length;

  begin
    new_line(file); put(file,"INTEGER LIFTING FUNCTION : ");
    case choice is
      when '0' => put_line(file,"random point-wise.");
                  Integer_Random_Point_Wise_Lifting(file,points,lifted);
      when '1' => put_line(file,"point-wise provided by user.");
                  Integer_User_Point_Wise_Lifting(points,lifted);
      when '2' => put_line(file,"random linear.");
                  Integer_Random_Linear_Lifting(file,n,points,lifted,lilifu);
      when '3' => put_line(file,"linear provided by user.");
                  Integer_User_Linear_Lifting(file,n,points,lifted,lilifu);
      when '4' => put_line(file,"polynomial system.");
                  lifted := Polynomial_Lift(p,points);
      when '5' => put_line(file,"polynomials provided by user.");
                  Integer_Polynomial_Lifting(file,points,lifted);
      when others => null;
    end case;
    new_line(file); put_line(file,"THE LIFTED SUPPORTS :"); new_line(file);
    put(file,lifted);
  end Dispatch_Integer_Lifting;

  procedure Dispatch_Float_Lifting
              ( file : in file_type; p : in Poly_Sys; choice : in character;
                b : in double_float;
                points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lilifu : in out Standard_Floating_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Calls the appropriate lifting function to perform the lifting.

    n : constant integer32 := p'length;

  begin
    new_line(file); put(file,"FLOATING-POINT LIFTING FUNCTION : ");
    case choice is
      when '0' => put_line(file,"random point-wise.");
                  Float_Random_Point_Wise_Lifting(file,points,lifted);
      when '1' => put_line(file,"point-wise provided by user.");
                  Float_User_Point_Wise_Lifting(points,lifted);
      when '2' => put_line(file,"random linear.");
                  Float_Random_Linear_Lifting(file,n,points,lifted,lilifu);
      when '3' => put_line(file,"linear provided by user.");
                  Float_User_Linear_Lifting(file,n,points,lifted,lilifu);
      when '4' => put_line(file,"polynomial system.");
                  lifted := Polynomial_Lift(p,points);
      when '5' => put_line(file,"polynomials provided by user.");
                  Float_Polynomial_Lifting(file,points,lifted);
      when '6' => put_line(file,"lifting for stable mixed volume.");
                  Float_Stable_Lifting(file,b,points,lifted);               
      when others => null;
    end case;
    new_line(file); put_line(file,"THE LIFTED SUPPORTS :"); new_line(file);
    Floating_Mixed_Subdivisions_io.put(file,lifted);
  end Dispatch_Float_Lifting;

  procedure Dispatch_Float_Lifting
              ( file : in file_type; p : in Laur_Sys; choice : in character;
                b : in double_float;
                points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lilifu : in out Standard_Floating_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Calls the appropriate lifting function to perform the lifting.
  --   This version allows Laurent polynomial systems on input.

    n : constant integer32 := p'length;

  begin
    new_line(file); put(file,"FLOATING-POINT LIFTING FUNCTION : ");
    case choice is
      when '0' => put_line(file,"random point-wise.");
                  Float_Random_Point_Wise_Lifting(file,points,lifted);
      when '1' => put_line(file,"point-wise provided by user.");
                  Float_User_Point_Wise_Lifting(points,lifted);
      when '2' => put_line(file,"random linear.");
                  Float_Random_Linear_Lifting(file,n,points,lifted,lilifu);
      when '3' => put_line(file,"linear provided by user.");
                  Float_User_Linear_Lifting(file,n,points,lifted,lilifu);
      when '4' => put_line(file,"polynomial system.");
                  lifted := Polynomial_Lift(p,points);
      when '5' => put_line(file,"polynomials provided by user.");
                  Float_Polynomial_Lifting(file,points,lifted);
      when '6' => put_line(file,"lifting for stable mixed volume.");
                  Float_Stable_Lifting(file,b,points,lifted);
      when others => null;
    end case;
    new_line(file); put_line(file,"THE LIFTED SUPPORTS :"); new_line(file);
    Floating_Mixed_Subdivisions_io.put(file,lifted);
  end Dispatch_Float_Lifting;

-- TARGET ROUTINES :

  procedure Integer_Lifting
            ( file : in file_type; p : in Poly_Sys;
              points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lilifu : in out Standard_Integer_VecVecs.Link_to_VecVec ) is

    liftfun : constant character := Menu_for_Lifting_Functions(false);

  begin
    Dispatch_Integer_Lifting(file,p,liftfun,points,lifted,lilifu);
  end Integer_Lifting;

  procedure Floating_Lifting
            ( file : in file_type; p : in Poly_Sys; b : in double_float;
              points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
              lilifu : in out Standard_Floating_VecVecs.Link_to_VecVec ) is

    liftfun : constant character := Menu_for_Lifting_Functions(false);

  begin
    Dispatch_Float_Lifting(file,p,liftfun,b,points,lifted,lilifu);
  end Floating_Lifting;

  procedure Main_Polynomial
              ( file : in file_type; p : in Poly_Sys;
                ipoints : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                fltlif : out boolean; stlb : out double_float;
                fpts : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ilftd : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flftd : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ililifu : in out Standard_Integer_VecVecs.Link_to_VecVec;
                flilifu : in out Standard_Floating_VecVecs.Link_to_VecVec ) is

    liftfun : constant character := Menu_for_Lifting_Functions(false);
    ans : character;

  begin
    if liftfun /= '6' then
      stlb := 0.0;
      put("Change default floating-point lifting to integer? (y/n) ");
      Ask_Yes_or_No(ans);
    else
      stlb := Lifting_Bound(p);
      ans := 'n'; -- only floating-point lifting for stable mixed volumes
    end if;
   -- put("Type i for integer or f for floating-point lifting : ");
   -- Ask_Alternative(ans,"if"); fltlif := (ans = 'f');
   -- if ans = 'i'
    if ans = 'y' then
      fltlif := false;
      Dispatch_Integer_Lifting(file,p,liftfun,ipoints,ilftd,ililifu);
    else
      fltlif := true;
      fpts := Floating_Integer_Convertors.Convert(ipoints);
      Dispatch_Float_Lifting(file,p,liftfun,stlb,fpts,flftd,flilifu);
    end if;
  end Main_Polynomial;

  procedure Main_Laurent
              ( file : in file_type; p : in Laur_Sys;
                ipoints : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                fltlif : out boolean; stlb : out double_float;
                fpts : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ilftd : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flftd : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ililifu : in out Standard_Integer_VecVecs.Link_to_VecVec;
                flilifu : in out Standard_Floating_VecVecs.Link_to_VecVec ) is

    gl : constant boolean := Is_Genuine_Laurent(p);
    liftfun : constant character := Menu_for_Lifting_Functions(gl);
    ans : character;

  begin
    if liftfun /= '6' then
      stlb := 0.0;
      put("Change default floating-point lifting to integer? (y/n) ");
      Ask_Yes_or_No(ans);
    else
      stlb := Lifting_Bound(p);
      ans := 'n'; -- only floating-point lifting for stable mixed volumes
    end if;
   -- put("Type i for integer or f for floating-point lifting : ");
   -- Ask_Alternative(ans,"if"); fltlif := (ans = 'f');
   -- if ans = 'i'
    if ans = 'y' then
      fltlif := false;
      Dispatch_Integer_Lifting(file,p,liftfun,ipoints,ilftd,ililifu);
    else
      fltlif := true;
      fpts := Floating_Integer_Convertors.Convert(ipoints);
      Dispatch_Float_Lifting(file,p,liftfun,stlb,fpts,flftd,flilifu);
    end if;
  end Main_Laurent;

end Main_Lifting_Functions;
