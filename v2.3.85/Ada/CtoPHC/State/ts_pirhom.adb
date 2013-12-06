with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Pieri_Homotopy;

procedure ts_pirhom is

  procedure Store_Random_Input_Planes ( m,p,q : in natural32 ) is

    n : constant integer32 := integer32(m*p + q*(m+p));
    planes : VecMat(1..n);

  begin
    for i in 1..n loop
      planes(i) := new Matrix'(Random_Orthogonal_Matrix(m+p,m));
    end loop;
    Pieri_Homotopy.Initialize_Input_Planes(planes);
  end Store_Random_Input_Planes;

  procedure Store_Random_Interpolation_Points ( m,p,q : in natural32 ) is

    n : constant integer32 := integer32(m*p + q*(m+p));
    points : Standard_Complex_Vectors.Vector(1..n)
           := Random_Vector(1,n);

  begin
    Pieri_Homotopy.Initialize_Interpolation_Points(points);
  end Store_Random_Interpolation_Points;

  function Start_Bracket
             ( p : integer32 ) return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(1..p);

  begin
    for i in 1..p loop
      res(i) := natural32(i);
    end loop;
    return res;
  end Start_Bracket;

  procedure Write_Bracket
              ( b : in Standard_Natural_Vectors.Vector ) is

  begin
    put("["); put(b); put(" ]");
  end Write_Bracket;

  procedure Write_Bracket_Pair
              ( b1,b2 : in Standard_Natural_Vectors.Vector ) is

  begin
    put("(");
    Write_Bracket(b1);
    put(",");
    Write_Bracket(b2);
    put(")");
  end Write_Bracket_Pair;

  procedure Call_Pieri_Homotopy ( m,p,q : in natural32 ) is

    nb : constant integer32 := integer32(m*p + q*(m+p));
    start_top,start_bottom
      : Standard_Natural_Vectors.Vector(1..integer32(p));
    target_top,target_bottom
      : Standard_Natural_Vectors.Vector(1..integer32(p));
    n : natural32;
    x : Standard_Complex_Vectors.Vector(1..nb);
    ans,output : character;
    res : double_float;
    
  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(output);
    Pieri_Homotopy.Initialize_Dimensions(m,p,q);
    Store_Random_Input_Planes(m,p,q);
    Store_Random_Interpolation_Points(m,p,q);
    new_line;
    start_top := Start_Bracket(integer32(p));
    start_bottom := Start_Bracket(integer32(p));
    put("top and bottom pivots at start : ");
    Write_Bracket_Pair(start_top,start_bottom); new_line;
    target_top := Start_Bracket(integer32(p));
    loop
      new_line;
     -- put("Give "); put(p,1); put(" top pivots at target : ");  
     -- get(target_top);
      put("Give "); put(p,1); put(" bottom pivots at target : ");  
      get(target_bottom);
      Pieri_Homotopy.Store_Start_Pivots(start_top,start_bottom);
      n := Pieri_Homotopy.Degree_of_Freedom(start_top,start_bottom);
      if n = 0
       then x := (1..nb => Create(0.0));
      end if;
      Pieri_Homotopy.Store_Start(x(1..integer32(n)));
      Pieri_Homotopy.Store_Target_Pivots(target_top,target_bottom);
      if output = 'y'
       then Pieri_Homotopy.Track_Path(Standard_Output);
       else Pieri_Homotopy.Track_Path;
      end if;
      n := Pieri_Homotopy.Degree_of_Freedom(target_top,target_bottom);
      Pieri_Homotopy.Retrieve_Target(x(1..integer32(n)));
      put_line("The solution : "); put_line(x(1..integer32(n)));
      if output = 'y'
       then res := Pieri_Homotopy.Verify_Determinants(Standard_Output);
       else res := Pieri_Homotopy.Verify_Determinants;
      end if;
      put("Residual at intersection : "); put(res,3); new_line;
      put("Continue to the next step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      start_top := target_top; start_bottom := target_bottom;
      put("Current pivots at start : ");
      Write_Bracket_Pair(start_top,start_bottom); new_line;
    end loop;
    Pieri_Homotopy.Clear;
  end Call_Pieri_Homotopy;

  procedure Main is

    m,p,q : natural32 := 0;

  begin
    new_line;
    put_line("Welcome to the Pieri Homotopy machine...");
    new_line;
    put("Give m (dimension of input planes)  : "); get(m);
    put("Give p (dimension of output planes) : "); get(p);
    put("Give q (degree of solution curves)  : "); get(q);
    Call_Pieri_Homotopy(m,p,q);
  end Main;

begin
  Main;
end ts_pirhom;
