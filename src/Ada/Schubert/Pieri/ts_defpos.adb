with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Symbol_Table;                       use Symbol_Table;
with Matrix_Indeterminates;
with Main_Poly_Continuation;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Localization_Posets;                use Localization_Posets;
with Localization_Posets_io;             use Localization_Posets_io;
with Curves_into_Grassmannian_io;        use Curves_into_Grassmannian_io;
with Deformation_Posets;                 use Deformation_Posets;
with Make_Input_Planes;

procedure ts_defpos is

-- DESCRIPTION :
--   Test on the deformation posets.

  procedure Add_t_Symbol is

  -- DESCRIPTION :
  --   Adds the symbol for the continuation parameter t to the symbol table.

    tsb : Symbol;

  begin
    Symbol_Table.Enlarge(1);
    tsb(1) := 't';
    for i in 2..tsb'last loop
      tsb(i) := ' ';
    end loop;
    Symbol_Table.Add(tsb);
  end Add_t_Symbol;

  procedure Set_Parameters ( file : in file_type; report : out boolean ) is

  -- DESCRIPTION :
  --   Interactive determination of the continuation and output parameters.

    oc : natural32;

    use Main_Poly_Continuation;

  begin
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := not (oc = 0);
    new_line;
    put_line("No more input expected.  See output file for results...");
    new_line;
    new_line(file);
  end Set_Parameters;

  function Random_Input_Planes
             ( m,p : natural32; k : Bracket ) return VecMat is

  -- DESCRIPTION :
  --   Returns a vector of random (m+1-k(i))-planes, for i in k'range.

    res : VecMat(k'range);
    dim : constant natural32 := m+p;

  begin
    for i in res'range loop
      res(i)
        := new Standard_Complex_Matrices.Matrix'(Random_Matrix(dim,m+1-k(i)));
    end loop;
    return res;
  end Random_Input_Planes;

  procedure Solve_Deformation_Poset
               ( file : in file_type; m,p,q : in natural32;
                 level_poset : in Array_of_Nodes;
                 index_poset : in Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates a deformation poset and applies the Solve operator.

    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    mpq : constant integer32 := integer32(m*p + q*(m+p));
    planes : constant VecMat(1..mpq)
           := Make_Input_Planes.Random_Complex_Planes(m,p,q);
    svals : constant Standard_Complex_Vectors.Vector := Random_Vector(1,mpq);
    ans : character;
    report,outlog : boolean;
    timer : Timing_Widget;
    target_level : natural32 := natural32(mpq);
    nbp : natural32 := 0;
    npaths : Standard_Natural_Vectors.Vector(1..integer32(target_level))
           := (1..integer32(target_level) => 0);
    timings : Duration_Array(1..integer32(target_level))
            := (1..integer32(target_level) => 0.0);
    nb : natural32 := 0;

  begin
    put_line("The size of the deformation poset : ");
    put_line(file,"The size of the deformation poset : ");
    put_roco(index_poset);
    put_roco(file,index_poset);
    new_line;
    put("Give target level <= "); put(target_level,1);
    put(" = root level : "); get(target_level);
    for i in 1..target_level loop
      nbp := nbp + Row_Root_Count_Sum(level_poset,i);
    end loop;
    put("The number of paths : "); put(nbp,1); new_line;
    put(file,"The number of paths : "); put(file,nbp,1); new_line(file);
    if q > 0 then
      new_line;
      put("How many solution maps do you want ? "); get(nb);
    end if;
    new_line;
    put("Do you want to have the homotopies on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    outlog := (ans = 'y');
    if q = 0 then
      Matrix_Indeterminates.Initialize_Symbols(m+p,p);
      Add_t_Symbol;
    else
      One_Set_up_Symbol_Table
        (m,p,q,index_poset(mpq)(1).top,index_poset(mpq)(1).bottom);
    end if;
    Set_Parameters(file,report);
    tstart(timer);
    for i in index_poset(integer32(target_level))'range loop
      declare
        root : constant Node := index_poset(integer32(target_level))(i).all;
      begin
        if q = 0 then
          Solve(file,m+p,deform_poset,root,planes,report,outlog,
                npaths,timings);
        else
          Solve(file,m+p,q,nb,deform_poset,root,planes,svals,
                report,outlog,npaths,timings);
        end if;     
      end;
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Solving along the deformation poset");
  end Solve_Deformation_Poset;

  procedure Solve_Deformation_Poset
               ( file : in file_type; m,p,q : in natural32; k : in Bracket;
                 index_poset : in Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Applies the solver to general intersection conditions.

    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    mpq : constant natural32 := m*p + q*(m+p);
    planes : constant VecMat(k'range) := Random_Input_Planes(m,p,k);
    svals : constant Standard_Complex_Vectors.Vector
          := Random_Vector(1,k'last);
    ans : character;
    report,outlog : boolean;
    timer : Timing_Widget;
    target_level : natural32 := mpq;
    npaths : Standard_Natural_Vectors.Vector(1..integer32(target_level))
           := (1..integer32(target_level) => 0);
    timings : Duration_Array(1..integer32(target_level))
            := (1..integer32(target_level) => 0.0);
    nb : natural32 := 0;

  begin
    put_line("The size of the deformation poset : ");
    put_line(file,"The size of the deformation poset : ");
    put_roco(index_poset);
    put_roco(file,index_poset);
    new_line;
    put("Give target level <= "); put(target_level,1);
    put(" = root level : "); get(target_level);
    if q > 0 then
      new_line;
      put("How many solution maps do you want ? "); get(nb);
    end if;
    new_line;
    put("Do you want to have the homotopies on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    outlog := (ans = 'y');
    Matrix_Indeterminates.Initialize_Symbols(m+p,p);
    Add_t_Symbol;
    Set_Parameters(file,report);
    tstart(timer);
    for i in index_poset(integer32(target_level))'range loop
      declare
        root : constant Node := index_poset(integer32(target_level))(i).all;
        cnt : natural32 := 0;
      begin
        if q = 0 then
          Solve(file,m+p,k,deform_poset,root,planes,report,outlog,
                npaths,timings);
        else
          One_Solve(file,m+p,q,nb,cnt,k,deform_poset,root,planes,svals,
                    report,outlog,npaths,timings);
        end if;
      end;
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Solving along the deformation poset");
  end Solve_Deformation_Poset;

  procedure Create_Top_Hypersurface_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing only top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    Top_Create(lnkroot,m+p);
    put_line(file,"  q = 0");
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,0,level_poset,index_poset);
  end Create_Top_Hypersurface_Poset;

  procedure Create_Top_Hypersurface_Poset
              ( file : in file_type; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing only top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : Link_to_Node := new Node'(root);
    mpq : constant integer32 := integer32(m*p + q*(m+p));
    level_poset : Array_of_Nodes(0..mpq);
    index_poset : Array_of_Array_of_Nodes(0..mpq);

  begin
    Q_Top_Create(lnkroot,root.bottom(integer32(p)),m+p);
    put(file,"  q = "); put(file,q,1); new_line(file);
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,q,level_poset,index_poset);
  end Create_Top_Hypersurface_Poset;

  procedure Create_Bottom_Hypersurface_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by decrementing only bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    Bottom_Create(lnkroot);
    put_line(file,"  q = 0");
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,0,level_poset,index_poset);
  end Create_Bottom_Hypersurface_Poset;

  procedure Create_Bottom_Hypersurface_Poset
              ( file : in file_type; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by decrementing only bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    mpq : constant integer32 := integer32(m*p + q*(m+p));
    level_poset : Array_of_Nodes(0..mpq);
    index_poset : Array_of_Array_of_Nodes(0..mpq);

  begin
    Q_Bottom_Create(lnkroot,m+p);
    put(file,"  q = "); put(file,q,1); new_line(file);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,q,level_poset,index_poset);
  end Create_Bottom_Hypersurface_Poset;

  procedure Create_Mixed_Hypersurface_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top and decrementing bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    Top_Bottom_Create(lnkroot,m+p);
    put_line(file,"  q = 0");
    put_line("The poset created in a mixed fashion : ");
    put_line(file,"The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,0,level_poset,index_poset);
  end Create_Mixed_Hypersurface_Poset;

  procedure Create_Mixed_Hypersurface_Poset
              ( file : in file_type; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top and decrementing bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    mpq : constant integer32 := integer32(m*p + q*(m+p));
    level_poset : Array_of_Nodes(0..mpq);
    index_poset : Array_of_Array_of_Nodes(0..mpq);

  begin
    Q_Top_Bottom_Create(lnkroot,root.bottom(integer32(p)),m+p);
    put(file,"  q = "); put(file,q,1); new_line(file);
    put_line("The poset created in a mixed fashion : ");
    put_line(file,"The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,q,level_poset,index_poset);
  end Create_Mixed_Hypersurface_Poset;

  --function Finite ( dim : Bracket; m,p : natural32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the codimensions yield a finite number of solutions.

  --  sum : natural32 := 0;

  --begin
  --  for i in dim'range loop
  --    sum := sum + dim(i);
  --  end loop;
  --  if sum = m*p
  --   then return true;
  --   else return false;
  --  end if;
  --end Finite;

  procedure Create_Top_General_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    put(file,"  q = 0");
    put(file,"  k = "); put(file,codim); new_line(file);
    Top_Create(lnkroot,codim,m+p);
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);                           
    Solve_Deformation_Poset(file,m,p,0,codim,index_poset);
  end Create_Top_General_Poset;

  procedure Create_Top_General_Poset
              ( file : in file_type; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);
    mpq : constant integer32 := integer32(m*p + q*(m+p));
    level_poset : Array_of_Nodes(0..mpq);
    index_poset : Array_of_Array_of_Nodes(0..mpq);

  begin
    put(file,"  q = "); put(file,q,1);
    put(file,"  k = "); put(file,codim); new_line(file);
    Q_Top_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);                           
    Solve_Deformation_Poset(file,m,p,q,codim,index_poset);
  end Create_Top_General_Poset;

  procedure Create_Bottom_General_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    put(file,"  q = 0");
    put(file,"  k = "); put(file,codim); new_line(file);
    Bottom_Create(lnkroot,codim);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,0,codim,index_poset);
  end Create_Bottom_General_Poset;

  procedure Create_Bottom_General_Poset
              ( file : in file_type; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);
    mpq : constant integer32 := integer32(m*p + q*(m+p));
    level_poset : Array_of_Nodes(0..mpq);
    index_poset : Array_of_Array_of_Nodes(0..mpq);

  begin
    put(file,"  q = "); put(file,q,1);
    put(file,"  k = "); put(file,codim); new_line(file);
    Q_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,q,codim,index_poset);
  end Create_Bottom_General_Poset;

  procedure Create_Mixed_General_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by incrementing the top and decrementing the bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    put(file,"  q = 0");
    put(file,"  k = "); put(file,codim); new_line(file);
    Top_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created in a mixed fashion : ");
    put_line(file,"The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,0,codim,index_poset);
  end Create_Mixed_General_Poset;

  procedure Create_Mixed_General_Poset
              ( file : in file_type; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by incrementing the top and decrementing the bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);
    mpq : constant integer32 := integer32(m*p + q*(m+p));
    level_poset : Array_of_Nodes(0..mpq);
    index_poset : Array_of_Array_of_Nodes(0..mpq);

  begin
    put(file,"  q = "); put(file,q,1);
    put(file,"  k = "); put(file,codim); new_line(file);
    Q_Top_Bottom_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created in a mixed fashion : ");
    put_line(file,"The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Deformation_Poset(file,m,p,q,codim,index_poset);
  end Create_Mixed_General_Poset;

  procedure Main is

    m,p,q : natural32 := 0;
    ans : character;
    file : file_type;

  begin
    new_line;
    put_line("MENU for posets for deforming p-planes in (m+p)-space : ");
    put_line("------ the case q = 0 ----------------------------------------");
    put_line("  1. k_i == 1 consistently incrementing the top pivots.");
    put_line("  2.          consistently decrementing the bottom pivots.");
    put_line("  3.          mixed top-bottom sequence for poset creation.");
    put_line("  4. k_i >= 1 consistently incrementing the top pivots.");
    put_line("  5.          consistently decrementing the bottom pivots.");
    put_line("  6.          mixed top-bottom sequence for poset creation.");
    put_line("------ the case q >= 0 ---------------------------------------");
    put_line("  7. k_i == 1 consistently incrementing the top pivots.");
    put_line("  8.          consistently decrementing the bottom pivots.");
    put_line("  9.          mixed top-bottom sequence for poset creation.");
    put_line("  A. k_i >= 1 consistently incrementing the top pivots.");
    put_line("  B.          consistently decrementing the bottom pivots.");
    put_line("  C.          mixed top-bottom sequence for poset creation.");
    put_line("--------------------------------------------------------------");
    put("Type 1,2,3,4,5,6,7,9,A,B,C to choose : ");
    Ask_Alternative(ans,"123456789ABC");
    new_line;
    put_line("Reading the name of the file for the deformations.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give p, the number of entries in bracket : "); get(p);
    put("Give m, the complementary dimension : "); get(m);
    put(file,"p = "); put(file,p,1); put(file,"  m = "); put(file,m,1);
    new_line;
    case ans is
      when '1' => Create_Top_Hypersurface_Poset(file,m,p);
      when '2' => Create_Bottom_Hypersurface_Poset(file,m,p);
      when '3' => Create_Mixed_Hypersurface_Poset(file,m,p);
      when '4' => Create_Top_General_Poset(file,m,p);
      when '5' => Create_Bottom_General_Poset(file,m,p);
      when '6' => Create_Mixed_General_Poset(file,m,p);
      when '7' => put("Give q, the degree of the maps : "); get(q);
                  Create_Top_Hypersurface_Poset(file,m,p,q);
      when '8' => put("Give q, the degree of the maps : "); get(q);
                  Create_Bottom_Hypersurface_Poset(file,m,p,q);
      when '9' => put("Give q, the degree of the maps : "); get(q);
                  Create_Mixed_Hypersurface_Poset(file,m,p,q);
      when 'A' => put("Give q, the degree of the maps : "); get(q);
                  Create_Top_General_Poset(file,m,p,q);
      when 'B' => put("Give q, the degree of the maps : "); get(q);
                  Create_Bottom_General_Poset(file,m,p,q);
      when 'C' => put("Give q, the degree of the maps : "); get(q);
                  Create_Mixed_General_Poset(file,m,p,q);
      when others => put_line("Option not recognized.  Please try again.");
    end case;
  end Main;

begin
  new_line;
  put_line("Test on deformation posets for linear subspace intersections.");
  Main;
end ts_defpos;
