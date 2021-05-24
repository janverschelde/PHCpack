with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Faces_of_Polytope;
with Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Floating_Lifting_Utilities;         use Floating_Lifting_Utilities;
with Pruning_Statistics;
with Integer_Pruning_Methods;            use Integer_Pruning_Methods;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Floating_Integer_Convertors;        use Floating_Integer_Convertors;
with Floating_Faces_of_Polytope;
with Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Floating_Pruning_Methods;           use Floating_Pruning_Methods;
with Symmetry_Group;                     use Symmetry_Group;
with Symmetry_Group_io;                  use Symmetry_Group_io;
with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Symbolic_Symmetry_Group_io;         use Symbolic_Symmetry_Group_io;
with Drivers_for_Symmetry_Group_io;      use Drivers_for_Symmetry_Group_io;
with Equivariant_Polynomial_Systems;     use Equivariant_Polynomial_Systems;
with Faces_of_Symmetric_Polytopes;       use Faces_of_Symmetric_Polytopes;
with Symmetric_Lifting_Functions;        use Symmetric_Lifting_Functions;
with Generating_Mixed_Cells;             use Generating_Mixed_Cells;
with Symmetric_Randomize;
with Symmetric_Polyhedral_Continuation;  use Symmetric_Polyhedral_Continuation;

package body Drivers_for_Symmetric_Lifting is

  procedure Symmetric_Lifting_Info is

    i : array(1..6) of string(1..65);

  begin
    i(1):="  Symmetric lifting allows one to exploit permutation  symmetries  in";
    i(2):="the  tuple  of  Newton  polytopes.   A  symmetric  subdivision is";
    i(3):="induced by lifting the points in the same orbit up  to  the  same";
    i(4):="height.   The  corresponding  random coefficient start system has";
    i(5):="the same symmetric structure, so that in the homotopy,  only  the";
    i(6):="generating solution paths need to be traced.                     ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Symmetric_Lifting_Info;

  procedure Driver_for_Symmetric_Mixed_Volume_Computation
              ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                q : out Poly_Sys; qsols : out Solution_List;
                mv : out natural32 ) is

    welcome : constant string
            := "Mixed-Volume Computation by Symmetric Lifting";

    solsft,gft,outsubft : file_type;

   -- SWITCHES :

    invariant : boolean;   -- true if the polynomials are invariant
    equivaria : boolean;   -- true if the system is equivariant
    compmisu  : boolean;   -- if a mixed subdivision has to be computed
    misufile  : boolean;   -- when the mixed subdivision has to be put on file
    signsym   : boolean;   -- there is sign symmetry 
    allperms  : boolean;   -- equi-invariant w.r.t. all permutations
    allsigns  : boolean;   -- equi-invariant w.r.t. all sign permutations
    tosolve   : boolean;   -- on when the system has to be solved
    torandq   : boolean;   -- on when random symmetric start system is solved

    procedure Read_Symmetry
                ( n : in natural32; pv,pw,fv,fw : in out List_of_Permutations;
                  fail : out boolean ) is

    -- DESCRIPTON :
    --   Reads and builds representations of the symmetry groups.

    -- ON ENTRY :
    --   n         dimension;

    -- ON RETURN :
    --   pv        permutation symmetry on the unknowns;
    --   pw        effect of the group actions in pv on the equations;
    --   fv        contains pv, plus eventually sign symmetry;
    --   fw        effect of the group actions in fv on the equations;
    --   fail      when the polynomial system is not symmetric.

      ans : character;
      nb : natural32;
      pg,fg : List_of_Permutations;

    begin
      Read_Permutation_Group(n,pg,pv,allperms);
      put("Is there a sign symmetry to take into account ? (y/n) ");
      Ask_Yes_or_No(ans);
      signsym := (ans = 'y');
      if signsym then
        put("Is the system invariant under all changes of signs ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          allsigns := true;
          fv := pv;
        else
          allsigns := false;
          signsym := false;  -- fv will contain these permutations
          allperms := false; -- fv will be used for the generating solutions
          put("The sign inversion of all elements is represented as ");
          for i in 1..n loop
            put('-');
            declare
              sb : Symbol;
            begin
              sb := (sb'range => ' ');
              sb := Symbol_Table.get(i);
              Symbol_Table_io.put(sb); put(" ");
            end;
          end loop;
          new_line;
          put("Give the number of generating elements in the group : ");
          Read_Natural(nb);
          put("Give "); put(nb,1);
          put_line(" vector representations of the generating elements :");
          Symbolic_Symmetry_Group_io.get(fg,n,nb);
          put("Do you want the generation of the group ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'y'
           then fv := Generate(Union(fg,pv));
           else fv := Union(fg,pv);
          end if;
        end if;
      else
        allsigns := false;
        fv := pv;
      end if;
      new_line(file);
      put_line(file,"THE SYMMETRY GROUP :");
      new_line(file);
      put_line(file,"v:"); Symbolic_Symmetry_Group_io.put(file,fv);
      new_line(file);
      Act(pv,p,pw,fail,invariant,equivaria);
      if not Is_Null(fg)
       then Act(fv,p,fw,fail,invariant,equivaria);
       else fw := pw;
      end if;
      new_line(file);
      put_line(file,"w:"); Symmetry_Group_io.put(file,fw); new_line(file);
      if allsigns then
        put_line(file,"The system is invariant under all changes of signs.");
      end if;
    end Read_Symmetry;

    procedure Data_Management
                 ( points : in out Arrays_of_Integer_Vector_Lists.
                                   Array_of_Lists;
                   pv,pw,fv,fw : in List_of_Permutations; fltlif : out boolean;
                   mix : out Standard_Integer_Vectors.Link_to_Vector;
                   imixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                   permp : out Poly_Sys ) is

    -- DESCRIPTION :
    --   Allows to read in a mixed subdivision, determines type of mixture.

    -- ON ENTRY :
    --   points    supports;
    --   pv        permutation symmetry on the unknowns;
    --   pw        effect of the group actions in pv on the equations.
    --   fv        contains pv, eventually with sign symmetry;
    --   fw        effect of the group actions in fv on the equations.

    -- ON RETURN :
    --   fltlif    true when floating-point lifting is used, false otherwise;
    --   mix       type of mixture;
    --   points    ordered according to mix;
    --   imixsub   integer mixed subdivision;
    --   permp     system ordered according to mix.

      use Standard_Integer_Vectors;
      use Arrays_of_Integer_Vector_Lists;
      use Integer_Mixed_Subdivisions;

      tmpmix,perm : Link_to_Vector;
      lifted_points : Array_of_Lists(p'range);
      m : natural32;
      mixsub : Mixed_Subdivision;
      ans : character;

      procedure Read_Subdivision is

        insubft : file_type;
        nn,bkk : natural32;

      begin
        put_line("Reading the name of the input file.");
        Read_Name_and_Open_File(insubft);
        get(insubft,nn,m,tmpmix,mixsub);
        Close(insubft);
        new_line(file);
        put_line(file,"MIXED SUBDIVISION :");
        new_line(file);
        put(file,nn,tmpmix.all,mixsub,bkk);
        new_line(file);
        compmisu := false;
        fltlif := false;
        imixsub := mixsub;
        permp := p;
      exception
        when DATA_ERROR
                => put_line("Data not in correct format.  Will ignore it...");
                   Close(insubft);
      end Read_Subdivision;

    begin
      new_line;
      put("Do you have already a mixed subdivision ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Read_Subdivision;
       else compmisu := true;
      end if;
      if compmisu then
        put("Do you want to enforce the type of mixture ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          put("Give the number of different supports : ");
          Read_Natural(m);
          put("Give the vector with occurencies : "); get(m,tmpmix);
          permp := p;
        else
          Compute_Mixture(points,tmpmix,perm);
          permp := Permute(p,perm); Clear(perm);
        end if;
        put("Do you want to have the subdivision on separate file ? (y/n) ");
        Ask_Yes_or_No(ans);
        misufile := (ans = 'y');
        if misufile
         then put_line("Reading the name of the output file.");
              Read_Name_and_Create_File(outsubft);
        end if;
        new_line(file);
        put_line(file,"THE TYPE OF MIXTURE :");
        new_line(file);
        put(file,"The number of different supports : ");
        put(file,tmpmix'last,1);
        new_line(file);
        put(file,"Vector indicating the occurrences : "); put(file,tmpmix);
        new_line(file);
      else
        misufile := false;
      end if;
      mix := tmpmix;
    end Data_Management;

    procedure Write_Generating_Cells
                 ( subfile : in file_type; n : in integer32;
                   mix : in Standard_Integer_Vectors.Vector;
                   lifted : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                   mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   Writes the list of generating cells on file, as the subdivision
    --   of the cell that contains all lifted points.
    --   By doing so, we will have no troubles recovering the lifting.

      use Arrays_of_Integer_Vector_Lists;
      use Integer_Mixed_Subdivisions;
  
      genmic : Mixed_Cell;
      genmixsub : Mixed_Subdivision;

    begin
      genmic.nor := new Standard_Integer_Vectors.Vector'(1..n+1 => 0);
      genmic.pts := new Array_of_Lists'(lifted);
      genmic.sub := new Mixed_Subdivision'(mixsub);
      Construct(genmic,genmixsub);
      put(subfile,natural32(n),mix,genmixsub);
    end Write_Generating_Cells;

    procedure Write_Generating_Cells
                ( subfile : in file_type; n : in integer32;
                  mix : in Standard_Integer_Vectors.Vector;
                  lifted : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                  mixsub : in Floating_Mixed_Subdivisions.Mixed_Subdivision ) is
 
    -- DESCRIPTION :
    --   Writes the list of generating cells on file, as the subdivision
    --   of the cell that contains all lifted points.
    --   By doing so, we will have no troubles recovering the lifting.

      use Arrays_of_Floating_Vector_Lists;
      use Floating_Mixed_Subdivisions;

      genmic : Mixed_Cell;
      genmixsub : Mixed_Subdivision;

    begin
      genmic.nor
        := new Standard_Floating_Vectors.Vector'(1..n+1 => 0.0);
      genmic.pts := new Array_of_Lists'(lifted);
      genmic.sub := new Mixed_Subdivision'(mixsub);
      Construct(genmic,genmixsub);
      put(subfile,natural32(n),mix,genmixsub);
    end Write_Generating_Cells;

    procedure Integer_Automatic_Lift_Orbits
                ( file : in file_type; norb : in natural32;
                  orbits
                    : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is

    -- DESCRIPTION :
    --   Random integer numbers are generated and given to each orbit.

      lower,upper : integer32 := 0;

    begin
      new_line(file);
      put_line(file,"AUTOMATIC RANDOM INTEGER LIFTING :");
      new_line(file);
      put("Give lower bound for random lifting : "); Read_Integer(lower);
      put("Give upper bound for random lifting : "); Read_Integer(upper);
      put(file,"  Lower bound for random lifting : ");
      put(file,lower,1); new_line(file);
      put(file,"  Upper bound for random lifting : ");
      put(file,upper,1); new_line(file);
      Integer_Random_Lift_Orbits(orbits,norb,lower,upper);
    end Integer_Automatic_Lift_Orbits;

    procedure Float_Automatic_Lift_Orbits
                ( file : in file_type; norb : in natural32;
                  orbits : in out Arrays_of_Floating_Vector_Lists.
                                  Array_of_Lists ) is

    -- DESCRIPTION :
    --   Random floating-point numbers are generated and given to each orbit.

      lower,upper : double_float;
  
    begin
      new_line(file);
      put_line(file,"AUTOMATIC RANDOM FLOATING-POINT LIFTING :");
      new_line(file);
      put("Give lower bound for random lifting : "); Read_Double_Float(lower);
      put("Give upper bound for random lifting : "); Read_Double_Float(upper);
      put(file,"  Lower bound for random lifting : ");
      put(file,lower,2,3,3); new_line(file);
      put(file,"  Upper bound for random lifting : ");
      put(file,upper,2,3,3); new_line(file);
      Float_Random_Lift_Orbits(orbits,norb,lower,upper);
    end Float_Automatic_Lift_Orbits;

    procedure Integer_Manual_Lift_Orbits
                ( file : in file_type; norb : in natural32;
                  orbits
                    : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is

    -- DESCRIPTION :
    --   The user can give integer lifting values for every orbit.

      rv : Standard_Integer_Vectors.Vector(1..integer32(norb))
         := (1..integer32(norb) => 0);

    begin
      new_line(file);
      put_line(file,"MANUAL INTEGER LIFTING :");
      new_line(file);
      put("Reading "); put(norb,1);
      put_line(" integer numbers to lift orbits");
      for i in rv'range loop
        put("  give lifting for orbit "); put(i,1); put(" : ");
        Read_Integer(rv(i));
      end loop;
      put_line(file,"  Lifting vector supplied by user :");
      put(file,rv); new_line(file);
      Integer_Lift_Orbits(orbits,rv);
    end Integer_Manual_Lift_Orbits;

    procedure Float_Manual_Lift_Orbits
                 ( file : in file_type; norb : in natural32;
                   orbits : in out Arrays_of_Floating_Vector_Lists.
                                   Array_of_Lists ) is

    -- DESCRIPTION :
    --   The user can give integer lifting values for every orbit.

      use Standard_Floating_Vectors;
      rv : Standard_Floating_Vectors.Vector(1..integer32(norb))
         := (1..integer32(norb) => 0.0);

    begin
      new_line(file);
      put_line(file,"MANUAL FLOATING-POINT LIFTING  :");
      new_line(file);
      put("Reading "); put(norb,1);
      put_line(" floating-point numbers to lift orbits");
      for i in rv'range loop
        put("  give lifting for orbit "); put(i,1); put(" : ");
        Read_Double_Float(rv(i));
      end loop;
      put_line(file,"  Lifting vector supplied by user :");
      Standard_Floating_Vectors_io.put(file,rv); new_line(file);
      Float_Lift_Orbits(orbits,rv);
    end Float_Manual_Lift_Orbits;

    procedure Classify_and_Lift_Orbits
                ( file : in file_type;
                  mix : in Standard_Integer_Vectors.Vector;
                  points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                  pv,pw : in List_of_Permutations; fltlif : out boolean;
                  ilft : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                  flft : out Arrays_of_Floating_Vector_Lists.Array_of_Lists ) is

    -- DESCRIPTION :
    --   Classifies the points into orbits and lifts orbits.

    -- ON ENTRY :
    --   file       output file;
    --   mix        type of mixture;
    --   points     support sets;
    --   pv         representation of permutation symmetry;
    --   pw         effect of group actions on the system.

    -- ON RETURN :
    --   fltlif     true when floating-point lifting, false otherwise;
    --   ilft       integer-valued lifted supports;
    --   flft       floating-point lifted supports;

      ans : character;
      norb : natural32;
      cnt : integer32;
      orbits : Arrays_of_Integer_Vector_Lists.Array_of_Lists(points'range);
      fltorb : Arrays_of_Floating_Vector_Lists.Array_of_Lists(points'range);

    begin
      new_line(file);
      put_line(file,"CLASSIFICATION OF POINTS INTO ORBITS :");
      new_line(file);
      Classify_Orbits(points,mix,pv,pw,norb,orbits);
      cnt := orbits'first;
      new_line; put("Classified orbits,");
      put_line(" last coordinate of vector is orbit number : ");
      for k in mix'range loop
        put("support no. "); put(cnt,1); put_line(" :"); put(orbits(cnt));
        put(file,"support no. "); put(file,cnt,1);
        put_line(file," :"); put(file,orbits(cnt));
        cnt := cnt + mix(k);
      end loop;
      new_line;
      put_line("MENU for Lifting of Orbits");
      put("  1. Integer Automatic : ");
      put(norb,1); put_line(" random integer numbers as lifting.");
      put("  2.         Manual    : ");
      put("you can give "); put(norb,1); put_line(" integer lifting values.");
      put("  3. Float   Automatic : ");
      put(norb,1); put_line(" random floating-point numbers as lifting.");
      put("  4.         Manual    : ");
      put("you can give "); put(norb,1);
      put_line(" floating-point lifting values.");
      put("Type 1, 2, 3, or 4 to select lifting : ");
      Ask_Alternative(ans,"1234");
      case ans is
        when '1' => Integer_Automatic_Lift_Orbits(file,norb,orbits);
                    ilft := orbits; fltlif := false;
        when '2' => Integer_Manual_Lift_Orbits(file,norb,orbits);
                    ilft := orbits; fltlif := false;
        when '3' => fltorb := Convert(orbits);
                    Float_Automatic_Lift_Orbits(file,norb,fltorb);
                    flft := fltorb; fltlif := true;
        when '4' => fltorb := Convert(orbits);
                    Float_Manual_Lift_Orbits(file,norb,fltorb);
                    flft := fltorb; fltlif := true;
        when others => null;
      end case;
    end Classify_and_Lift_Orbits;

    procedure Integer_Prune_for_Mixed_Cells
                 ( file : in file_type; n : in integer32;
                   mix : in Standard_Integer_Vectors.Vector;
                   lifpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                   pv : in List_of_Permutations;
                   mixsub
                     : in out Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   Given the lifted supports, the mixed cells will be computed.

    -- ON ENTRY :
    --   file      output file;
    --   n         dimension;
    --   mix       type of mixture;
    --   lifpts    lifted supports;
    --   pv        representation of permutation group.

    -- ON RETURN :
    --   mixsub    mixed cells;

      use Integer_Faces_of_Polytope;
      use Integer_Mixed_Subdivisions;

      fa : Array_of_Faces(mix'range);
      nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                    := (mix'range => 0.0);
      timer : Timing_Widget;

    begin
      tstart(timer);
      for k in fa'range loop
        fa(k) := Create_Lower(mix(k),n+1,lifpts(k));
      end loop;
      if invariant then
        if allperms
         then fa(fa'first) := Generating_Lifted_Faces(fa(fa'first));
         else fa(fa'first) := Generating_Lifted_faces(pv,fa(fa'first));
        end if;
      end if;
      tstop(timer);
      new_line(file);
      put_line(file,"CARDINALITIES OF THE LIFTED FACES :");
      new_line(file);
      for i in fa'range loop
        put(file,"  # lifted "); put(file,mix(i),1);
        put(file,"-faces of polytope "); put(file,i,1);
        put(file," : "); put(file,Length_Of(fa(i)),1); new_line(file);
      end loop;
      new_line(file);
      print_times(file,timer,"Creation of the Lower Faces");
      new_line(file);
      tstart(timer);
      Create_CS(n,mix,fa,lifpts,nbsucc,nbfail,mixsub);
      tstop(timer);
      Pruning_Statistics(file,nbsucc,nbfail);
      new_line(file);
      print_times(file,timer,"Pruning for Mixed Cells");
    end Integer_Prune_for_Mixed_Cells;

    procedure Float_Prune_for_Mixed_Cells
                 ( file : in file_type; n : in integer32;
                   mix : in Standard_Integer_Vectors.Vector;
                   lifpts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                   pv : in List_of_Permutations;
                   mixsub
                     : in out Floating_Mixed_Subdivisions.Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   Given the lifted supports, the mixed cells will be computed.

    -- ON ENTRY :
    --   file      output file;
    --   n         dimension;
    --   mix       type of mixture;
    --   lifpts    lifted supports;
    --   pv        representation of permutation group.

    -- ON RETURN :
    --   mixsub    mixed cells;

      use Floating_Faces_of_Polytope;
      use Floating_Mixed_Subdivisions;

      tol : constant double_float := 1.0E-10;
      fa : Array_of_Faces(mix'range);
      nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                    := (mix'range => 0.0);
      timer : Timing_Widget;

    begin
      tstart(timer);
      for k in fa'range loop
        fa(k) := Create_Lower(mix(k),n+1,lifpts(k),tol);
      end loop;
--    if invariant
--     then if allperms
--           then fa(fa'first) := Generating_Lifted_Faces(fa(fa'first));
--           else fa(fa'first) := Generating_Lifted_faces(pv,fa(fa'first));
--          end if;
--    end if;
      tstop(timer);
      new_line(file);
      put_line(file,"CARDINALITIES OF THE LIFTED FACES :");
      new_line(file);
      for i in fa'range loop
        put(file,"  # lifted "); put(file,mix(i),1);
        put(file,"-faces of polytope "); put(file,i,1);
        put(file," : "); put(file,Length_Of(fa(i)),1); new_line(file);
      end loop;
      new_line(file);
      print_times(file,timer,"Creation of the Lower Faces");
      new_line(file);
      tstart(timer);
      Create(n,mix,fa,lifpts,tol,nbsucc,nbfail,mixsub);
      tstop(timer);
      Pruning_Statistics(file,nbsucc,nbfail);
      new_line(file);
      print_times(file,timer,"Pruning for Mixed Cells");
    end Float_Prune_for_Mixed_Cells;

    procedure Compute_Mixed_Volume
                 ( file : in file_type; n : in integer32;
                   mix : in Standard_Integer_Vectors.Vector;
                   lifpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                   mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                   mv : out natural32 ) is
  
    -- DESCRIPTION :
    --   Computes the mixed volume of the supports in lifpts and checks
    --   on the mixed cells in the subdivision.
    --   Results are written on the output file.

    -- ON ENTRY :
    --   file      output file;
    --   n         dimension;
    --   mix       type of mixture;
    --   lifpts    lifted supports;
    --   mixsub    list of mixed cells.

    -- ON RETURN :
    --   mixsub    can contain refinements of cells;
    --   mv        mixed volume.

      use Integer_Mixed_Subdivisions;
  
      timer : timing_widget;
      bkk : natural32;

    begin
      new_line(file);
      put_line(file,"THE LIFTED SUPPORTS :");
      new_line(file);
      put(file,lifpts);
      new_line(file);
      put_line(file,"VOLUME OF MIXED CELLS :");
      new_line(file);
      tstart(timer);
      put(file,natural32(n),mix,mixsub,bkk);
      tstop(timer);
      new_line(file);
      put(file,"The mixed volume : "); put(file,bkk,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Volume computation of mixed cells");
      new_line(file);
      mixsub := Create(lifpts,mixsub);
      new_line(file);
      put_line(file,"COMPUTING AGAIN AFTER CHECKING :");
      new_line(file);
      tstart(timer);
      put(file,natural32(n),mix,mixsub,bkk);
      tstop(timer);
      mv := bkk;
      put(file,"The mixed volume : "); put(file,bkk,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Checking the Mixed Volume Computation");
    end Compute_Mixed_Volume;

    procedure Compute_Mixed_Volume
                ( file : in file_type; n : in integer32;
                  mix : in Standard_Integer_Vectors.Vector;
                  lifpts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                  mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                  mv : out natural32 ) is

    -- DESCRIPTION :
    --   Computes the mixed volume of the supports in lifpts.
    --   Results are written on the output file.

    -- ON ENTRY :
    --   file      output file;
    --   n         dimension;
    --   mix       type of mixture;
    --   lifpts    lifted supports;
    --   mixsub    list of mixed cells.

    -- ON RETURN :
    --   mixsub    can contain refinements of cells;
    --   mv        mixed volume.

      timer : timing_widget;
      bkk : natural32;

    begin
      new_line(file);
      put_line(file,"THE LIFTED SUPPORTS :");
      new_line(file);
      put(file,lifpts);
      new_line(file);
      put_line(file,"VOLUME OF MIXED CELLS :");
      new_line(file);
      tstart(timer);
      put(file,natural32(n),mix,mixsub,bkk);
      mv := bkk;
      tstop(timer);
      new_line(file);
      put(file,"The mixed volume : "); put(file,bkk,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Volume computation of mixed cells");
    end Compute_Mixed_Volume;

    procedure Generating_Mixed_Cells
                   ( file : in file_type; n : in integer32;
                     points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                     pv,pw : in List_of_Permutations;
                     mix : in Standard_Integer_Vectors.Vector;
                     lorb : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                     lifted : out Arrays_of_Integer_Vector_Lists.
                                  Link_to_Array_of_Lists;
                     gensub : in out Integer_Mixed_Subdivisions.
                                     Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   Computes a list of generating mixed cells.

    -- ON ENTRY :
    --   file        output file;
    --   n           dimension;
    --   points      support sets;
    --   pv          representation of permutation symmetry;
    --   pw          effect of group actions in pv on the system;
    --   mix         type of mixture;
    --   lorb        classified and lifted supports.

    -- ON RETURN :
    --   lifted      lifted points;
    --   gensub      generating cells.

      use Arrays_of_Integer_Vector_Lists;
      use Integer_Mixed_Subdivisions;

      lifpts : Array_of_Lists(mix'range);
      mixsub,genmixsub : Mixed_Subdivision;
      bkk : natural32;
      index : integer32;
      timer : timing_widget;

    begin
      if compmisu then
        index := lorb'first;
        for k in lifpts'range loop
          lifpts(k) := lorb(index);
          index := index + mix(k);
        end loop;
        Integer_Prune_for_Mixed_Cells(file,n,mix,lifpts,pv,mixsub);
        Compute_Mixed_Volume(file,n,mix,lifpts,mixsub,mv);
        tstart(timer);
        if allperms
         then genmixsub := Generating_Cells(mixsub);
         else genmixsub := Generating_Cells(pv,pw,mix,mixsub);
        end if;
        new_line(file);
        put_line(file,"THE GENERATING CELLS :");
        new_line(file);
        put(file,natural32(n),mix,genmixsub,bkk);
        tstop(timer);
        if misufile
         then Write_Generating_Cells(outsubft,n,mix,lifpts,genmixsub);
        end if;
        new_line(file);
        put(file,"Number of generating solutions : "); put(file,bkk,1);
        new_line(file); new_line(file);
        print_times(file,timer,"Computing generating cells");
        gensub := genmixsub;
      else
        lifpts := Induced_Lifting(n,mix,points,gensub);
        new_line(file);
        put_line(file,"THE LIFTED SUPPORTS :");
        new_line(file);
        put(file,lifpts);
        new_line(file);
      end if;
      lifted := new Array_of_Lists'(lifpts);
    end Generating_Mixed_Cells;

    procedure Generating_Mixed_Cells
                   ( file : in file_type; n : in integer32;
                     points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                     pv,pw : in List_of_Permutations;
                     mix : in Standard_Integer_Vectors.Vector;
                     lorb : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                     lifted : out Arrays_of_Floating_Vector_Lists.
                                  Link_to_Array_of_Lists;
                     gensub : out Floating_Mixed_Subdivisions.
                                  Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   Computes a list of generating mixed cells.
  
    -- ON ENTRY :
    --   file        output file;
    --   n           dimension;
    --   points      support sets;
    --   pv          representation of permutation symmetry;
    --   pw          effect of group actions in pv on the system;
    --   mix         type of mixture;
    --   lorb        classified and lifted supports.

      use Arrays_of_Floating_Vector_Lists;
      use Floating_Mixed_Subdivisions;

      fpoints,lifpts : Array_of_Lists(mix'range);
      mixsub,genmixsub : Mixed_Subdivision;
      bkk : natural32;
      index : integer32;
      timer : timing_widget;

    begin
      if compmisu then
        index := lorb'first;
        for k in lifpts'range loop
          lifpts(k) := lorb(index);
          index := index + mix(k);
        end loop;
        Float_Prune_for_Mixed_Cells(file,n,mix,lifpts,pv,mixsub);
        Compute_Mixed_Volume(file,n,mix,lifpts,mixsub,mv);
      else
        fpoints := Convert(points);
        lifpts := Induced_Lifting(n,mix,fpoints,mixsub);
        new_line(file);
        put_line(file,"THE LIFTED SUPPORTS : ");
        new_line(file);
        put(file,lifpts);
        new_line(file);
      end if;
      tstart(timer);
      if allperms
       then genmixsub := Generating_Cells(mixsub);
       else genmixsub := Generating_Cells(pv,pw,mix,mixsub);
      end if;
      new_line(file);
      put_line(file,"THE GENERATING CELLS  :");
      new_line(file);
      put(file,natural32(n),mix,genmixsub,bkk);
      tstop(timer);
      if misufile
       then Write_Generating_Cells(outsubft,n,mix,lifpts,genmixsub);
      end if;
      new_line(file);
      put(file,"Number of generating solutions : "); put(file,bkk,1);
      new_line(file); new_line(file);
      print_times(file,timer,"Computing generating cells");
      lifted := new Array_of_Lists'(lifpts);
      gensub := genmixsub;
    end Generating_Mixed_Cells;

    procedure Settings_for_Polyhedral_Continuation is

    -- DESCRIPTION :
    --   Displays the menu and allows the user to set the parameters.

      ans : character;
      oc : natural32;

    begin
      new_line;
      put_line("MENU for Symmetric Polyhedral Continuation : ");
      put_line("  0. No polyhedral continuation, leave the menu.");
      put_line("  1. Solve given system by polyhedral continuation.");
      put_line("  2. Create and solve random coefficient system.");
      put("Type 0,1, or 2 to choose : "); Ask_Alternative(ans,"012");
      tosolve := not (ans = '0');
      torandq := (ans = '2');
      if tosolve then
        if torandq then
          put_line("Reading a name of a file to write start system on.");
          Read_Name_and_Create_File(gft);
        else 
          new_line;
          put_line("Reading a name of a file to write start solutions on.");
          Read_Name_and_Create_File(solsft);
        end if;
        new_line;
        Driver_for_Continuation_Parameters(file);
        new_line;
        Driver_for_Process_io(file,oc);
      end if;
      if byebye then
        new_line;
        put_line("No more input expected.  See output file for results.");
        new_line;
      else
        new_line;
        put_line("Starting the Computations ...");
        new_line;
      end if;
    end Settings_for_Polyhedral_Continuation;

    procedure Symmetric_Polyhedral_Continuation
                  ( file : in file_type;
                    pp : in Poly_Sys; fv,fw : List_of_Permutations;
                    n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                    lifpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                    gensub : in Integer_Mixed_Subdivisions.
                                Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   Constructs and solves a symmetric random-coefficient start system.  

    -- ON ENTRY :
    --   file       output file;
    --   pp         polynomial system, ordered according to mixture;
    --   fv         representation of symmetry group;
    --   fw         effect on group actions on system;
    --   n          dimension;
    --   mix        type of mixture;
    --   lifpts     lifted points;
    --   gensub     generating mixed cells.

      use Arrays_of_Integer_Vector_Lists;
      use Integer_Mixed_Subdivisions;
  
      lq,lifted_lq : Laur_Sys(p'range);
      lp : Laur_Sys(p'range);   
      qq : Poly_Sys(p'range);
      qqsols : Solution_List;
      timer : timing_widget;

    begin
      new_line(file);
      put_line(file,"COMPUTING THE GENERATING SOLUTIONS :");
      new_line(file);
      tstart(timer);
      lp := Polynomial_to_Laurent_System(pp);
      if torandq
       then lq := Symmetric_Randomize(lp,fv,fw);
            qq := Laurent_to_Polynomial_System(lq); q := qq;
            put(gft,qq);
       else lq := lp; q := p;
      end if;
      lifted_lq := Perform_Lifting(n,mix,lifpts,lq);
      if allperms
       then qqsols := Symmetric_Mixed_Solve
                        (file,signsym,lifted_lq,gensub,n,mix);
       else qqsols := Symmetric_Mixed_Solve
                        (file,fv,signsym,lifted_lq,gensub,n,mix);
      end if;
      tstop(timer);
      new_line(file);
      put(file,qqsols);
      if torandq
       then new_line(gft);
            put_line(gft,"THE GENERATING SOLUTIONS :"); new_line(gft);
            put(gft,Length_Of(qqsols),natural32(n),qqsols);
            Close(gft);
       else put(solsft,Length_Of(qqsols),natural32(n),qqsols);
            Close(solsft);
      end if;
      new_line(file);
      print_times(file,timer,"Symmetric polyhedral continuation");
      qsols := qqsols;
    end Symmetric_Polyhedral_Continuation;

    procedure Main_Driver is

      n : constant integer32 := p'length;
      timer : timing_widget;
      notsym : boolean;
      pv,pw,fv,fw : List_of_Permutations;
      igencells : Integer_Mixed_Subdivisions.Mixed_Subdivision;
      fgencells : Floating_Mixed_Subdivisions.Mixed_Subdivision;
      mix : Standard_Integer_Vectors.Link_to_Vector;
      permp : Poly_Sys(p'range);
      points,iliforb : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
      fliforb : Arrays_of_Floating_Vector_Lists.Array_of_Lists(p'range);
      fltlif : boolean;
      ilifpts : Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
      flifpts : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;

    begin
      new_line; put_line(welcome);
      tstart(timer);
      Read_Symmetry(natural32(n),pv,pw,fv,fw,notsym);
      if notsym then
        put_line("The given system is not symmetric !");
        put_line(file,"The given system is not symmetric !");
      else
        points := Create(p);
        Data_Management(points,pv,pw,fv,fw,fltlif,mix,igencells,permp);
        if compmisu
         then Classify_and_Lift_Orbits
                (file,mix.all,points,pv,pw,fltlif,iliforb,fliforb);
        end if;
        Settings_for_Polyhedral_Continuation;
        tstart(timer);
        if fltlif then
          Generating_Mixed_Cells
            (file,n,points,pv,pw,mix.all,fliforb,flifpts,fgencells);
        else
          Generating_Mixed_Cells
            (file,n,points,pv,pw,mix.all,iliforb,ilifpts,igencells);
          if tosolve
           then Symmetric_Polyhedral_Continuation
                  (file,permp,fv,fw,n,mix.all,ilifpts.all,igencells);
          end if;
        end if;
      end if;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"All Computations");
    end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_Symmetric_Mixed_Volume_Computation;

  procedure Driver_for_Symmetric_Mixed_Volume_Computation
               ( file : in file_type; p : in Laur_Sys; byebye : in boolean;
                 q : out Laur_Sys; qsols : out Solution_List;
                 mv : out natural32 ) is

    pp,pq : Poly_Sys(p'range);

  begin
    pp := Laurent_to_Polynomial_System(p);
    Driver_for_Symmetric_Mixed_Volume_Computation(file,pp,byebye,pq,qsols,mv);
    q := Polynomial_to_Laurent_System(pq);
  end Driver_for_Symmetric_Mixed_Volume_Computation;

end Drivers_for_Symmetric_Lifting;
