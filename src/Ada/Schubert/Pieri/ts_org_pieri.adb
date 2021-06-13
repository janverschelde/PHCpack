with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_VecMats_io;        use Standard_Complex_VecMats_io;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Matrix_Indeterminates;              use Matrix_Indeterminates;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Bracket_Monomials;                  use Bracket_Monomials;
with Bracket_Monomials_io;               use Bracket_Monomials_io;
with Bracket_Expansions;                 use Bracket_Expansions;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Polynomials_io;    use Standard_Bracket_Polynomials_io;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;
with Pieri_Trees,Pieri_Trees_io;         use Pieri_Trees,Pieri_Trees_io;
with Pieri_Root_Counts;                  use Pieri_Root_Counts;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Numeric_Minor_Equations;            use Numeric_Minor_Equations;
with Solve_Pieri_Leaves;
with Specialization_of_Planes;           use Specialization_of_Planes;
with Pieri_Deformations;                 use Pieri_Deformations;

procedure ts_org_pieri is

-- DESCRIPTION :
--   This program executes the original implementation of the Pieri homotopy
--   algorithm, implemented strictly along the lines of the description in
--   the paper "Numerical Schubert Calculus" of Birkett Huber, Frank Sottile
--   and Bernd Sturmfels.

  function Maximum ( n1,n2 : natural32 ) return natural32 is
  begin
    if n1 >= n2
     then return n1;
     else return n2;
    end if;
  end Maximum;

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

-- DISPLAYING THE REPRESENTATIONS OF THE PLANES :

  procedure Display_Polynomial_Pattern
              ( file : in file_type; n : in integer32; b1,b2 : in Bracket ) is

  -- DESCRIPTION :
  --   Displays the pattern by a polynomial matrix.

    pm : Standard_Complex_Poly_Matrices.Matrix(1..n,b1'range)
       := Schubert_Pattern(natural32(n),b1,b2);

  begin
    put(file,pm);
    Standard_Complex_Poly_Matrices.Clear(pm);
  end Display_Polynomial_Pattern;

-- AUXILIARIES TO SET UP EQUATIONS :

  procedure Expand_Minors ( file : in file_type;
                            mat : in Standard_Complex_Poly_Matrices.Matrix;
                            bm : in Bracket_Monomial ) is

  -- DESCRIPTION :
  --   Expands the quadratic bracket monomial.

    first : boolean := true;
    lb : Link_to_Bracket;

    procedure Visit_Bracket ( b : in Bracket; continue : out boolean ) is

      p : Poly := Null_Poly;

    begin
      if first
       then lb := new Bracket'(b);
            first := false;
       else p := Expanded_Minor(mat,b);
            if p /= Null_Poly
             then put(file,lb.all);
                  put(file,"*("); put(file,p); put(file,")");
            end if;
            Clear(lb); Clear(p);
      end if;
      continue := true;
    end Visit_Bracket;
    procedure Visit_Brackets is new Enumerate_Brackets(Visit_Bracket);

  begin
    Visit_Brackets(bm);
  end Expand_Minors;

  procedure Expand_Minors ( file : in file_type;
                            mat : in Standard_Complex_Poly_Matrices.Matrix;
                            bp : Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Writes the expansion of the matrix, using the bracket polynomial which
  --   is a list of quadratic monomials that represent the Laplace expansion.

    procedure Visit_Term ( t : in Bracket_Term; continue : out boolean ) is
    begin
      if REAL_PART(t.coeff) > 0.0
       then put("+");
       else put("-");
      end if;
      Expand_Minors(file,mat,t.monom);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Enumerate_Terms(Visit_Term);

  begin
    Visit_Terms(bp);
  end Expand_Minors;

  procedure Pieri_Equations
              ( file : in file_type; n,d : in integer32;
                bs : in Bracket_System; b1,b2 : in Bracket ) is

  -- DESCRIPTION :
  --   Writes the Pieri equations corresponding to the pair of brackets.

    cffmat : constant Standard_Complex_Matrices.Matrix(1..n,1..(n-d))
           := Random_Matrix(natural32(n),natural32(n-d));
    polmat : Standard_Complex_Poly_Matrices.Matrix(1..n,1..d)
           := Schubert_Pattern(natural32(n),b1,b2);
    sys : Poly_Sys(bs'first+1..bs'last);
    sol : Standard_Complex_Matrices.Matrix(1..n,1..d);

  begin
    put(file,"Plane X for (");
    put(file,b1); put(file,","); put(file,b2); put_line(file,") :");
    put(file,polmat);
    put_line(file,"The system with expanded minors of X : ");
    for i in 1..bs'last loop
      Expand_Minors(file,polmat,bs(i)); new_line(file);
    end loop;
   -- put("Give a "); put(n,1); put("x"); put(n-d,1);
   -- put_line("-matrix of complex numbers : ");
   -- get(cffmat);
   -- put("Your "); put(n-d,1); put_line("-plane : "); put(cffmat);
   -- put("A random "); put(n-d,1); put_line("-plane : "); put(cffmat);
    put("The minors evaluated at a random ");
    put(n-d,1); put_line("-plane :");
    sys := Expanded_Minors(cffmat,polmat,bs);
    put_line(sys);
    Standard_Complex_Poly_Matrices.Clear(polmat);
    if Pieri_Condition(natural32(n),b1,b2) then
      put("The "); put(d,1); put_line("-plane at the leaves :");
      sol := Solve_Pieri_Leaves(Standard_Output,b1,b2,cffmat); put(sol);
      put_line("The solution evaluated at the system : ");
      put(Evaluate(sys,sol)); new_line;
    else
      put_line("Pair of leaves does not satisfy Pieri's condition.");
    end if;
  end Pieri_Equations;

  procedure Expand_Minors ( n,d : in integer32; bs : in Bracket_System ) is

  -- DESCRIPTION :
  --   Expands the minors to obtain a symbolic formulation of the equations.

    b1,b2 : Bracket(1..d);

  begin
    put("Give 1st bracket : "); get(b1);
    put("Give 2nd bracket : "); get(b2);
    Pieri_Equations(Standard_Output,n,d,bs,b1,b2);
  end Expand_Minors;

  procedure Pieri_Equations_for_Paired_Chains
              ( file : in file_type; n,d : in integer32;
                bs : in Bracket_System; b1,b2 : in bracket_Array ) is

  -- DESCRIPTION :
  --   Writes the equations for each node along a pair of chains.

    maxlen : constant natural32
           := Maximum(natural32(b1'last),natural32(b2'last));
    lb1,lb2 : Link_to_Bracket;

  begin
    for i in b1'first..integer32(maxlen) loop
      if i <= b1'last
       then lb1 := b1(i);
       else lb1 := b1(b1'last);
      end if;
      if i <= b2'last
       then lb2 := b2(i);
       else lb2 := b2(b2'last);
      end if;
      Pieri_Equations(file,n,d,bs,lb1.all,lb2.all);
    end loop;
  end Pieri_Equations_for_Paired_Chains;

  procedure Write_Pieri_Equations
              ( file : in file_type; n,d : in integer32;
                t1,t2 : in Pieri_Tree; bs : in Bracket_System ) is

  -- DESCRIPTION :
  --   Writes the Pieri Equations for all pairs of chains.

    procedure Visit_Chain ( b1,b2 : in Bracket_Array; cont : out boolean ) is
    begin
      Pieri_Equations_for_Paired_Chains(file,n,d,bs,b1,b2);
      cont := true;
    end Visit_Chain;
    procedure Visit_Chains is new Enumerate_Paired_Chains(Visit_Chain);

  begin
    Visit_Chains(t1,t2);
  end Write_Pieri_Equations;

  procedure Write_Pieri_Equations
              ( n,d : in integer32; t1,t2 : in Pieri_Tree ) is

    k,kd : integer32 := 0;
    bm : Bracket_Monomial;
    file : file_type;

  begin
    new_line; skip_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    put("Give k to determine (m-k+1)-plane : "); get(k);
    kd := n-k+1;
    bm := Maximal_Minors(natural32(n),natural32(kd));      -- because n = m+p
    put(file,"All maximal minors : "); put(file,bm); new_line(file);
    declare
      bs : constant Bracket_System
         := Minor_Equations(natural32(kd),natural32(kd-d),bm);
    begin
      put_line(file,"The generic equation in the Laplace expansion : ");
      put(file,bs(0));
      put_line(file,"The specific equations in the system : ");
      for i in 1..bs'last loop
        put(file,bs(i));
      end loop;
      Write_Pieri_Equations(file,n,d,t1,t2,bs);
    end;
    Close(file);
    Clear(bm);
  end Write_Pieri_Equations;

-- AUXILIARIES TO REPRESENT PIERI TREES :

  procedure Write_a_Chain ( file : in file_type; b : in Bracket_Array ) is
  begin
    put(b(b'first).all);
    for i in b'first+1..b'last loop
      put(" < "); put(b(i).all);
    end loop;
    new_line;
  end Write_a_Chain;

  procedure Write_Chains ( t : in Pieri_Tree ) is

    procedure Write_Chain ( b : in Bracket_Array; cont : out boolean ) is
    begin
      Write_a_Chain(Standard_Output,b);
      cont := true;
    end Write_Chain;
    procedure Write is new Enumerate_Chains(Write_Chain);

  begin
    Write(t);
  end Write_Chains;

-- AUXILIARIES TO COUNT THE ROOTS :

  procedure Write_Polynomial_Patterns
              ( file : in file_type;
                n : in integer32; b1,b2 : in Bracket_Array ) is

  -- DESCRIPTION :
  --   Writes the two chains in a paired fashion.  If they have unequal
  --   length, then the last element of the shortest chain appears repeated.

    maxlen : constant natural32
           := Maximum(natural32(b1'last),natural32(b2'last));
    lb1,lb2 : Link_to_Bracket;

  begin
    for i in b1'first..integer32(maxlen) loop
      put(file,"(");
      if i <= b1'last
       then lb1 := b1(i);
       else lb1 := b1(b1'last);
      end if;
      put(file,lb1.all); put(file,",");
      if i <= b2'last
       then lb2 := b2(i);
       else lb2 := b2(b2'last);
      end if;
      put(file,lb2.all); put_line(file,") has pattern : ");
      Display_Polynomial_Pattern(file,n,lb1.all,lb2.all);
    end loop;
  end Write_Polynomial_Patterns;

  procedure Write_Paired_Chain
              ( file : in file_type;
                n : in integer32; b1,b2 : in Bracket_Array ) is

  -- DESCRIPTION :
  --   Writes the two chains in a paired fashion.  If they have unequal
  --   length, then the last element of the shortest chain appears repeated.

    maxlen : constant natural32
           := Maximum(natural32(b1'last),natural32(b2'last));

  begin
    put(file,"("); put(file,b1(b1'first).all); put(file,",");
                   put(file,b2(b2'first).all); put(file,")");
    for i in b1'first+1..integer32(maxlen) loop
      put(file," < "); put(file,"(");
      if i <= b1'last
       then put(file,b1(i).all);
       else put(file,b1(b1'last).all);
      end if;
      put(file,",");
      if i <= b2'last
       then put(file,b2(i).all);
       else put(file,b2(b2'last).all);
      end if;
      put(file,")");
    end loop;
    new_line(file);
    Write_Polynomial_Patterns(Standard_Output,n,b1,b2);
    if Pieri_Condition(natural32(n),b1(b1'last).all,b2(b2'last).all)
     then put_line("Leaves satisfy Pieri's condition.");
     else put_line("Leaves do not satisfy Pieri's condition.");
    end if;
  end Write_Paired_Chain;

  procedure Write_Pieri_Chains
              ( n : in integer32; t1,t2 : in Pieri_Tree ) is

  -- DESCRIPTION :
  --   Enumerates all pairs of chains and checks Pieri's condition
  --   at the leaves.

    procedure Visit_Chain ( b1,b2 : in Bracket_Array; cont : out boolean ) is
    begin
      Write_Paired_Chain(Standard_Output,n,b1,b2);
      cont := true;
    end Visit_Chain;
    procedure Visit_Chains is new Enumerate_Paired_Chains(Visit_Chain);

  begin
    Visit_Chains(t1,t2);
  end Write_Pieri_Chains;

  function First_Standard_Plane
              ( n,m,r : integer32 ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the plane spanned by the first m+1-r standard basis vectors.

    res : Standard_Complex_Matrices.Matrix
        := Random_Matrix(natural32(n),natural32(m+1-r));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := Create(1.0);
         else res(i,j) := Create(0.0);
        end if;
      end loop;
    end loop;
    return res;
  end First_Standard_Plane;

  function Last_Standard_Plane
              ( n,m,r : integer32 ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the plane spanned by the first m+1-r standard basis vectors.

    res : Standard_Complex_Matrices.Matrix
        := Random_Matrix(natural32(n),natural32(m+1-r));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = res'last(2) + 1 - j
         then res(i,j) := Create(1.0);
         else res(i,j) := Create(0.0);
        end if;
      end loop;
    end loop;
    return res;
  end Last_Standard_Plane;

  function First_Random_Input_Sequence
             ( n,m,a : integer32; kp : Vector ) return VecMat is

  -- DESCRIPTION :
  --   Returns the first sequence of random input planes.  The first plane
  --   in the sequence is spanned by the first standard basis vectors.
  --   The dimensions of the planes are m+1-kp(i), for the appropriate i.

    res : VecMat(0..a-1);

  begin
    res(0) := new Standard_Complex_Matrices.Matrix'
                    (First_Standard_Plane(n,m,integer32(kp(1))));
    for i in 1..res'last loop
      res(i) := new Standard_Complex_Matrices.Matrix'
                      (Random_Matrix(natural32(n),
                                     natural32(m)+1-natural32(kp(i+1))));
    end loop;
    return res;
  end First_Random_Input_Sequence;

  function Second_Random_Input_Sequence
             ( n,m,a : integer32; kp : Vector ) return VecMat is

  -- DESCRIPTION :
  --   Returns the second sequence of random input planes.  The first plane
  --   in the sequence is spanned by the last standard basis vectors.

    res : VecMat(0..kp'last-a-2);

  begin
    res(0) := new Standard_Complex_Matrices.Matrix'
                    (Last_Standard_Plane(n,m,integer32(kp(1))));
    for i in 1..res'last loop
      res(i) := new Standard_Complex_Matrices.Matrix'
                      (Random_Matrix(natural32(n),
                                     natural32(m)+1-natural32(kp(a+1+i))));
    end loop;
    return res;
  end Second_Random_Input_Sequence;

  procedure Reallify ( c : in out Complex_Number ) is

  -- DESCRIPTION :
  --   Sets the imaginary part of c to zero.

  begin
    c := Create(REAL_PART(c),0.0);
  end Reallify;

  procedure Reallify ( m : in out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Sets the imaginary part of every entry in the matrix to zero.

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        Reallify(m(i,j));
      end loop;
    end loop;
  end Reallify;

  procedure Reallify ( v : in out VecMat ) is

  -- DESCRIPTION :
  --   Sets the imaginary part of every entry of every matrix in v to zero
  --   and makes the matrices orthonormal.

  begin
    for i in v'range loop
      Reallify(v(i).all);
      v(i).all := Orthogonalize(v(i).all);
    end loop;
  end Reallify;

  procedure Orthogonalize ( v : in out VecMat ) is

  -- DESCRIPTION :
  --   Orthonormalizes every matrix in the array.

  begin
    for i in v'range loop
      v(i).all := Orthogonalize(v(i).all);
    end loop;
  end Orthogonalize;

  procedure Set_Parameters ( file : in file_type; report : out boolean ) is

  -- DESCRIPTION :
  --   Interactive determination of the continuation and output parameters.

    oc : natural32;

  begin
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := not (oc = 0);
  end Set_Parameters;

  function Select_Pairs ( lps : List_of_Paired_Nodes )
                        return List_of_Paired_Nodes is

  -- DESCRIPTION :
  --   Returns a selection of the list of pairs.

    res,res_last : List_of_Paired_Nodes;
    k : integer32 := 0;
    tmp : List_of_Paired_Nodes := lps;

  begin
    put("Give the number of pairs : "); get(k);
    declare
      sel : Standard_Natural_Vectors.Vector(1..k);
      ind : integer32 := 1;
    begin
      put("Give an increasing sequence of "); put(k,1); put(" numbers : ");
      get(sel);
      for i in 1..Length_Of(lps) loop
        if i = sel(ind) then
          ind := ind+1;
          Append(res,res_last,Head_Of(tmp));
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end;
    skip_line;
    return res;
  end Select_Pairs;

  procedure put ( file : in file_type; lp : in List_of_Paired_Nodes ) is

    tmp : List_of_Paired_Nodes := lp;
    pnd : Paired_Nodes;

  begin
    put_line(file,"The pairs of leaves that satisfy Pieri's condition :");
    for i in 1..Length_Of(lp) loop
      pnd := Head_Of(tmp);
      put(file,"Leaf "); put(file,i,3); put(file," : ");
      put(file,"("); put(file,pnd.left.node); put(file,",");
      put(file,pnd.right.node); put_line(file,")");
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure Root_Count ( n,d,a : in integer32; kp : in Vector ) is

  -- DESCRIPTION :
  --   Set up of Pieri trees from a partition of the planes.

  -- ON ENTRY :
  --  n          the dimension of the space we are working in;
  --  d          dimension of the brackets, of the output planes;
  --  a          number of planes in the first set of the partition;
  --  kp         the dimensions of the input planes are m+1-kp(i),
  --             where m = n-d.

    file : file_type;
    timer : Timing_Widget;
    ans : character;
    m : constant integer32 := n-d;
    v1 : constant Vector(0..a-1) := kp(1..a);
    t1 : Pieri_Tree(d,a-1) := Create(natural32(n),natural32(d),v1);
    a2 : constant integer32 := kp'last-a-1;
    v2 : constant Vector(0..a2-1) := kp(a+1..kp'last-1);
    t2 : Pieri_Tree(d,a2) := Create(natural32(n),natural32(d),v2);
    lp : List_of_Paired_Nodes := Create(n,d,t1,t2);
    np : constant Nodal_Pair(d) := Create(d,lp);
    sel_lp : List_of_Paired_Nodes;
    rc : constant natural32 := Length_Of(lp);
    nb : constant natural32 := Number_of_Paths(np);
    l1 : VecMat(0..a-1) := First_Random_Input_Sequence(n,m,a,kp);
    l2 : VecMat(0..kp'last-a-2) := Second_Random_Input_Sequence(n,m,a,kp);
    ln : Standard_Complex_Matrices.Matrix(1..n,1..m+1-integer32(kp(kp'last)))
       := Random_Matrix(natural32(n),
                        natural32(m)+1-natural32(kp(kp'last)));
    report,outlog : boolean;
    sols : VecMat(1..integer32(rc));

  begin
    skip_line;
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results.");
    new_line;
    put(file,"(m = "); put(file,m,1); put(file,",p = "); put(file,d,1);
    put(file,")");
    put(file," k = {");
    for i in 1..a-1 loop
      put(file,kp(i),1); put(file,",");
    end loop;
    put(file,kp(a),1);
    put(file,"}");
    put(file,"{");
    for i in a+1..kp'last-2 loop
      put(file,kp(i),1); put(file,",");
    end loop;
    put(file,kp(kp'last-1),1);
    put(file,"} ");
    put(file,kp(kp'last),1);
    new_line(file);
    put(file,"The first Pieri tree has height "); put(file,Height(t1),1);
    put_line(file," :"); Write_Tree(file,t1); -- put(file,t1);
    put(file,"The second Pieri tree has height "); put(file,Height(t2),1);
    put_line(file," :"); Write_Tree(file,t2); -- put(file,t2);
    put(file,"The root count equals : "); put(file,rc,1); new_line(file);
    put("The root count equals : "); put(rc,1); new_line;
    put(file,lp);
    put_line(file,"The tree of paired nodes :");
    Write(file,np);
    put("The number of paths : "); put(nb,1); new_line;
    put(file,"The number of paths : "); put(file,nb,1); new_line(file);
    new_line;
    put("Do you want real random input planes ? (y/n) "); get(ans);
    if ans = 'y'
     then Reallify(l1); Reallify(l2); Reallify(ln);
     else Orthogonalize(l1); Orthogonalize(l2);
          ln := Orthogonalize(ln);
    end if;
    put("Do you want moving cycles/polynomial systems on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    outlog := (ans = 'y');
    put("Do you want to select only certain pairs ? (y/n) ? ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then sel_lp := Select_Pairs(lp);
     else sel_lp := lp;
    end if;
    Set_Parameters(file,report);
    new_line(file);
    put_line(file,"The first sequence of random input planes :");
    put(file,l1,2);
    put_line(file,"The second sequence of random input planes :");
    put(file,l2,2);
    put_line(file,"The last random input plane :"); put(file,ln,2);
    put_line(file,"Starting the deformations at the paired leaves.");
    Matrix_Indeterminates.Initialize_Symbols(natural32(n),natural32(d));
    Add_t_Symbol;
    tstart(timer);
    Deform_Pairs(file,natural32(n),natural32(d),sel_lp,l1,l2,ln,report,
                 outlog,sols);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Pieri Deformations");
    Matrix_Indeterminates.Clear_Symbols;
    Clear(t1); Clear(t2); Clear(lp);
    Clear(l1); Clear(l2); Clear(sols);
    Close(file);
  end Root_Count;

-- FIVE MAJOR TEST PROGRAMS :

  procedure Test_Pieri_Condition ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Reads two brackets and tests Pieri's condition.

    b1,b2 : Bracket(1..d);
    ans : character;

  begin
    new_line;
    put_line("Interactive test on Pieri's condition for pair of brackets.");
    Matrix_Indeterminates.Initialize_Symbols(natural32(n),natural32(d));
    loop
      new_line;
      put("Give first bracket : "); get(b1);
      put("Give second bracket : "); get(b2);
      put("("); put(b1); put(","); put(b2); put(")");
      if Pieri_Condition(natural32(n),b1,b2)
       then put_line(" satisfies Pieri's condition.");
       else put_line(" does not satisfy Pieri's condition.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
    Matrix_Indeterminates.Clear_Symbols;
  end Test_Pieri_Condition;

  procedure Test_Pieri_Tree ( n,d : in integer32 ) is

  -- DESCRIPTION :  Constructs T(r_0,r_1,..,r_a).

    a : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Interactive test on the construction of Pieri trees");
    loop
      new_line;
      put("Give number of jumped-branching levels : "); get(a);
      declare
        v : Vector(0..a);
        t : Pieri_Tree(d,a);
      begin
        put("Give "); put(a+1,1); put(" natural numbers : "); get(v);
        t := Create(natural32(n),natural32(d),v);
        put_line("The Pieri tree : "); Write_Tree(t); --put(t);
        put_line("The chains in the Pieri tree :"); Write_Chains(t);
        Clear(t);
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Pieri_Tree;

  procedure Test_Root_Count ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Reads the input planes and sets up a pair of trees to perform
  --   the combinatorial root count.

    m : constant integer32 := n-d;
    sum : natural32;
    np,a : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Counting roots combinatorially by Pieri trees");
    loop
      new_line;
      put("Give the number of planes to intersect : "); get(np);
      declare 
        kp : Vector(1..np);
      begin
        put("Give "); put(np,1); put(" co-dimensions of the planes : ");
        get(kp);
        put("m = "); put(m,1); put(" p = "); put(d,1);
        put(" Sum : "); put(kp(kp'first),1);
        sum := kp(kp'first);
        for i in kp'first+1..kp'last loop
          put(" + "); put(kp(i),1);
          sum := sum + kp(i);
        end loop;
        put(" = "); put(sum,1);
        if integer32(sum) = m*d then
          put_line(" = m*p");
          loop
            put("Give number of elements in first set : "); get(a);
            Root_Count(n,d,a,kp);
            new_line;
            put("Do you want to test other partitions ? (y/n) "); get(ans);
            exit when ans /= 'y';
          end loop;
        else
          put_line(" /= m*p");
        end if;
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Root_Count;

  procedure Test_Pieri_Equations ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Set up of the expansions of the maximal minors.

    k,kd,nb : integer32 := 0;
    bm : Bracket_Monomial;
    ans : character;
    b1,b2 : Bracket(1..d);

  begin
    new_line;
    put_line("Set up equations for certain Schubert conditions.");
    Matrix_Indeterminates.Initialize_Symbols(natural32(n),natural32(d));
    loop
      new_line;
      put("Give k to determine (m-k+1)-plane : "); get(k);
      kd := n-k+1;
      bm := Maximal_Minors(natural32(n),natural32(kd));  -- because n = m+p
      put("All maximal minors : "); put(bm); new_line;
      declare
        bs : constant Bracket_System 
           := Minor_Equations(natural32(kd),natural32(kd-d),bm);
      begin
        put_line("The generic equation in the Laplace expansion : ");
        put(bs(0));
        put_line("The specific equations in the system : ");
        for i in 1..bs'last loop
          put(bs(i));
        end loop;
        Expand_Minors(n,d,bs);
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      Clear(bm);
      exit when (ans /= 'y');
    end loop;
    Matrix_Indeterminates.Clear_Symbols;
    new_line;
    put("Do you want to test memory consumption ? (y/n) "); get(ans);
    if ans = 'y' then
      put("Give number of times : "); get(nb);
      k := 1;
      for i in 1..d loop
        b1(i) := natural32(i);
        b2(i) := natural32(i);
      end loop;
      for i in 1..nb loop
        kd := n-k+1;
        bm := Maximal_Minors(natural32(n),natural32(kd));
        declare
          nva : constant integer32 := n*d+1;
          bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
             := Minor_Equations(natural32(kd),natural32(kd-d),bm);
          cffmat : constant Standard_Complex_Matrices.Matrix(1..n,1..(n-d))
                 := Random_Matrix(natural32(n),natural32(n-d));
          stamat : constant Standard_Complex_Matrices.Matrix(1..n,1..(n-d))
                 := Random_Matrix(natural32(n),natural32(n-d));
          movmat : Standard_Complex_Poly_Matrices.Matrix
                     (cffmat'range(1),cffmat'range(2))
                 := Moving_U_Matrix(nva,stamat,cffmat);
          polmat : Standard_Complex_Poly_Matrices.Matrix(1..n,1..d)
                 := Schubert_Pattern(natural32(n),b1,b2);  
          sys : Poly_Sys(1..bs'last)
              := Expanded_Minors(cffmat,polmat,bs);
          movcyc : Poly_Sys(1..bs'last)
                 := Lifted_Expanded_Minors(movmat,polmat,bs);
        begin
          Clear(bm); Clear(bs);
          Standard_Complex_Poly_Matrices.Clear(movmat);
          Standard_Complex_Poly_Matrices.Clear(polmat);
          Clear(sys); Clear(movcyc);
        end;
      end loop;
    end if;
  end Test_Pieri_Equations;

  procedure Test_Pieri_Homotopies ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Set up of the parametric families in the Pieri homotopy algorithm.

    b1,b2 : Bracket(1..d);
    m : constant integer32 := n-d;
    nva : constant integer32 := n*d+1;
    k,kd : integer32 := 0;
    bm : Bracket_Monomial;
    F1 : constant Standard_Complex_Matrices.Matrix
       := Random_Upper_Triangular(n);
    F2 : constant Standard_Complex_Matrices.Matrix
       := Random_Lower_Triangular(n);

  begin
    new_line;
    put_line("Set up the moving cycles in the Pieri homotopy algorithm.");
    new_line;
    Matrix_Indeterminates.Initialize_Symbols(natural32(n),natural32(d));
    Add_t_Symbol;
    put_line("The general upper triangular matrix F : "); put(F1,2);
    put_line("The general lower triangular matrix F' : "); put(F2,2);
    put("Give first "); put(d,1); put("-bracket : "); get(b1);
    put("Give second "); put(d,1); put("-bracket : "); get(b2);
    put("The pair of brackets : ");
    put("("); put(b1); put(","); put(b2); put_line(")");
    put("Give k to determine (m-k+1)-plane : "); get(k);
    kd := n-k+1;
    bm := Maximal_Minors(natural32(n),natural32(kd));
    declare
      L : constant Standard_Complex_Matrices.Matrix
        := Random_Matrix(natural32(n),natural32(m+1-k));
      X : Standard_Complex_Poly_Matrices.Matrix(1..n,b1'range);
      U1 : constant Standard_Complex_Matrices.Matrix := U_Matrix(F1,b1);
      movU1 : Standard_Complex_Poly_Matrices.Matrix(L'range(1),L'range(2));
      U2 : constant Standard_Complex_Matrices.Matrix := U_Matrix(F2,b2);
      movU2 : Standard_Complex_Poly_Matrices.Matrix(L'range(1),L'range(2));
      bs : constant Bracket_System
         := Minor_Equations(natural32(kd),natural32(kd-d),bm);
      movcyc1,movcyc2 : Poly_Sys(1..bs'last);
    begin
      put("The U-matrix for F and "); put(b1); put_line(" :"); put(U1,2);
      put("The U-matrix for F' and "); put(b2); put_line(" :"); put(U2,2);
      put("A random "); put(m+1-k,1); put_line("-plane :"); put(L,2);
      movU1 := Moving_U_Matrix(nva,U1,L);
      put_line("The first moving U-matrix :"); put(movU1);
      movU2 := Moving_U_Matrix(nva,U2,L);
      put_line("The second moving U-matrix :"); put(movU2);
      X := Schubert_Pattern(natural32(n),b1,b2);
      put("The unknown "); put(d,1); put_line("-plane X :"); put(X);
      movcyc1 := Lifted_Expanded_Minors(movU1,X,bs);
      put_line("The polynomial system for the first moving U-matrix :");
      put_line(movcyc1);
      movcyc2 := Lifted_Expanded_Minors(movU2,X,bs);
      put_line("The polynomial system for the second moving U-matrix :");
      put_line(movcyc2);
      put_line("Target system at t = 1 :");
      put_line(Eval(movcyc2,Create(1.0),
                    integer32(Number_of_Unknowns(movcyc2(1)))));
      put_line("System that must be solved :");
      put_line(Expanded_Minors(L,X,bs));
    end;
    Matrix_Indeterminates.Clear_Symbols;
  end Test_Pieri_Homotopies;

  procedure Main is

    m,p,n : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Interactive test for setting up of Pieri homotopies.");
    new_line;
    put("Give m, m+p = n, dimension of space : "); get(m);
    put("Give p, number of entries in brackets : "); get(p);
    n := m+p;
    new_line;
    put_line("MENU for testing the Pieri Homotopy Algorithm : ");
    put_line("  1. Construction of a Pieri tree. ");
    put_line("  2. Interactively test Pieri's condition on pairs of brackets.");
    put_line("  3. Set up equations for certain Schubert conditions.");
    put_line("  4. Set up the moving cycles in the Pieri homotopy algorithm.");
    put_line("  5. Test combinatorial root count and start deforming.");
    put("Make your choice (1, 2, 3, 4 or 5) : "); get(ans);
    case ans is
      when '1'    => Test_Pieri_Tree(n,p);
      when '2'    => Test_Pieri_Condition(n,p);
      when '3'    => Test_Pieri_Equations(n,p);
      when '4'    => Test_Pieri_Homotopies(n,p);
      when '5'    => Test_Root_Count(n,p);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_org_pieri;
