with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Bracket_Monomials;                  use Bracket_Monomials;
with Bracket_Monomials_io;               use Bracket_Monomials_io;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Polynomials_io;    use Standard_Bracket_Polynomials_io;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;
with Pieri_Trees,Pieri_Trees_io;         use Pieri_Trees,Pieri_Trees_io;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Numeric_Minor_Equations;            use Numeric_Minor_Equations;
with Determinantal_Systems;              use Determinantal_Systems;
with Solve_Pieri_Leaves;
with Plane_Representations;              use Plane_Representations;
with Specialization_of_Planes;           use Specialization_of_Planes;
with Pieri_Continuation;                 use Pieri_Continuation;

package body Pieri_Deformations is

  function Coordinate_Bracket
             ( nd : Pieri_Node; jmp : integer32 ) return Bracket is

  -- DESCRIPTION :
  --   On return we find the bracket determines the local coordinates.
  --   This bracket depends whether jmp = nd.node'last or not.

  begin
    if Is_Leaf(nd) or jmp < nd.node'last
     then return nd.node;
     else return Upper_Jump_Decrease(nd);
    end if;
  end Coordinate_Bracket;

  procedure Test_Solution
              ( file : in file_type; nd1,nd2 : in Pieri_Node;
                id : in natural32; l1,l2 : in VecMat; l : in Matrix;
                x : Standard_Complex_Poly_Matrices.Matrix;
                solpla : in Matrix ) is

  -- DESCRIPTION :
  --   Evaluates the system that represents the intersection conditions
  --   for the planes in l1,l2 and the plane l in the given solution plane.

    sys : Link_to_Poly_Sys := new Poly_Sys'(Polynomial_Equations(l,x));
    eva : Standard_Complex_Vectors.Vector(sys'range);
    fst : integer32;

  begin
    if Is_Leaf(nd1) then
      fst := l1'last+1;
    elsif nd1.i = 0 then
      fst := integer32(nd1.c);
    else 
      fst := integer32(nd1.c)+1;
    end if;
    for i in fst..l1'last loop
      Concat(sys,Polynomial_Equations(l1(i).all,x));
    end loop;
    if Is_Leaf(nd2) then
      fst := l2'last+1;
    elsif nd2.i = 0 then
      fst := integer32(nd2.c);
    else 
      fst := integer32(nd2.c)+1;
    end if;
    for i in fst..l2'last loop
      Concat(sys,Polynomial_Equations(l2(i).all,x));
    end loop;
   -- put_line(file,"The polynomial system that has been solved :");
   -- put_line(file,sys.all);
    eva := Eval(sys.all,Vector_Rep(solpla));
    put(file,"Residual at the target system :");
    put(file,Max_Norm(eva),3);
    Clear(sys);
    if (((nd1.i = 0) and (nd1.c = 1))
      or else ((nd2.i = 0) and (nd2.c = 1)))
     then put(file," PAIR "); put(file,id,1); put_line(file," DONE.");
     else put_line(file,".");
    end if;
  end Test_Solution;

  procedure Test_Solutions
              ( file : in file_type;
                l1,l2 : in VecMat; l : in Matrix;
                 x : Standard_Complex_Poly_Matrices.Matrix;
                 solpla : in VecMat ) is

  -- DESCRIPTION :
  --   Evaluates the system that represents the intersection conditions
  --   for the planes in l1,l2 and the plane l in the given solution plane.

    sys : Link_to_Poly_Sys := new Poly_Sys'(Polynomial_Equations(l,x));

  begin
    for i in l1'range loop
      Concat(sys,Polynomial_Equations(l1(i).all,x));
    end loop;
    for i in l2'range loop
      Concat(sys,Polynomial_Equations(l2(i).all,x));
    end loop;
    put_line(file,"The polynomial system that has been solved :");
    put_line(file,sys.all);
    for i in solpla'range loop
      exit when (solpla(i) = null);
      put(file,"Solution plane no. "); put(file,i,1); put_line(file," :");
      put(file,solpla(i).all,2);
      put_line(file,"The system evaluated at the plane :");
      put_line(file,Eval(sys.all,Vector_Rep(solpla(i).all)));
    end loop;
    Clear(sys);
  end Test_Solutions;

  procedure Start_at_Leaves
              ( file : in file_type; pnd : in Paired_Nodes;
                ln : in Matrix; sol : in out Matrix ) is

  -- DESCRIPTION :
  --   Computes the start solutions at the leaves.

   -- xpm : Standard_Complex_Poly_Matrices.Matrix(sol'range(1),sol'range(2))
   --     := Schubert_Pattern(sol'last(1),pnd.left.node,pnd.right.node);
   -- sys : Link_to_Poly_Sys := new Poly_Sys'(Polynomial_Equations(ln,xpm));

  begin
    sol := Solve_Pieri_Leaves(file,pnd.left.node,pnd.right.node,ln);
    put_line(file,"The solution plane at the leaves : ");
    put(file,sol,2);
   -- put_line(file,"The polynomial equations : "); put(file,sys.all);
   -- put(file,"The solution evaluated at the polynomial equations : ");
   -- put(file,Evaluate(sys.all,sol),3); new_line(file);
   -- Standard_Complex_Poly_Matrices.Clear(xpm);
   -- Clear(sys);
  end Start_at_Leaves;

  function Moving_U_Matrix
             ( file : in file_type; outlog : in boolean;
               nd : Pieri_Node; lb : Bracket;
               f : Standard_Complex_Matrices.Matrix; l : VecMat )
             return Standard_Complex_Poly_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the moving U-matrix.

  -- REQUIRED : nd.c > 0.

    n : constant integer32 := f'length(1);
    p : constant integer32 := lb'length;
    nva : constant integer32 := n*p+1;
    m : constant integer32 := n-p;
    lc : constant integer32 := l(integer32(nd.c))'length(2);
    r : integer32;
    U : constant Standard_Complex_Matrices.Matrix := U_Matrix(f,lb);

  begin
    if outlog
     then put_line(file,"The U-matrix : "); put(file,U,2);
    end if;
    if nd.i = 0 then
      return Moving_U_Matrix(nva,U,l(integer32(nd.c)).all);
    else
      r := m+1 - lc;
      return Moving_U_Matrix(U,integer32(nd.i),r,lb);
    end if;
  end Moving_U_Matrix;

  function Moving_Cycle ( movU,x : Standard_Complex_Poly_Matrices.Matrix )
                        return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the moving cycle expaned as polynomial system.

    n : constant natural32 := natural32(x'length(1));
    p : constant natural32 := natural32(x'length(2));
    kd : constant natural32 := p + natural32(movU'length(2));
    bm : Bracket_Monomial := Maximal_Minors(n,kd);
    bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(kd,kd-p,bm);
    res : constant Poly_Sys := Lifted_Expanded_Minors(movU,x,bs);

  begin
    Clear(bm); Clear(bs);
    return res;
  end Moving_Cycle;

  function One_Moving_Cycle
              ( file : in file_type; nd : Pieri_Node; lb : Bracket; 
                f : Standard_Complex_Matrices.Matrix; l : VecMat;
                x : Standard_Complex_Poly_Matrices.Matrix; outlog : boolean )
              return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the moving part for one node in case i = 0.

  -- REQUIRED : nd.c > 0.

    lc : constant integer32 := l(integer32(nd.c))'length(2);
    movU : Standard_Complex_Poly_Matrices.Matrix(x'range(1),1..lc)
         := Moving_U_Matrix(file,outlog,nd,lb,f,l);
    res : constant Poly_Sys := Moving_Cycle(movU,x);

  begin
    if outlog then
      put_line(file,"The moving U matrix : "); put(file,MovU);
      put_line(file,"The equations of the moving cycle : "); put(file,res);
    end if;
    Standard_Complex_Poly_Matrices.Clear(movU);
    return res;
  end One_Moving_Cycle;

  function Left_Moving_Cycle
              ( file : file_type; nd : Pieri_Node; lb : Bracket;
                jmp : integer32;
                f : Standard_Complex_Matrices.Matrix; l : VecMat;
                x : Standard_Complex_Poly_Matrices.Matrix; outlog : boolean )
              return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the moving part for the case i > 0 and jmp < p,
  --   for the node from the first Pieri tree.

  -- REQUIRED : nd.c > 0.

    n : constant integer32 := x'length(1);
   -- p : constant natural32 := natural32(x'length(2));
    lc : constant integer32 := integer32(l(integer32(nd.c))'length(2));
    movU : Standard_Complex_Poly_Matrices.Matrix(x'range(1),1..lc)
         := Moving_U_Matrix(file,outlog,nd,lb,f,l);
    movUsec : constant Standard_Complex_Poly_Matrices.Matrix
            := Lower_Section(movU,n+1-integer32(lb(jmp)));
    res : constant Poly_Sys := Moving_Cycle(movUsec,x);

  begin
    if outlog then
      put_line(file,"The left moving U matrix : "); put(file,movU);
      put_line(file,"The left moving cycle : "); put(file,movUsec);
      put_line(file,"The equations of the moving cycle : "); put(file,res);
    end if;
    Standard_Complex_Poly_Matrices.Clear(movU);
    return res;
  end Left_Moving_Cycle;

  function Right_Moving_Cycle
              ( file : file_type; nd : Pieri_Node; lb : Bracket;
                jmp : integer32;
                f : Standard_Complex_Matrices.Matrix; l : VecMat;
                x : Standard_Complex_Poly_Matrices.Matrix; outlog : boolean )
              return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the moving part for the case i > 0 and jmp < p,
  --   for the node from the second Pieri tree.

  -- REQUIRED : nd.c > 0.

    lc : constant integer32 := l(integer32(nd.c))'length(2);
    movU : Standard_Complex_Poly_Matrices.Matrix(x'range(1),1..lc)
         := Moving_U_Matrix(file,outlog,nd,lb,f,l);
    movUsec : constant Standard_Complex_Poly_Matrices.Matrix
            := Upper_Section(movU,integer32(lb(jmp)));
    res : constant Poly_Sys := Moving_Cycle(movUsec,x);

  begin
    if outlog then
      put_line(file,"The right moving U matrix : "); put(file,movU);
      put_line(file,"The right moving cycle : "); put(file,movUsec);
      put_line(file,"The equations of the moving cycle :"); put(file,res);
    end if;
    Standard_Complex_Poly_Matrices.Clear(movU);
    return res;
  end Right_Moving_Cycle;

  procedure Moving_Cycles
              ( file : in file_type; pnd : in Paired_Nodes;
                id : in natural32;
                jmp1,jmp2 : in integer32; b1,b2 : in Bracket;
                f1,f2 : in Standard_Complex_Matrices.Matrix;
                l1,l2 : in VecMat; ln : in Matrix; report,outlog : in boolean;
                sol : in out Matrix ) is

  -- DESCRIPTION :
  --   Set up the moving cycles for the current pair of nodes down to
  --   the nodes (b1,b2), jumping along the indices (jmp1,jmp2).
  --   The other parameters are the same as in Deform_Pair.

    cb1 : constant Bracket := Coordinate_Bracket(pnd.left.all,jmp1);
    cb2 : constant Bracket := Coordinate_Bracket(pnd.right.all,jmp2);
    xpm : Standard_Complex_Poly_Matrices.Matrix(sol'range(1),sol'range(2))
        := Schubert_Pattern(natural32(sol'last(1)),cb1,cb2);
    homotopy : Link_to_Poly_Sys;
    locmap : constant Standard_Natural_Matrices.Matrix
           := Standard_Coordinate_Frame(xpm,sol);

  begin
    if outlog
     then
       put_line(file,"The localization map of the solution at the pair :");
       put(file,xpm);
    end if;
    homotopy := new Poly_Sys'(Polynomial_Equations(ln,xpm));
    for i in integer32(pnd.left.c)+1..l1'last loop
      Concat(homotopy,Polynomial_Equations(l1(i).all,xpm));
    end loop;
    for i in integer32(pnd.right.c)+1..l2'last loop
      Concat(homotopy,Polynomial_Equations(l2(i).all,xpm));
    end loop;
    for i in homotopy'range loop
      Embed(homotopy(i));
    end loop;
    if not Is_Leaf(pnd.left.all) and (pnd.left.c > 0) then
      if pnd.left.i = 0 then
        Concat(homotopy,
               One_Moving_Cycle(file,pnd.left.all,b1,f1,l1,xpm,outlog));
      elsif jmp1 < b1'last then
        Concat(homotopy,Left_Moving_Cycle
                          (file,pnd.left.all,b1,jmp1,f1,l1,xpm,outlog));
      end if;
    end if;
    if not Is_Leaf(pnd.right.all) and (pnd.right.c > 0) then
      if pnd.right.i = 0 then
        Concat(homotopy,
               One_Moving_Cycle(file,pnd.right.all,b2,f2,l2,xpm,outlog));
      elsif jmp2 < b2'last then
        Concat(homotopy,Right_Moving_Cycle
                          (file,pnd.right.all,b2,jmp2,f2,l2,xpm,outlog));
      end if;
    end if;
    if outlog then
      put_line(file,"All equations in the homotopy :");
      put_line(file,homotopy.all);
    end if;
    Trace_Paths(file,homotopy.all,locmap,report,outlog,sol);
    Test_Solution(file,pnd.left.all,pnd.right.all,id,l1,l2,ln,xpm,sol);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Clear(homotopy);
  end Moving_Cycles;

-- TARGET PROCEDURES :

  procedure Deform_Pair
              ( file : in file_type; pnd : in Paired_Nodes; id : in natural32;
                f1,f2 : in Standard_Complex_Matrices.Matrix;
                l1,l2 : in VecMat; ln : in Matrix; report,outlog : in boolean;
                sol : in out Matrix ) is

    down : Paired_Nodes;
    jmp1 : constant integer32 := Jump(pnd.left.all);
    jmp2 : constant integer32 := Jump(pnd.right.all);
    lb1 : constant Bracket := Lower_Jump_Decrease(pnd.left.all);
    lb2 : constant Bracket := Lower_Jump_Decrease(pnd.right.all);

  begin
    put(file,"JUMP @(");
    put(file,pnd.left.all); put(file,",");
    put(file,pnd.right.all); put(file,")"); put(file," = (");
    put(file,jmp1,1); put(file,","); put(file,jmp2,1); put(file,") -> (");
    put(file,lb1); put(file,","); put(file,lb2); put_line(file,").");
    if (Is_Leaf(pnd.left.all) and Is_Leaf(pnd.right.all))
     then Start_at_Leaves(file,pnd,ln,sol);
     else Moving_Cycles
            (file,pnd,id,jmp1,jmp2,lb1,lb2,f1,f2,l1,l2,ln,report,outlog,sol);
    end if;
    if not At_First_Branch_Point(pnd)
     then down := Ancestor(pnd); 
          Deform_Pair(file,down,id,f1,f2,l1,l2,ln,report,outlog,sol);
    end if;
   -- if pnd.left.h = pnd.right.h
   --  then if (((pnd.left.c <= 1) and (pnd.right.c <=1))
   --         and then (((pnd.left.i = 0) and (pnd.left.c = 1))
   --          or else ((pnd.right.i = 0) and (pnd.right.c = 1))))
   --        then null;
   --        else down.left := pnd.left.ancestor;
   --             down.right := pnd.right.ancestor;
   --             Deform_Pair(file,down,id,f1,f2,l1,l2,ln,report,outlog,sol);
   --       end if;
   --  elsif pnd.left.h > pnd.right.h
   --      then down.left := pnd.left.ancestor;
   --           down.right := pnd.right;
   --           Deform_Pair(file,down,id,f1,f2,l1,l2,ln,report,outlog,sol);
   --      else down.left := pnd.left;
   --           down.right := pnd.right.ancestor;
   --           Deform_Pair(file,down,id,f1,f2,l1,l2,ln,report,outlog,sol);
   -- end if;
  end Deform_Pair;

  procedure Write_Paired_Chains
              ( file : in file_type;
                pnd : in Paired_Nodes; ind : in integer32 ) is

  -- DESCRIPTION :
  --   Writes the chains downwards that start at the leaves of pnd.
  --   The paired nodes have index number ind.

  begin
    put(file,"DESCENDING THE PAIRED CHAINS "); put(file,ind,1);
    put_line(file," :");
    put(file,pnd.left);  new_line(file);
    put(file,pnd.right); new_line(file);
  end Write_Paired_Chains;

  procedure Deform_Pairs
              ( file : in file_type; n,d : in natural32;
                lp : in List_of_Paired_Nodes; l1,l2 : in VecMat;
                ln : in Matrix; report,outlog : in boolean;
                sols : out VecMat ) is

    tmp : List_of_Paired_Nodes := lp;
    pnd : Paired_Nodes;
    f1 : constant Standard_Complex_Matrices.Matrix
       := Random_Upper_Triangular(integer32(n));
    f2 : constant Standard_Complex_Matrices.Matrix
       := Random_Lower_Triangular(integer32(n));
    firstpair : constant Paired_Nodes := Head_Of(lp);
    lb1 : constant Bracket := Lowest_Jump_Decrease(firstpair.left.all);
    lb2 : constant Bracket := Lowest_Jump_Decrease(firstpair.right.all);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(ln'range(1),lb1'range)
        := Schubert_Pattern(natural32(ln'last(1)),lb1,lb2);

  begin
    for i in 1..integer32(Length_Of(lp)) loop
      pnd := Head_Of(tmp);
      Write_Paired_Chains(file,pnd,i);
      sols(i) := new Standard_Complex_Matrices.Matrix
                       (1..integer32(n),1..integer32(d));
      Deform_Pair(file,pnd,natural32(i),f1,f2,l1,l2,ln,report,outlog,
                  sols(i).all);
      tmp := Tail_Of(tmp);
    end loop;
    Test_Solutions(file,l1,l2,ln,xpm,sols);
  end Deform_Pairs;

end Pieri_Deformations;
