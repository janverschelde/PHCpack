-- NOTICE :
--   This is the Pieri solver which applies a random orthonormal transformation
--   on the input planes.

with text_io,integer_io;                 use text_io,integer_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_VecMats_io;        use Standard_Complex_VecMats_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Homotopy;
with Continuation_Parameters;
with Increment_and_Fix_Continuation;     use Increment_and_Fix_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Brackets;                           use Brackets;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;
with Curves_into_Grassmannian_io;        use Curves_into_Grassmannian_io;
with Localization_Posets;                use Localization_Posets;
with Deformation_Posets;                 use Deformation_Posets;
with Drivers_for_Input_Planes;           use Drivers_for_Input_Planes;
with Determinantal_Systems;              use Determinantal_Systems;
with Complex_Polynomial_Matrices;        use Complex_Polynomial_Matrices;
with Complex_Polynomial_Matrices_io;     use Complex_Polynomial_Matrices_io;
with Interfaces.C;
with C_Integer_Arrays,C_Double_Arrays;   use C_Integer_Arrays,C_Double_Arrays;
with C_Double_Arrays_io;                 use C_Double_Arrays_io;
with C_to_Ada_Arrays;                    use C_to_Ada_Arrays;

function Pieri_Solver ( m,p,q : integer;
                        points,planes : C_dblarrs.Pointer ) return integer is

  result : integer := 0;
  nq : constant natural := m*p + q*(m+p);
  dimpts : constant natural := 2*nq;
  valpts : C_Double_Array(0..Interfaces.C.size_T(dimpts-1))
         := C_dblarrs.Value(planes,Interfaces.C.ptrdiff_T(dimpts));
  dimpla : constant natural := 2*nq*m*p;
  valpla : C_Double_Array(0..Interfaces.C.size_T(dimpla-1))
         := C_dblarrs.Value(planes,Interfaces.C.ptrdiff_T(dimpla));

  function pieri_sols ( m,p,q,nbsols : integer;
                        ns : integer; s : C_Integer_Array;
                        nc : integer; c : C_Double_Array ) return integer;
  pragma Import(C,pieri_sols,"pieri_sols");

  -- DESCRIPTION :
  --   The function pieri_sols takes the solution maps as input from Ada to C.

  procedure Output_Results
               ( solmaps : in Array_of_Polynomial_Matrices;
                 orthmat : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Converts the solution maps into a basic format, ready for
  --   processing by a C function.

    d : constant Standard_Integer_Vectors.Vector := Degrees(solmaps);
    k : constant natural := d'length + Standard_Integer_Vectors.Sum(d);
    c : constant Standard_Complex_Vectors.Vector := Coefficients(k,solmaps);
    c_d : constant C_Integer_Array := Convert(d);
    c_c : constant C_Double_Array := Convert(c);
    c_o : constant C_Double_Array := Convert(orthmat);
    c_co : constant C_Double_Array := Concat(c_c,c_o);
    result_of_call : integer;

  begin
   -- put("The degrees : "); put(d); new_line;
   -- put_line("The coefficients : "); put_line(c);
    result_of_call := pieri_sols(m,p,q,solmaps'length,
                                 c_d'length,c_d,c_co'length,c_co);
  end Output_Results;

  procedure Refine_Roots ( file : in file_type;
                           p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the root refiner, for square polynomial systems only.

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural := 3;
    numit : natural := 0;

  begin
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,false,false);
  end Refine_Roots;

  procedure Refine_Roots ( p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the root refiner, for square polynomial systems only.

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural := 3;
    numit : natural := 0;

  begin
    Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max,false);
  end Refine_Roots;

  function Convert ( m,p,n : natural; cff : Standard_Complex_Vectors.Vector )
                   return VecMat is

  -- DESCRIPTION :
  --   The coefficient vector cff represents n complex p-by-m matrices.
  --   The vector on return contains n complex (m+p)-by-m matrices, where
  --   the identity matrix is added to the given matrices.

    res : VecMat(1..n);
    ind : integer := cff'first;

  begin
    for i in 1..n loop
      res(i) := new Standard_Complex_Matrices.Matrix(1..(m+p),1..m);
      for j1 in 1..p loop
        for j2 in 1..m loop
          res(i)(j1,j2) := cff(ind);
          ind := ind + 1;
        end loop;
      end loop;
      for j1 in 1..m loop
        for j2 in 1..m loop
          if j1 = j2
           then res(i)(p+j1,j2) := Create(1.0);
           else res(i)(p+j1,j2) := Create(0.0);
          end if;
        end loop;
      end loop;
    end loop;
    return res;
  end Convert;

  function Degree ( locmap : Standard_Natural_Matrices.Matrix;
                    i,j : natural ) return natural is

  -- DESCRIPTION :
  --   Returns the degree of the (i,j)-th entry in the polynomial matrix
  --   following the given localization map.  Note that the zero polynomial
  --   is represented by the constant zero.

    res : natural := 0;
    mp : constant natural := m+p;
    ind : natural;

  begin
    for k in 1..q loop
      ind := i + k*mp;
      exit when (ind > locmap'last(1));
      exit when (locmap(ind,j) = 0);
      res := res + 1;
    end loop;
    return res; 
  end Degree;

  function Solution_Map
              ( locmap : Standard_Natural_Matrices.Matrix;
                sol : Solution ) return Polynomial_Matrix is

  -- DESCRIPTION :
  --   Returns the polynomial matrix to represent a curve of degree q,
  --   following the localization map and using the coefficients in the
  --   solution.  Note that the coefficients have been stored columnwise.

    pm : Polynomial_Matrix(1..m+p,1..p);
    d,ind,row : natural;

  begin
   -- put_line("degrees : ");
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        d := Degree(locmap,i,j);
       -- put(" "); put(d,1);
        pm(i,j) := new Standard_Complex_Vectors.Vector'(0..d => Create(0.0));
        if d = 0
         then if locmap(i,j) = 0
               then pm(i,j)(0) := Create(0.0);
               else pm(i,j)(0) := Create(1.0);
              end if;
        end if;
      end loop;
     -- new_line;
    end loop;
    ind := 0;
    for j in locmap'range(2) loop
      d := 0; row := 0;
      for i in locmap'range(1) loop
        row := row + 1;
        if row > pm'last(1)
         then row := pm'first(1);
              d := d + 1;
        end if;
        if locmap(i,j) = 1
         then pm(row,j)(0) := Create(1.0);
         elsif locmap(i,j) = 2
             then ind := ind + 1;
                  pm(row,j)(d) := sol.v(ind);
        end if;
      end loop;
    end loop;
    return pm;
  end Solution_Map;
 
  procedure Create_Solution_Maps
              ( locmap : in Standard_Natural_Matrices.Matrix;
                sols : in Solution_List;
                orthmat : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Converts the list of solutions into maps of p-planes, 
  --   following the localization map.

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    solmaps : Array_of_Polynomial_Matrices(1..Length_Of(sols));

  begin
    for i in solmaps'range loop
      ls := Head_Of(tmp);
      solmaps(i) := new Polynomial_Matrix'(Solution_Map(locmap,ls.all));
      tmp := Tail_Of(tmp);
    end loop;
   -- put_line("The solution maps : "); put(solmaps);
    Output_Results(solmaps,orthmat);
  end Create_Solution_Maps;

  function Multiply ( mat : Standard_Complex_Matrices.Matrix;
                      vm : VecMat ) return VecMat is

  -- DESCRIPTION :
  --   The array of matrices on return contains the products of mat
  --   with all matrices in vm.

    res : VecMat(vm'range);

    use Standard_Complex_Matrices;

  begin
    for i in res'range loop
      res(i) := new Standard_Complex_Matrices.Matrix'(mat*vm(i).all);
    end loop;
    return res;
  end Multiply;

  procedure Solve_Target_System
              ( start_points : in Standard_Complex_Vectors.Vector;
                start_planes : in VecMat;
                index_poset : in Array_of_Array_of_Nodes;
                deform_poset : in Array_of_Array_of_VecMats ) is

  -- DESCRIPTION :
  --   Application of Cheater's homotopy to solve the given problem,
  --   starting at the solutions to the generic problem.

    orthmat : Standard_Complex_Matrices.Matrix(1..m+p,1..m+p)
            := Random_Orthogonal_Matrix(m+p,m+p);
    target_vals : constant Standard_Complex_Vectors.Vector := Convert(valpts);
    target_cffs : constant Standard_Complex_Vectors.Vector := Convert(valpla);
    converted_planes : constant VecMat := Convert(m,p,nq,target_cffs);
    target_planes : constant VecMat := Multiply(orthmat,converted_planes);
    solplanes : constant VecMat := deform_poset(nq)(1).all;
    top : constant Bracket := index_poset(nq)(1).top;
    bot : constant Bracket := index_poset(nq)(1).bottom;
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..m+p,1..p)
        := Symbolic_Create(m,p,q,top,bot);
    locmap : Standard_Natural_Matrices.Matrix(1..(m+p)*(q+1),1..p)
           := Standard_Coordinate_Frame(m,p,q,top,bot,solplanes(1).all);
    start : Poly_Sys(start_planes'range)
          := Create_Polynomial_System
                (top,bot,locmap,xpm,start_points,start_planes);
    target : Poly_Sys(target_planes'range)
           := Create_Polynomial_System
                (top,bot,locmap,xpm,target_vals,target_planes);
    sols : Solution_List := Solution_Planes(top,bot,locmap,solplanes);
    n : constant natural := target'last;
    a : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    b : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
   -- One_Set_up_Symbol_Table(m,p,q,top,bot);
   -- Reduce_Symbols(top,bot,locmap);
    Refine_Roots(Standard_Output,start,sols);
   -- Refine_Roots(start,sols);
   -- put_line("The localization map : "); put(locmap);
   -- put_line("The target interpolation points :");
   -- put_line(target_vals);
   -- put_line("The coefficients of the target planes :");
   -- put_line(target_cffs);
   -- put_line("The target planes : "); put(target_planes);
    Standard_Homotopy.Create(target,start,1,a,b,true);      -- linear cheater
    Set_Continuation_Parameter(sols,Create(0.0));
    Sil_Cont(sols,false,Create(1.0));
    Refine_Roots(Standard_Output,target,sols);
   -- Refine_Roots(target,sols);
    Create_Solution_Maps(locmap,sols,orthmat);
  end Solve_Target_System;

  procedure Solve_Deformation_Poset
              ( index_poset : in out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   This procedure calls the Pieri homotopies to solve a random
  --   complex instance, based on the poset.

    root : Node := index_poset(nq)(1).all;
    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    input : VecMat(1..nq) := Random_Complex_Planes(m,p,q);
    svals : Standard_Complex_Vectors.Vector(1..nq) := Random_Vector(1,nq);
    npaths : Standard_Natural_Vectors.Vector(1..nq) := (1..nq => 0);
    timings : Duration_Array(1..nq) := (1..nq => 0.0);

  begin
    Continuation_Parameters.Tune(0);
   -- Solve(Standard_Output,m+p,q,deform_poset,root,input,svals,true,false,
   --       npaths,timings);
    Solve(m+p,q,deform_poset,root,input,svals,npaths,timings);
    Solve_Target_System(svals,input,index_poset,deform_poset);
  end Solve_Deformation_Poset;

  procedure Combinatorial_Root_Count is

  -- DESCRIPTION :
  --   This procedure creates the posets for the combinatorial root count
  --   and calls the Pieri homotopies.  The only input parameters to this
  --   procedure is the triplet (m,p,q) which is known globally.

    root : Node(p) := Trivial_Root(m,p,q);
    lnkroot : Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..nq);
    index_poset : Array_of_Array_of_Nodes(0..nq);

  begin
    Q_Top_Bottom_Create(lnkroot,root.bottom(p),m+p);
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    Solve_Deformation_Poset(index_poset);
    result := level_poset(level_poset'last).roco;
  end Combinatorial_Root_Count;

begin
 -- put("The number of coefficients in the points : "); 
 -- put(dimpts,1); new_line;
 -- put_line("The coefficients of the interpolation points : ");
 -- put(dimpts,valpts);
 -- put("The number of coefficients in the planes : "); 
 -- put(dimpla,1); new_line;
 -- put_line("The coefficients of the input planes : ");
 -- put(dimpla,valpla);
  Combinatorial_Root_Count;
  return result;
end Pieri_Solver;
