with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Complex_Norms_Equals;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Homotopy;
with Continuation_Parameters;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Brackets;                           use Brackets;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;
with Curves_into_Grassmannian_io;        use Curves_into_Grassmannian_io;
with Localization_Posets;                use Localization_Posets;
with Deformation_Posets;                 use Deformation_Posets;
with Make_Input_Planes;
with Determinantal_Systems;              use Determinantal_Systems;
with Complex_Polynomial_Matrices;        use Complex_Polynomial_Matrices;
with Complex_Polynomial_Matrices_io;     use Complex_Polynomial_Matrices_io;
with Verification_with_Determinants;     use Verification_with_Determinants;
with Verify_Solution_Maps;
with Interfaces.C;
with C_Integer_Arrays;                   use C_Integer_Arrays;
with C_to_Ada_Arrays;                    use C_to_Ada_Arrays;

function Pieri_Solver ( m,p,q,nb,output_level : integer32;
                        points,planes : C_dblarrs.Pointer;
                        filename : in chars_ptr ) return integer32 is

  result : integer32 := 0;
  timer : Timing_Widget;
  nq : constant integer32 := m*p + q*(m+p);
  wnb : integer32;
  dimpts : constant integer32 := 2*nq;
  valpts : constant C_Double_Array(0..Interfaces.C.size_T(dimpts-1))
         := C_dblarrs.Value(points,Interfaces.C.ptrdiff_T(dimpts));
  dimpla : constant integer32 := 2*nq*m*p;
  valpla : constant C_Double_Array(0..Interfaces.C.size_T(dimpla-1))
         := C_dblarrs.Value(planes,Interfaces.C.ptrdiff_T(dimpla));
  vs : constant string := Value(filename);

  function pieri_sols ( m,p,q,nbsols : integer32;
                        ns : integer32; s : C_Integer_Array;
                        nc : integer32; c : C_Double_Array;
                        npts : integer32; pts : C_Double_Array;
                        npla : integer32; pla : C_Double_Array;
                        name : chars_ptr ) return integer32;
  pragma Import(C,pieri_sols,"pieri_sols");

  -- DESCRIPTION :
  --   The function pieri_sols takes the solution maps as input from Ada to C.
  ---  It returns 0 if the function call ended well.

  procedure Output_Results ( solmaps : in Array_of_Polynomial_Matrices ) is

  -- DESCRIPTION :
  --   Converts the solution maps into a basic format, ready for
  --   processing by a C function.

    d : constant Standard_Integer_Vectors.Vector
      := Complex_Polynomial_Matrices.Degrees(solmaps);
    k : constant integer32 := d'length + Standard_Integer_Vectors.Sum(d);
    c : constant Standard_Complex_Vectors.Vector := Coefficients(k,solmaps);
    c_d : constant C_Integer_Array := Convert(d);
    c_c : constant C_Double_Array := Convert(c);
    result_of_call : integer32;

  begin
   -- put("The degrees : "); put(d); new_line;
   -- put_line("The coefficients : "); put_line(c);
    result_of_call := pieri_sols(m,p,q,solmaps'length,
                                 c_d'length,c_d,c_c'length,c_c,
                                 dimpts,valpts,dimpla,valpla,filename);
    if result_of_call /= 0
     then result := -1;
    end if;
  end Output_Results;

  procedure Refine_Roots ( file : in file_type;
                           p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the root refiner, for square polynomial systems only.

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,deflate,false);
  end Refine_Roots;

  procedure Refine_Roots ( p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the root refiner, for square polynomial systems only.

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max,deflate);
  end Refine_Roots;

  function Convert ( m,p,n : integer32; cff : Standard_Complex_Vectors.Vector )
                   return VecMat is

  -- DESCRIPTION :
  --   The coefficient vector cff represents n complex p-by-m matrices.
  --   The vector on return contains n complex (m+p)-by-m matrices, where
  --   the identity matrix is added to the given matrices.

    res : VecMat(1..n);
    ind : integer32 := cff'first;

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

  function Left_Multiply ( a : Standard_Floating_Matrices.Matrix;
                           b : VecMat ) return VecMat is

  -- DESCRIPTION :
  --   Multiplies every matrix in b with a from the left.

    res : VecMat(b'range);
    use Standard_Complex_Matrices;

  begin
    for i in res'range loop
      declare
        bi : constant Matrix(b(i)'range(1),b(i)'range(2)) := b(i).all;
        ab : Matrix(b(i)'range(1),b(i)'range(2));
      begin
        for j1 in ab'range(1) loop
          for j2 in ab'range(2) loop
            ab(j1,j2) := Create(0.0);
            for k in a'range(2) loop 
              ab(j1,j2) := ab(j1,j2) + Create(a(j1,k))*bi(k,j2);
            end loop;
          end loop;
        end loop;
        res(i) := new Standard_Complex_Matrices.Matrix'(ab);
      end;
    end loop;
    return res;
  end Left_Multiply;

--  function Left_Multiply ( a : Standard_Complex_Matrices.Matrix;
--                           b : VecMat ) return VecMat is
--
--  -- DESCRIPTION :
--  --   Multiplies every matrix in b with a from the left.
--
--    res : VecMat(b'range);
--    use Standard_Complex_Matrices;
--
--  begin
--    for i in res'range loop
--      declare
--        bi : Matrix(b(i)'range(1),b(i)'range(2)) := b(i).all;
--        ab : Matrix(b(i)'range(1),b(i)'range(2));
--      begin
--        for j1 in ab'range(1) loop
--          for j2 in ab'range(2) loop
--            ab(j1,j2) := Create(0.0);
--            for k in a'range(2) loop 
--              ab(j1,j2) := ab(j1,j2) + a(j1,k)*Conjugate(bi(k,j2));
--            end loop;
--          end loop;
--        end loop;
--        res(i) := new Standard_Complex_Matrices.Matrix'(ab);
--      end;
--    end loop;
--    return res;
--  end Left_Multiply;
 
  function Degree ( locmap : Standard_Natural_Matrices.Matrix;
                    i,j : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the degree of the (i,j)-th entry in the polynomial matrix
  --   following the given localization map.  Note that the zero polynomial
  --   is represented by the constant zero.

    res : integer32:= 0;
    mp : constant integer32 := m+p;
    ind : integer32;

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
    d,ind,row : integer32;

  begin
   -- put_line("degrees : ");
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        d := Degree(locmap,i,j);
       -- put(" "); put(d,1);
        pm(i,j) := new Standard_Complex_Vectors.Vector'(0..d => Create(0.0));
        if d = 0 then
          if locmap(i,j) = 0
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
         then row := pm'first(1); d := d + 1;
        end if;
        if locmap(i,j) = 1 then
          pm(row,j)(0) := Create(1.0);
        elsif locmap(i,j) = 2 then
          ind := ind + 1;
          pm(row,j)(d) := sol.v(ind);
        end if;
      end loop;
    end loop;
    return pm;
  end Solution_Map;
 
  function Create_Solution_Maps
             ( locmap : Standard_Natural_Matrices.Matrix;
               sols : Solution_List ) return Array_of_Polynomial_Matrices is

  -- DESCRIPTION :
  --   Converts the list of solutions into maps of p-planes, 
  --   following the localization map.

    res : Array_of_Polynomial_Matrices(1..integer32(Length_Of(sols)));
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
   -- put_line("The localization map : "); put(locmap);
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := new Polynomial_Matrix'(Solution_Map(locmap,ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create_Solution_Maps;

  procedure Solve_Target_System
              ( start_points : in Standard_Complex_Vectors.Vector;
                start_planes : in VecMat;
                index_poset : in Array_of_Array_of_Nodes;
                deform_poset : in Array_of_Array_of_VecMats ) is

  -- DESCRIPTION :
  --   Application of Cheater's homotopy to solve the given problem,
  --   starting at the solutions to the generic problem.

    target_vals : constant Standard_Complex_Vectors.Vector := Convert(valpts);
    target_cffs : constant Standard_Complex_Vectors.Vector := Convert(valpla);
    ct : constant Standard_Floating_Matrices.Matrix(1..m+p,1..m+p)
       := Random_Orthogonal_Matrix(natural32(m+p),natural32(m+p));
    ict : constant Standard_Floating_Matrices.Matrix(1..m+p,1..m+p)
        := Standard_Floating_Matrices.Transpose(ct);
    target_planes : constant VecMat := Convert(m,p,nq,target_cffs);
    ct_target_planes : constant VecMat := Left_Multiply(ct,target_planes);
    solplanes : constant VecMat := deform_poset(nq)(1).all;
    top : constant Bracket := index_poset(nq)(1).top;
    bot : constant Bracket := index_poset(nq)(1).bottom;
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..m+p,1..p)
        := Symbolic_Create(natural32(m),natural32(p),natural32(q),top,bot);
    locmap : constant Standard_Natural_Matrices.Matrix(1..(m+p)*(q+1),1..p)
           := Standard_Coordinate_Frame
                (natural32(m),natural32(p),natural32(q),
                 top,bot,solplanes(1).all);
    start : constant Poly_Sys(start_planes'range)
          := Create_Polynomial_System
                (top,bot,locmap,xpm,start_points,start_planes);
    target : constant Poly_Sys(ct_target_planes'range)
           := Create_Polynomial_System
                (top,bot,locmap,xpm,target_vals,ct_target_planes);
    sols : Solution_List
         := Solution_Planes(top,bot,locmap,solplanes(1..wnb));
    n : constant integer32 := target'last;
    a : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    b : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    solmaps,ict_solmaps : Array_of_Polynomial_Matrices(1..wnb);
    tol : constant double_float := 1.0E-8;
    fail : boolean;

    procedure Sil_Cont is
      new Silent_Continue(Standard_Complex_Norms_Equals.Max_Norm,
                          Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

    procedure Rep_Cont is
      new Reporting_Continue(Standard_Complex_Norms_Equals.Max_Norm,
                             Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    if output_level > 0 then
      One_Set_up_Symbol_Table(natural32(m),natural32(p),natural32(q),top,bot);
      Reduce_Symbols(top,bot,locmap);
     -- put_line("The matrix xpm : "); put(xpm);
      declare
        lxpm : constant Standard_Complex_Poly_Matrices.Matrix
                 (xpm'range(1),xpm'range(2))
             := Column_Localize(top,bot,locmap,xpm);
      begin
        put_line("The localized matrix xpm : "); put(lxpm);
      end;
    end if;
    if output_level > 1 then
      Refine_Roots(Standard_Output,start,sols);
      Verify_Determinants(Standard_Output,index_poset(nq)(1).all,
                          xpm,locmap,sols,start_points,start_planes);
      solmaps := Create_Solution_Maps(locmap,sols);
      put_line("Verifying the solution maps at the start...");
      Verify_Solution_Maps
        (Standard_Output,start_points,start_planes,solmaps,tol,fail);
    else
      Refine_Roots(start,sols);
    end if;  
   -- put_line("The localization map : "); put(locmap);
   -- put_line("The target interpolation points :");
   -- put_line(target_vals);
   -- put_line("The coefficients of the target planes :");
   -- put_line(target_cffs);
   -- put_line("The target planes : "); put(target_planes);
    Standard_Homotopy.Create(target,start,1,a,b,true);      -- linear cheater
    Set_Continuation_Parameter(sols,Create(0.0));
    if output_level > 2
     then Rep_Cont(Standard_Output,sols,false,target=>Create(1.0));
     else Sil_Cont(sols,false,target=>Create(1.0));
    end if;
    if output_level > 1 then
      Refine_Roots(Standard_Output,target,sols);
      Verify_Determinants(Standard_Output,index_poset(nq)(1).all,
                          xpm,locmap,sols,target_vals,ct_target_planes);
    elsif output_level > 0 then
      Refine_Roots(Standard_Output,target,sols);
    else
      Refine_Roots(target,sols);
    end if;
   -- declare
   --   id : Standard_Complex_Matrices.Matrix(ct'range(1),ct'range(2)) := ct;
   --   use Standard_Complex_Matrices;
   -- begin
   --   put_line("The identity matrix?");
   --   for i in id'range(1) loop
   --     for j in id'range(2) loop
   --       id(i,j) := Create(0.0);
   --       for k in ct'range(1) loop 
   --         id(i,j) := id(i,j) + ict(i,k)*Conjugate(ct(k,j));
   --       end loop;
   --     end loop;
   --   end loop;
   --   put(id,3);
   -- end;
    solmaps := Create_Solution_Maps(locmap,sols);
    put_line("The coordinate change matrix : "); put(ct,3);
    put_line("Transpose of the coordinate change matrix : "); put(ict,3);
    ict_solmaps := Left_Multiply(ict,solmaps);
    if output_level > 0 then
      put_line("The solution maps : "); put(solmaps);
      put_line("Verifying the solution maps at the target...");
      Verify_Solution_Maps
        (Standard_Output,target_vals,ct_target_planes,solmaps,tol,fail);
      put_line("The transformed solution maps : "); put(ict_solmaps);
      put_line("Verifying the transformed solution maps at the target...");
      Verify_Solution_Maps
        (Standard_Output,target_vals,target_planes,ict_solmaps,tol,fail);
    end if;
    Output_Results(ict_solmaps);
  end Solve_Target_System;

  procedure Solve_Deformation_Poset
              ( index_poset : in out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   This procedure calls the Pieri homotopies to solve a random
  --   complex instance, based on the poset.

    use Make_Input_Planes;

    root : constant Node := index_poset(nq)(1).all;
    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    input : constant VecMat(1..nq)
          := Random_Complex_Planes(natural32(m),natural32(p),natural32(q));
    svals : constant Standard_Complex_Vectors.Vector(1..nq)
          := Random_Vector(1,nq);
    npaths : Standard_Natural_Vectors.Vector(1..nq) := (1..nq => 0);
    timings : Duration_Array(1..nq) := (1..nq => 0.0);

  begin
    Continuation_Parameters.Tune(0);
    if output_level > 1
     then Solve(Standard_Output,natural32(m+p),natural32(q),natural32(wnb),
                deform_poset,root,input,svals,true,true,npaths,timings);
     else Solve(natural32(m+p),natural32(q),natural32(wnb),deform_poset,root,
                input,svals,npaths,timings);
    end if;
    Solve_Target_System(svals,input,index_poset,deform_poset);
  end Solve_Deformation_Poset;

  procedure Combinatorial_Root_Count is

  -- DESCRIPTION :
  --   This procedure creates the posets for the combinatorial root count
  --   and calls the Pieri homotopies.  The only input parameters to this
  --   procedure is the triplet (m,p,q) which is known globally.

    root : constant Node(p)
         := Trivial_Root(natural32(m),natural32(p),natural32(q));
    lnkroot : constant Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..nq);
    index_poset : Array_of_Array_of_Nodes(0..nq);

  begin
    Q_Bottom_Create(lnkroot,natural32(m+p));
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    if nb /= 0 then
      if nb < 0 or nb > 2**16
      -- the nb > 2**16 is a patch for conversion of a negative integer in C
      -- to a 32-bit integer in Ada
       then wnb := level_poset(level_poset'last).roco;
       else wnb := nb;
      end if;
      Solve_Deformation_Poset(index_poset);
    end if;
    if result /= -1
     then result := level_poset(level_poset'last).roco;
    end if;
  end Combinatorial_Root_Count;

  procedure Main is
  begin
    tstart(timer);
    put("The name of the file : "); put(vs); new_line;
   -- put("The number of coefficients in the points : "); 
   -- put(dimpts,1); new_line;
   -- put_line("The coefficients of the interpolation points : ");
   -- put(dimpts,valpts);
   -- put("The number of coefficients in the planes : "); 
   -- put(dimpla,1); new_line;
   -- put_line("The coefficients of the input planes : ");
   -- put(dimpla,valpla);
    Combinatorial_Root_Count;
    tstop(timer);
    if output_level > 0
     then print_times(Standard_Output,timer,"Pieri Solver");
    end if;
  end Main;

begin
  Main;
  return result;
end Pieri_Solver;
