with unchecked_deallocation;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Brackets_io;                        use Brackets_io;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Determinantal_Systems;              use Determinantal_Systems;
with Specialization_of_Planes;           use Specialization_of_Planes;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;
with Curves_into_Grassmannian_io;        use Curves_into_Grassmannian_io;
with Pieri_Homotopies;                   use Pieri_Homotopies;
with Pieri_Continuation;                 use Pieri_Continuation;

package body Deformation_Posets is

-- BRACKET AUXILIARITIES TO DETERMINE PIVOTS :

  function Complement ( n : natural32; b : Bracket ) return Bracket is

  -- DESCRIPTION :
  --   Returns the complement of the bracket b, defined as a bracket
  --   of range 1..n-b'length as an ordered subset of {1..n} \ b.

    res : Bracket(1..integer32(n)-b'last);
    cnt : integer32 := 0;
    ind : integer32 := 1;

  begin
    for i in 1..n loop
      if ((ind > b'last) or else (i < b(ind))) then
        cnt := cnt + 1;
        res(cnt) := i;
      elsif i = b(ind) then
        ind := ind + 1;
      end if;
    end loop;
    return res;
  end Complement;

  function Remove ( b : Bracket; ell : natural32 ) return Bracket is

  -- DESCRIPTION :
  --   Returns a smaller bracket that does not contain ell.

  -- REQUIRED : there exists a k: b(k) = ell.

    res : Bracket(1..b'last-1);
    cnt : integer32 := 0;

  begin
    for i in b'range loop
      if b(i) /= ell then
        cnt := cnt+1;
        res(cnt) := b(i);
      end if;
    end loop;
    return res;
  end Remove;

  function Is_In  ( b : Bracket; ell : natural32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there exists an index k such that b(k) = l.

  begin
    for k in b'range loop
      if b(k) = ell
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Remove ( b1,b2 : Bracket ) return Bracket is

  -- DESCRIPTION :
  --   Returns b1 minus the first element in b2 that also occurs in b1.

  begin
    for i in b2'range loop
      if Is_In(b1,b2(i))
       then return Remove(b1,b2(i));
      end if;
    end loop;
    return b1;
  end Remove;

  function Remove ( cols,b,subb : Bracket ) return Bracket is

  -- DESCRIPTION :
  --   The indices in cols correspond to the entries in b.
  --   The bracket subb is a sub-bracket of b, with only one entry removed.
  --   The indices on return correspond to the entries in subb.

    res : Bracket(subb'range);

  begin
    for i in subb'range loop
      if b(i) = subb(i)
       then res(i) := cols(i);
       else res(i) := cols(i+1);
      end if;
    end loop;
    return res;
  end Remove;

-- POSET-ORIENTED PIERI DEFORMATIONS :

  function Leaf_Plane ( n : natural32; nd : Node )
                      return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the solution plane that corresponds to a leaf of the poset.

    res : Standard_Complex_Matrices.Matrix(1..integer32(n),nd.top'range);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for i in nd.top'range loop
      res(integer32(nd.top(i)),i) := Create(1.0);
    end loop;
    return res;
  end Leaf_Plane;

-- DEFORMATIONS :

  function Path_Coordinates ( level,label,child : integer32 ) 
                            return string is
  begin
    return "tracing (User time) at node("
           & Convert(level) & ")("
           & Convert(label) & ") from child " 
           & Convert(child);
  end Path_Coordinates;

--  procedure Write_Path_Coordinates
--                 ( file : in file_type;
--                   level,label,path,child,childpath : in natural ) is
--
--  -- DESCRIPTION :
--  --   Writes all coordinates from the current path that is to be traced.
--
--  begin
--    put(file,"Tracing at node("); put(file,level,1); put(file,")(");
--    put(file,label,1); put(file,") path "); put(file,path,1);
--    put(file," as path "); put(file,childpath,1);
--    put(file," from child "); put(file,child,1); new_line(file);
--  end Write_Path_Coordinates;

  procedure Write_Path_Coordinates
                 ( file : in file_type;
                   level,label,child : in integer32 ) is

  -- DESCRIPTION :
  --   Writes all coordinates from the current path that is to be traced.

  begin
    put(file,"Tracing paths at node("); put(file,level,1); put(file,")(");
    put(file,label,1);
    put(file,") from child "); put(file,child,1); new_line(file);
  end Write_Path_Coordinates;

  procedure Deform_from_Children
                 ( poset : in out Array_of_Array_of_VecMats;
                   nd : in Node; n,uplevel : in integer32;
                   homotopy : in Poly_Sys;
                   x : in Standard_Complex_Poly_Matrices.Matrix;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Deforms from the i-th non-empty child of nd at uplevel.

  -- ON ENTRY :
  --   poset       poset of solution p-planes;
  --   nd          current node in the localization poset;
  --   n           dimension of the working space;
  --   uplevel     level where to find the start planes in the poset;
  --   homotopy    family of moving planes;
  --   x           matrix of unknowns according to a localization pattern.

  -- ON RETURN :
  --   poset       updated poset of solution p-planes;
  --   npaths      updated number of paths traced at that level;
  --   timings     updated elapsed user timings.

    locmap : Standard_Natural_Matrices.Matrix(1..n,1..nd.p);
    solcnt,label : integer32 := 0;

  begin
    for i in nd.child_labels'range loop
      label := integer32(nd.child_labels(i));
      if not Empty(poset,uplevel,label) then -- child.roco > 0
        declare
          planes : VecMat(poset(uplevel)(label)'range);
          timer : Timing_Widget;
        begin
          tstart(timer);
          for i in planes'range loop  -- create to avoid sharing
            planes(i) := new Standard_Complex_Matrices.Matrix'
                               (poset(uplevel)(label)(i).all);
          end loop;
          locmap := Standard_Coordinate_Frame(x,planes(planes'first).all);
          Trace_Paths(standard_output,homotopy,locmap,false,false,planes);
          for i in planes'range loop       
            solcnt := solcnt+1;
            poset(nd.level)(nd.label)(solcnt) := planes(i);
          end loop;
          tstop(timer);
          timings(nd.level) := timings(nd.level) + Elapsed_User_Time(timer);
        end;
      end if;
    end loop;
    npaths(nd.level) := npaths(nd.level) + natural32(solcnt);
  end Deform_from_Children;

  procedure Deform_from_Children
                 ( file : in file_type;
                   poset : in out Array_of_Array_of_VecMats;
                   nd : in Node; n,uplevel : in integer32;
                   homotopy : in Poly_Sys; report,outlog : in boolean;
                   x : in Standard_Complex_Poly_Matrices.Matrix;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Deforms from the i-th non-empty child of nd at uplevel.

  -- ON ENTRY :
  --   file        to write intermediate results on;
  --   poset       poset of solution p-planes;
  --   nd          current node in the localization poset;
  --   n           dimension of the working space;
  --   uplevel     level where to find the start planes in the poset;
  --   homotopy    family of moving planes;
  --   report      indicates whether intermediate output during continuation;
  --   outlog      flag to write homotopies on file if set to true. 
  --   x           matrix of unknowns according to a localization pattern.

  -- ON RETURN :
  --   poset       updated poset of solution p-planes;
  --   npaths      updated number of paths traced at that level;
  --   timings     updated elapsed user timings.

    locmap : Standard_Natural_Matrices.Matrix(1..n,1..nd.p);
    solcnt,label : integer32 := 0;

  begin
    for i in nd.child_labels'range loop
      label := integer32(nd.child_labels(i));
      if not Empty(poset,uplevel,label) then -- child.roco > 0
        Write_Path_Coordinates(file,nd.level,nd.label,label);
        declare
          planes : VecMat(poset(uplevel)(label)'range);
          timer : Timing_Widget;
        begin
          tstart(timer);
          for i in planes'range loop  -- create to avoid sharing
            planes(i) := new Standard_Complex_Matrices.Matrix'
                               (poset(uplevel)(label)(i).all);
          end loop;
          locmap := Standard_Coordinate_Frame(x,planes(planes'first).all);
          Trace_Paths(file,homotopy,locmap,report,outlog,planes);
          for i in planes'range loop       
            solcnt := solcnt+1;
            poset(nd.level)(nd.label)(solcnt) := planes(i);
          end loop;
          tstop(timer);
          new_line(file);
          print_times(file,timer,Path_Coordinates(nd.level,nd.label,label));
          new_line(file);
          timings(nd.level) := timings(nd.level) + Elapsed_User_Time(timer);
        end;
      end if;
    end loop;
    npaths(nd.level) := npaths(nd.level) + natural32(solcnt);
  end Deform_from_Children;

  procedure Quantum_Deform_from_Children
                 ( poset : in out Array_of_Array_of_VecMats; nd : in Node;
                   n,q,nb : in natural32; uplevel : in integer32;
                   cnt : in out natural32;
                   homotopy : in Poly_Sys; conpar,s_mode : in natural32;
                   x : in Standard_Complex_Poly_Matrices.Matrix;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Deforms from the i-th non-empty child of nd at uplevel.
  --   This is the quantum analogue to the hypersurface Pieri.

  -- ON ENTRY :
  --   poset       poset of solution p-planes;
  --   nd          current node in the localization poset;
  --   n           dimension of the working space;
  --   q           degree of the curve;
  --   nb          number of solution curves wanted;
  --   uplevel     level where to find the start planes in the poset;
  --   cnt         number of solution curves already computed;
  --   homotopy    family of moving planes;
  --   conpar      number of the continuation parameter;
  --   s_mode      if = 0, then s = 0, otherwise s = 1 at start;
  --   x           symbolic representation of the curve matrix of polynomials.

  -- ON RETURN :
  --   poset       updated poset of solution curves;
  --   cnt         updated counter of computed solution curves;
  --   npaths      updated number of paths at each level;
  --   timings     updated CPU user timings for each level.

    m : constant natural32 := n - natural32(nd.p);
    rws : constant integer32 := integer32(n*(q+1));
    locmap : Standard_Natural_Matrices.Matrix(1..rws,1..nd.p);
    solcnt,label : integer32 := 0;
    len : natural32;

  begin
    for i in nd.child_labels'range loop
      label := integer32(nd.child_labels(i));
      if not Empty(poset,uplevel,label) then    -- child.roco > 0 then
        len := poset(uplevel)(label)'length;
        if len > nb - cnt
         then len := nb - cnt;
        end if;
        declare
          ind : constant integer32 := poset(uplevel)(label)'first;
          planes : VecMat(ind..ind+integer32(len)-1);
          timer : Timing_Widget;
        begin
          tstart(timer);
          for i in planes'range loop  -- create to avoid sharing
            planes(i) := new Standard_Complex_Matrices.Matrix'
                               (poset(uplevel)(label)(i).all);
          end loop;
          locmap := Standard_Coordinate_Frame
                       (m,natural32(nd.p),q,nd.top,nd.bottom,
                        planes(planes'first).all);
          Quantum_Trace_Paths
            (m,natural32(nd.p),q,nd,homotopy,conpar,s_mode,locmap,planes);
          for i in planes'range loop       
            solcnt := solcnt+1;
            poset(nd.level)(nd.label)(solcnt) := planes(i);
          end loop;
          tstop(timer);
          timings(nd.level) := timings(nd.level) + Elapsed_User_Time(timer);
          cnt := cnt + len;
        end;
      end if;
      exit when (cnt = nb);
    end loop;
    npaths(nd.level) := npaths(nd.level) + natural32(solcnt);
  end Quantum_Deform_from_Children;

  procedure Quantum_Deform_from_Children
                 ( file : in file_type;
                   poset : in out Array_of_Array_of_VecMats; nd : in Node;
                   n,q,nb : in natural32; uplevel : in integer32;
                   cnt : in out natural32;
                   homotopy : in Poly_Sys; conpar,s_mode : in natural32;
                   report,outlog : in boolean;
                   x : in Standard_Complex_Poly_Matrices.Matrix;
                   s : in Standard_Complex_Vectors.Vector; ip : in VecMat;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Deforms from the i-th non-empty child of nd at uplevel.
  --   This is the quantum analogue to the hypersurface Pieri.

  -- ON ENTRY :
  --   file        to write intermediate results on;
  --   poset       poset of solution p-planes;
  --   nd          current node in the localization poset;
  --   n           dimension of the working space;
  --   q           degree of the curve;
  --   nb          number of solution maps wanted;
  --   uplevel     level where to find the start planes in the poset;
  --   cnt         counts the number of maps already computed;
  --   homotopy    family of moving planes;
  --   conpar      number of the continuation parameter;
  --   s_mode      if = 0, then s = 0, otherwise s = 1 at start;
  --   report      indicates whether intermediate output during continuation;
  --   outlog      flag to write homotopies on file if set to true;
  --   x           symbolic representation of the curve matrix of polynomials,
  --               not in localized form,
  --   s           interpolation points;
  --   ip          planes where the curve has to meet at the points.

  -- ON RETURN :
  --   poset       updated poset of solution maps;
  --   cnt         updated number of computed solution maps;
  --   npaths      updated number of paths at each level;
  --   timings     updated CPU user timings for each level.

    m : constant natural32 := n - natural32(nd.p);
    rws : constant natural32 := n*(q+1);
    locmap : Standard_Natural_Matrices.Matrix(1..integer32(rws),1..nd.p);
    solcnt,label : integer32 := 0;
    len : natural32;

  begin
    for i in nd.child_labels'range loop
      label := integer32(nd.child_labels(i));
      if not Empty(poset,uplevel,label) then -- child.roco > 0
        Write_Path_Coordinates(file,nd.level,nd.label,label);
        len := poset(uplevel)(label)'length;
        if len > nb - cnt
         then len := nb - cnt;
        end if;
        declare
          ind : constant integer32 := poset(uplevel)(label)'first;
          planes : VecMat(ind..ind+integer32(len)-1);
          timer : Timing_Widget;
        begin
          tstart(timer);
          for i in planes'range loop  -- create to avoid sharing
            planes(i) := new Standard_Complex_Matrices.Matrix'
                               (poset(uplevel)(label)(i).all);
          end loop;
          locmap := Standard_Coordinate_Frame
                      (m,natural32(nd.p),q,nd.top,nd.bottom,
                       planes(planes'first).all);
          Quantum_Trace_Paths
            (file,m,natural32(nd.p),q,nd,x,s,ip,homotopy,conpar,s_mode,locmap,
             report,outlog,planes);
          for i in planes'range loop       
            solcnt := solcnt+1;
            poset(nd.level)(nd.label)(solcnt) := planes(i);
          end loop;
          tstop(timer);
          new_line(file);
          print_times(file,timer,Path_Coordinates(nd.level,nd.label,label));
          new_line(file);
          timings(nd.level) := timings(nd.level) + Elapsed_User_Time(timer);
          cnt := cnt + len;
        end;
      end if;
      exit when (cnt = nb);
    end loop;
    npaths(nd.level) := npaths(nd.level) + natural32(solcnt);
  end Quantum_Deform_from_Children;

  procedure Hypersurface_Deform
                 ( n : in natural32;
                   poset : in out Array_of_Array_of_VecMats;
                   nd : in Node; expbp : in Bracket_Polynomial;
                   planes : in VecMat;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Does the Pieri deformations to the node, same specifications as Solve
  --   for the hypersurface case.

  -- REQUIRED : nd.level > 0, solutions at children of nd are in poset.

    xpm : Standard_Complex_Poly_Matrices.Matrix
            (1..integer32(n),1..nd.p)
        := Localization_Pattern(n,nd.top,nd.bottom);
    homsys : Poly_Sys(1..nd.level);

  begin
    if nd.tp = mixed then
      homsys := Two_Hypersurface_Pieri_Homotopy
                  (integer32(n),nd,expbp,xpm,planes);
      Deform_from_Children(poset,nd,integer32(n),nd.level-2,
                           homsys,xpm,npaths,timings);
    else
      homsys := One_Hypersurface_Pieri_Homotopy
                  (integer32(n),nd,expbp,xpm,planes);
      Deform_from_Children(poset,nd,integer32(n),nd.level-1,
                           homsys,xpm,npaths,timings);
    end if;
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Clear(homsys);
  end Hypersurface_Deform;

  procedure Hypersurface_Deform
                 ( file : in file_type; n : in natural32;
                   poset : in out Array_of_Array_of_VecMats;
                   nd : in Node; expbp : in Bracket_Polynomial;
                   planes : in VecMat; report,outlog : in boolean;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Does the Pieri deformations to the node, same specifications as Solve
  --   for the hypersurface case.

  -- REQUIRED : nd.level > 0, solutions at children of nd are in poset.

    xpm : Standard_Complex_Poly_Matrices.Matrix(1..integer32(n),1..nd.p)
        := Localization_Pattern(n,nd.top,nd.bottom);
    homsys : Poly_Sys(1..nd.level);

  begin
    if nd.tp = mixed then
      homsys := Two_Hypersurface_Pieri_Homotopy
                  (integer32(n),nd,expbp,xpm,planes);
      Deform_from_Children
        (file,poset,nd,integer32(n),nd.level-2,homsys,report,outlog,
         xpm,npaths,timings);
    else
      homsys := One_Hypersurface_Pieri_Homotopy
                  (integer32(n),nd,expbp,xpm,planes);
      Deform_from_Children
        (file,poset,nd,integer32(n),nd.level-1,homsys,report,outlog,
         xpm,npaths,timings);
    end if;
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Clear(homsys);
  end Hypersurface_Deform;

  procedure One_General_Deform
                  ( file : in file_type; n,ind : in natural32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    start,target : in Standard_Complex_Matrices.Matrix;
                    planes : in VecMat; bs : in Bracket_System;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Does the Pieri deformations to the node, same specifications as Solve
  --   for the general case.

  -- REQUIRED : nd.level > 0, solutions at children of nd are in poset.

  -- ON ENTRY :
  --   file         to write intermediate output to;
  --   n            number of rows in the matrices;
  --   ind          planes(ind) is currently being folded in with this chain;
  --   poset        contains solution planes at higher levels;
  --   nd           current node in the localization poset;
  --   start        start (m+1-k)-plane for pivots;
  --   target       target (m+1-k)-plane for pivots;
  --   planes       target planes;
  --   bs           structure to expand the minors;
  --   report       switch to determine output during continuation;
  --   outlog       flag to write homotopies on file if set to true. 

  -- ON RETURN :
  --   poset        solution planes at (nd.level)(nd.label) are determined;
  --   npaths       number of paths followed at each level;
  --   timings      CPU user time at each level.

    xpm : Standard_Complex_Poly_Matrices.Matrix(1..integer32(n),1..nd.p)
        := Localization_Pattern(n,nd.top,nd.bottom);
    hom : Link_to_Poly_Sys
        := One_General_Pieri_Homotopy
             (integer32(n),integer32(ind),nd,bs,start,target,xpm,planes);

  begin
    Deform_from_Children
      (file,poset,nd,integer32(n),nd.level-1,hom.all,report,outlog,
       xpm,npaths,timings);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Clear(hom);
  end One_General_Deform;

  procedure Two_General_Deform
                  ( file : in file_type; n,ind : in integer32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    top_start,top_target,bot_start,bot_target
                      : in Standard_Complex_Matrices.Matrix;
                    planes : in VecMat; top_bs,bot_bs : in Bracket_System;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Does the Pieri deformations to the node, same specifications as Solve
  --   for the general case.

  -- REQUIRED : nd.level > 0, solutions at children of nd are in poset.

  -- ON ENTRY :
  --   file         to write intermediate output to;
  --   n            number of rows in the matrices;
  --   ind          planes(ind) is currently being folded in with this chain;
  --   poset        contains solution planes at higher levels;
  --   nd           current node in the localization poset;
  --   top_start    start (m+1-k)-plane for top pivots;
  --   top_target   target (m+1-k)-plane for top pivots;
  --   bot_start    start (m+1-k)-plane for bottom pivots;
  --   bot_target   target (m+1-k)-plane for bottom pivots;
  --   planes       target planes;
  --   top_bs       structure to expand the minors for top pivots;
  --   bot_bs       structure to expand the minors for bottom pivots;
  --   report       switch to determine output during continuation;
  --   outlog       flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   poset        solution planes at (nd.level)(nd.label) are determined;
  --   npaths       number of paths traced at each level;
  --   timings      updated CPU user times for each level.

    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..nd.p)
        := Localization_Pattern(natural32(n),nd.top,nd.bottom);
    homotopy : Link_to_Poly_Sys;

  begin
    case nd.tp is
      when top
        => homotopy := One_General_Pieri_Homotopy
                         (n,ind,nd,top_bs,top_start,top_target,xpm,planes);
           Deform_from_Children
             (file,poset,nd,n,nd.level-1,homotopy.all,
              report,outlog,xpm,npaths,timings);
      when bottom
        => homotopy := One_General_Pieri_Homotopy
                         (n,ind,nd,bot_bs,bot_start,bot_target,xpm,planes);
           Deform_from_Children
             (file,poset,nd,n,nd.level-1,homotopy.all,
              report,outlog,xpm,npaths,timings);
      when mixed
        => homotopy := Two_General_Pieri_Homotopy
                         (n,ind,nd,top_bs,bot_bs,top_start,top_target,
                          bot_start,bot_target,xpm,planes);
           Deform_from_Children
             (file,poset,nd,n,nd.level-2,homotopy.all,
              report,outlog,xpm,npaths,timings);
    end case;
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Clear(homotopy);
  end Two_General_Deform;

  procedure Quantum_Deform
             ( n,q,nb : in natural32; cnt : in out natural32;
               poset : in out Array_of_Array_of_VecMats;
               nd : in Node; expbp : in Bracket_Polynomial;
               planes : in VecMat; s : Standard_Complex_Vectors.Vector;
               npaths : in out Standard_Natural_Vectors.Vector;
               timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   This is the q-analogue to the Hypersurface Deform.

  -- REQUIRED : nd.level > 0, solutions at children of nd are in poset.

    m : constant natural32 := n - natural32(nd.p);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..integer32(n),1..nd.p)
        := Symbolic_Create(m,natural32(nd.p),q,nd.top,nd.bottom);

  begin
    if nd.tp = mixed then
      declare
        homsys : Poly_Sys(1..nd.level+2);
      begin
        homsys := Two_Quantum_Pieri_Homotopy
                    (integer32(n),nd,expbp,xpm,planes,s);
        Quantum_Deform_from_Children
          (poset,nd,n,q,nb,nd.level-2,cnt,
           homsys,natural32(nd.level+3),1,xpm,npaths,timings);
        Clear(homsys);
      end;
    else
      declare
        homsys : Poly_Sys(1..nd.level+1);
      begin
        homsys := One_Quantum_Pieri_Homotopy
                    (integer32(n),nd,expbp,xpm,planes,s);
        Quantum_Deform_from_Children
          (poset,nd,n,q,nb,nd.level-1,cnt,
           homsys,natural32(nd.level+2),1,xpm,npaths,timings);
        Clear(homsys);
      end;
    end if;
    Standard_Complex_Poly_Matrices.Clear(xpm);
  end Quantum_Deform;

  procedure Quantum_Deform
             ( file : in file_type;
               n,q,nb : in natural32; cnt : in out natural32;
               poset : in out Array_of_Array_of_VecMats;
               nd : in Node; expbp : in Bracket_Polynomial;
               planes : in VecMat; s : Standard_Complex_Vectors.Vector;
               report,outlog : in boolean;
               npaths : in out Standard_Natural_Vectors.Vector;
               timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   This is the q-analogue to the Hypersurface Deform.

  -- REQUIRED : nd.level > 0, solutions at children of nd are in poset.

    m : constant natural32 := n - natural32(nd.p);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..integer32(n),1..nd.p)
        := Symbolic_Create(m,natural32(nd.p),q,nd.top,nd.bottom);

  begin
    if outlog then
      put(file,"Curve at node(");
      put(file,nd.level,1); put(file,")("); put(file,nd.label,1);
      put(file,") for pivots (");
      put(file,nd.top); put(file,","); put(file,nd.bottom);
      put_line(file,") :");
      One_Set_up_Symbol_Table(m,natural32(nd.p),q,nd.top,nd.bottom);
      put(file,xpm);
    end if;
    if nd.tp = mixed then
      declare
        homsys : Poly_Sys(1..nd.level+2);
      begin
        homsys := Two_Quantum_Pieri_Homotopy
                    (integer32(n),nd,expbp,xpm,planes,s);
        if outlog then
          Two_Set_up_Symbol_Table(m,natural32(nd.p),q,nd.top,nd.bottom);
          put_line(file,"The homotopy : "); put_line(file,homsys);
        end if;
        Quantum_Deform_from_Children
          (file,poset,nd,n,q,nb,nd.level-2,cnt,homsys,
           natural32(nd.level+3),1,report,outlog,xpm,s,planes,npaths,timings);
        Clear(homsys);
      end;
    else
      declare
        homsys : Poly_Sys(1..nd.level+1);
      begin
        homsys := One_Quantum_Pieri_Homotopy
                    (integer32(n),nd,expbp,xpm,planes,s);
        if outlog
         then put_line(file,"The homotopy : "); put_line(file,homsys);
        end if;
        Quantum_Deform_from_Children
          (file,poset,nd,n,q,nb,nd.level-1,cnt,homsys,
           natural32(nd.level+2),1,report,outlog,xpm,s,planes,npaths,timings);
        Clear(homsys);
      end;
    end if;
    Standard_Complex_Poly_Matrices.Clear(xpm);
  end Quantum_Deform;

  function Moving_Point_Mode
             ( l,k : natural32; modpiv : Bracket;
               start : Standard_Complex_Matrices.Matrix ) return natural32 is

  -- DESCRIPTION : 
  --   Returns a natural number that indicates the moving of the 
  --   interpolation point.  The value on return means the following
  --   when = 0 : s goes from 0 to 1;
  --        = 1 : s remains constant at 1;
  --        = 2 : s goes from 1 to a target value.

    tol : constant double_float := (10.0**(-12));

  begin
    if l = 0 then
      return 2;
    elsif l = k-1 and modpiv(1) > 1 then
      if AbsVal(start(2,1)) > tol
       then return 0;
       else return 1;
      end if;
      else return 1;
    end if;
  end Moving_Point_Mode;

  procedure One_General_Quantum_Deform
                  ( file : in file_type;
                    n,q,nb,l,k,ind : in natural32; cnt : in out natural32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    modpiv : in Bracket;
                    start,target : in Standard_Complex_Matrices.Matrix;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    bs : in Bracket_System; report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   This is the quantum analogue to One_General_Deform.

  -- REQUIRED : nd.level > 0, solutions at children of nd are in poset.

  -- ON ENTRY :
  --   file         to write intermediate output to;
  --   n            number of rows in the matrices;
  --   q            degree of the map;
  --   nb           number of solution maps wanted;
  --   l            runs from k-1 downto 0;
  --   k            co-dimension currently being treated;
  --   ind          planes(ind) is currently being folded in with this chain;
  --   cnt          counts the number of already computed maps;
  --   poset        contains solution planes at higher levels;
  --   nd           current node in the localization poset;
  --   modpiv       bottom or top pivots modulo n;
  --   start        start (m+1-k)-plane for pivots;
  --   target       target (m+1-k)-plane for pivots;
  --   planes       target planes;
  --   s            interpolation points where the maps meets the planes;
  --   bs           structure to expand the minors;
  --   report       switch to determine output during continuation;
  --   outlog       flag to write homotopies on file if set to true. 

  -- ON RETURN :
  --   cnt          updated counter of solution maps;
  --   poset        solution planes at (nd.level)(nd.label) are determined;
  --   npaths       number of paths followed at each level;
  --   timings      CPU user time at each level.

    m : constant natural32 := n - natural32(nd.p);
    xpm : Standard_Complex_Poly_Matrices.Matrix
            (1..integer32(n),1..nd.p)
        := Symbolic_Create(m,natural32(nd.p),q,nd.top,nd.bottom);
    s_mode : constant natural32 := Moving_Point_Mode(l,k,modpiv,start);
    hom : Link_to_Poly_Sys
        := One_General_Quantum_Pieri_Homotopy
             (integer32(n),integer32(ind),nd,s_mode,bs,
              start,target,xpm,planes,s);

  begin
    if outlog then
      put(file,"level l : "); put(file,l,1); put(file,"  ");
      put(file,"codim k : "); put(file,k,1); put(file,"  ");
      put(file,"modpiv : "); put(file,modpiv); put(file,"  ");
      put(file,"s_mode : "); put(file,s_mode,1); new_line(file);
      put(file,"Curve at node(");
      put(file,nd.level,1); put(file,")("); put(file,nd.label,1);
      put(file,") for pivots (");
      put(file,nd.top); put(file,","); put(file,nd.bottom);
      put_line(file,") :");
      One_Set_up_Symbol_Table(m,natural32(nd.p),q,nd.top,nd.bottom);
      put(file,xpm);
      put_line(file,"The homotopy : "); put_line(file,hom.all);
    end if;
    Quantum_Deform_from_Children
      (file,poset,nd,n,q,nb,nd.level-1,cnt,hom.all,
       natural32(nd.level+2),s_mode,report,outlog,xpm,s,planes,npaths,timings);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Clear(hom);
  end One_General_Quantum_Deform;

-- CREATORS :

  function Create ( index_poset : Array_of_Array_of_Nodes )
                  return Array_of_Array_of_VecMats is

    res : Array_of_Array_of_VecMats(index_poset'range);
    lnd : Link_to_Node;

  begin
    for i in index_poset'range loop
      if index_poset(i) /= null then
        res(i) := new Array_of_VecMats(index_poset(i)'range);
        for j in res(i)'range loop
          lnd := index_poset(i)(j);
          if lnd.roco /= 0
           then res(i)(j) := new VecMat(1..lnd.roco);
          end if;
        end loop;
      end if;
    end loop;
    return res;
  end Create;

-- SELECTORS :

  function Empty ( poset : Array_of_Array_of_VecMats;
                   level,label : integer32 ) return boolean is

    use Standard_Complex_Matrices;

  begin
    if poset(level) = null then
      return true;
    elsif poset(level)(label) = null then
      return true;
    else
      declare
        lavm : constant Link_to_VecMat := poset(level)(label);
      begin
        return (lavm(lavm'first) = null);
      end;
    end if;
  end Empty;

-- ANALOGUES TO THE ROOT COUNTERS :

  procedure Recursive_Hypersurface_Solve
               ( n : in natural32;
                 nd : in Node; expbp : in Bracket_Polynomial;
                 poset : in out Array_of_Array_of_VecMats;
                 planes : in VecMat;
                 npaths : in out Standard_Natural_Vectors.Vector;
                 timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   This additional layer is added to avoid the repeated construction
  --   of the structure of the equations, that is now in expbp.

  begin
    if Empty(poset,nd.level,nd.label) then
      if nd.level = 0 then
        poset(nd.level)(nd.label)(1)
          := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n,nd));
      else
        for i in nd.children'range(1) loop
          for j in nd.children'range(2) loop
            if nd.children(i,j) /= null then
              Recursive_Hypersurface_Solve
                (n,nd.children(i,j).all,expbp,poset,planes,npaths,timings);
            end if;
          end loop;
        end loop;
        Hypersurface_Deform(n,poset,nd,expbp,planes,npaths,timings);
       end if;
    end if;
  end Recursive_Hypersurface_Solve;

  procedure Recursive_Hypersurface_Solve
               ( file : in file_type; n : in natural32;
                 nd : in Node; expbp : in Bracket_Polynomial;
                 poset : in out Array_of_Array_of_VecMats;
                 planes : in VecMat; report,outlog : in boolean;
                 npaths : in out Standard_Natural_Vectors.Vector;
                 timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   This additional layer is added to avoid the repeated construction
  --   of the structure of the equations, that is now in expbp.

  begin
    if Empty(poset,nd.level,nd.label) then
      if nd.level = 0 then
        poset(nd.level)(nd.label)(1)
          := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n,nd));
      else
        for i in nd.children'range(1) loop
          for j in nd.children'range(2) loop
            if nd.children(i,j) /= null then
              Recursive_Hypersurface_Solve
                (file,n,nd.children(i,j).all,expbp,
                 poset,planes,report,outlog,npaths,timings);
            end if;
          end loop;
        end loop;
        Hypersurface_Deform
          (file,n,poset,nd,expbp,planes,report,outlog,npaths,timings);
       end if;
    end if;
  end Recursive_Hypersurface_Solve;

  procedure Solve ( n : in natural32;
                    poset : in out Array_of_Array_of_VecMats;
                    nd : in Node; planes : in VecMat;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

    bm : Bracket_Monomial := Maximal_Minors(n,n);
    bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(n,n-natural32(nd.p),bm);

  begin
    Recursive_Hypersurface_Solve(n,nd,bs(1),poset,planes,npaths,timings);
    Clear(bm); Clear(bs);
  end Solve;

  procedure Solve ( file : in file_type; n : in natural32;
                    poset : in out Array_of_Array_of_VecMats;
                    nd : in Node; planes : in VecMat;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

    bm : Bracket_Monomial := Maximal_Minors(n,n);
    bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(n,n-natural32(nd.p),bm);

  begin
    Recursive_Hypersurface_Solve
      (file,n,nd,bs(1),poset,planes,report,outlog,npaths,timings);
    Clear(bm); Clear(bs);
  end Solve;

  procedure One_Solve_along_Chains
                 ( file : in file_type; nd : in Node; n,l,k,ind : in natural32;
                   poset : in out Array_of_Array_of_VecMats;
                   pivots,columns : in Bracket; bs : in Bracket_System;
                   special,start,target : in Standard_Complex_Matrices.Matrix;
                   planes : in VecMat; report,outlog : in boolean;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Applies the general solver along the nodes in all chains that end at
  --   the current node.  This is the analogue to the hypersurface solver,
  --   for use in connection with the grandchildren first recursive solving.
  --   This procedure is only called in "One_Solve".

  -- ON ENTRY :
  --   file        to write intermediate results on;
  --   nd          current node in the localization poset;
  --   n           working dimension, equation m+p;
  --   l           runs from 0 to k-1;
  --   k           current codimension condition;
  --   poset       structure with all solution p-planes;
  --   ind         ind-1 planes are already folded in;
  --   pivots      pivot elements used for the special m-plane;
  --   columns     which columns of the special m-plane are used;
  --   bs          Laplace expansion of the polynomial equations;
  --   special     special m-plane for top pivots;
  --   start       (m+1-k)-plane used at the start of the deformation;
  --   target      (m+1-k)-plane used as target;
  --   planes      sequence of (m+1-k)-planes;
  --   report      indicates whether intermediate output during continuation;
  --   outlog      flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   poset       updated structure of all solution p-planes;
  --   npaths      updated numbers of paths traced at each level;
  --   timings     updated CPU user times at each level.

    m : constant natural32 := n - natural32(nd.p);
    new_piv,new_col : Bracket(1..pivots'last-1);
    new_start : Standard_Complex_Matrices.Matrix
                  (1..integer32(n),start'range(2));

  begin
    if empty(poset,nd.level,nd.label) then
      if l < k-1 then
        for i in nd.children'range(1) loop
          for j in nd.children'range(2) loop
            if ((nd.children(i,j) /= null)
              and then (nd.children(i,j).roco > 0))
             then
               if nd.children(i,j).tp = top
                then new_piv := Remove(pivots,nd.children(i,j).top);
                else new_piv := Remove(pivots,nd.children(i,j).bottom);
               end if;
               new_col := Remove(columns,pivots,new_piv);
               new_start := Special_Plane
                              (integer32(n),integer32(m),integer32(k),
                               new_col,special);
               One_Solve_along_Chains
                 (file,nd.children(i,j).all,n,l+1,k,ind,poset,new_piv,
                  new_col,bs,special,new_start,start,planes,report,outlog,
                  npaths,timings);
            end if;
          end loop;
        end loop;
      end if;
      One_General_Deform
        (file,n,ind,poset,nd,start,target,planes,bs,report,outlog,
         npaths,timings);
    end if;
  end One_Solve_along_Chains;

  procedure One_Quantum_Solve_along_Chains
                 ( file : in file_type; nd : in Node;
                   n,q,nb,l,k,ind : in natural32; cnt : in out natural32;
                   poset : in out Array_of_Array_of_VecMats;
                   pivots,columns : in Bracket; bs : in Bracket_System;
                   special,start,target : in Standard_Complex_Matrices.Matrix;
                   planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                   report,outlog : in boolean;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Applies the general solver along the nodes in all chains that end at
  --   the current node.  This is the analogue to the hypersurface solver,
  --   for use in connection with the grandchildren first recursive solving.
  --   This procedure is only called in the q-analogue of "One_Solve".

  -- ON ENTRY :
  --   file        to write intermediate results on;
  --   nd          current node in the localization poset;
  --   n           working dimension, equation m+p;
  --   q           degree of the map;
  --   nb          number of solution maps wanted;
  --   l           runs from 0 to k-1;
  --   k           current codimension condition;
  --   cnt         counts the number of already computed solution maps;
  --   poset       structure with all solution p-planes;
  --   ind         ind-1 planes are already folded in;
  --   pivots      pivot elements used for the special m-plane;
  --   columns     which columns of the special m-plane are used;
  --   bs          Laplace expansion of the polynomial equations;
  --   special     special m-plane for top pivots;
  --   start       (m+1-k)-plane used at the start of the deformation;
  --   target      (m+1-k)-plane used as target;
  --   planes      sequence of (m+1-k)-planes;
  --   s           interpolation points where the map meets the planes;
  --   report      indicates whether intermediate output during continuation;
  --   outlog      flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   cnt         updated counter of solution maps;
  --   poset       updated structure of all solution p-planes;
  --   npaths      updated numbers of paths traced at each level;
  --   timings     updated CPU user times at each level.

    m : constant natural32 := n - natural32(nd.p);
    new_piv,new_col : Bracket(1..pivots'last-1);
    new_start : Standard_Complex_Matrices.Matrix
                  (1..integer32(n),start'range(2));
    mod_piv : Bracket(1..nd.p);

  begin
    if empty(poset,nd.level,nd.label) then
      if l < k-1 then
        for i in nd.children'range(1) loop
          for j in nd.children'range(2) loop
            if ((nd.children(i,j) /= null)
              and then (nd.children(i,j).roco > 0)) then
               if nd.children(i,j).tp = top then
                 mod_piv := Modulo(nd.children(i,j).top,n);
                 new_piv := Remove(pivots,mod_piv);
                 put(file,"Top pivots at node : ");
                 put(file,nd.top);
                 put(file,"  child top pivots : ");
                 put(file,nd.children(i,j).top); new_line(file);
               else
                 mod_piv := Modulo(nd.children(i,j).bottom,n);
                 new_piv := Remove(pivots,mod_piv);
                 put(file,"Bottom pivots at node : ");
                 put(file,nd.bottom);
                 put(file,"  child bottom pivots : ");
                 put(file,nd.children(i,j).bottom); new_line(file);
               end if;
               put(file,"Modular pivots : "); put(file,mod_piv);
               put(file,"  new pivots : "); put(file,new_piv);
               new_line(file);
               put(file,"Pivot columns : "); put(file,columns);
               put(file,"  new columns : "); put(file,new_piv);
               new_line(file);
               new_col := Remove(columns,pivots,new_piv);
               new_start := Special_Plane
                              (integer32(n),integer32(m),integer32(k),
                               new_col,special);
               One_Quantum_Solve_along_Chains
                 (file,nd.children(i,j).all,n,q,nb,l+1,k,ind,cnt,
                  poset,new_piv,new_col,bs,special,new_start,start,
                  planes,s,report,outlog,npaths,timings);
            end if;
          end loop;
        end loop;
      end if;
      if nd.tp = top
       then mod_piv := Modulo(nd.top,n);
       else mod_piv := Modulo(nd.bottom,n);
      end if;
      One_General_Quantum_Deform
        (file,n,q,nb,l,k,ind,cnt,poset,nd,mod_piv,start,target,planes,s,bs,
         report,outlog,npaths,timings);
    end if;
  end One_Quantum_Solve_along_Chains;

  procedure Solve_along_One_Chain
                 ( file : in file_type; nd : in Node; n,l,k,ind : in natural32;
                   cod : in Bracket; poset : in out Array_of_Array_of_VecMats;
                   pivots,columns : in Bracket; bs : in Bracket_System;
                   special,start,target : in Standard_Complex_Matrices.Matrix;
                   planes : in VecMat; report,outlog : in boolean;
                   npaths : in out Standard_Natural_Vectors.Vector;
                   timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Applies the general solver along the nodes in all chains that end at
  --   the current node.  This is the analogue to the hypersurface solver,
  --   which can be used in a general control structure.

  -- ON ENTRY :
  --   file        to write intermediate results on;
  --   nd          current node in the localization poset;
  --   n           working dimension, equation m+p;
  --   l           runs from 0 to k-1;
  --   k           current codimension condition;
  --   poset       structure with all solution p-planes;
  --   ind         ind-1 planes are already folded in;
  --   pivots      pivot elements used for the special m-plane;
  --   columns     which columns of the special m-plane are used;
  --   bs          Laplace expansion of the polynomial equations;
  --   special     special m-plane for top pivots;
  --   start       (m+1-k)-plane used at the start of the deformation;
  --   target      (m+1-k)-plane used as target;
  --   planes      sequence of (m+1-k)-planes;
  --   report      indicates whether intermediate output during continuation;
  --   outlog      file to write homotopies on if set to true.

  -- ON RETURN :
  --   poset       updated structure of solution p-planes;
  --   npaths      updated number of paths traced at each level;
  --   timings     updated CPU user times for each level.

    m : constant natural32 := n - natural32(nd.p);

  begin
    if empty(poset,nd.level,nd.label) then
      if nd.level = 0 then
        poset(nd.level)(nd.label)(1)
         := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n,nd));
      elsif nd.roco > 0 then
        if l = k then
          if cod'last > cod'first then
            declare
              kk : constant natural32 := cod(cod'last-1);
              kd : constant natural32 := n+1-kk;
              new_piv,new_col : Bracket(1..integer32(m));
              new_special : Standard_Complex_Matrices.Matrix
                              (1..integer32(n),1..integer32(m));
              new_target : constant Standard_Complex_Matrices.Matrix
                         := planes(cod'last-1).all;
              new_start : Standard_Complex_Matrices.Matrix
                            (1..integer32(n),1..integer32(m+1-kk));
              new_bm : Bracket_Monomial := Maximal_Minors(n,kd);
              new_bs : Bracket_System(0..integer32(Number_of_Brackets(new_bm)))
                     := Minor_Equations(kd,kd-natural32(nd.p),new_bm); 
            begin
              for i in new_col'range loop
                new_col(i) := natural32(i);
              end loop;
              if nd.tp = top  then
                new_piv := Complement(n,nd.top);
                new_special := Special_Top_Plane(integer32(m),nd.top);
              else
                new_piv := Complement(n,nd.bottom);
                new_special := Special_Bottom_Plane(integer32(m),nd.bottom);
              end if;
              new_start := Special_Plane(integer32(n),integer32(m),
                                         integer32(kk),new_col,new_special);
              Solve_along_One_Chain
                (file,nd,n,0,kk,natural32(cod'last-1),
                 cod(cod'first..cod'last-1),
                 poset,new_piv,new_col,new_bs,new_special,new_start,
                 new_target,planes,report,outlog,npaths,timings);
              Clear(new_bm); Clear(new_bs);  
            end;
          end if;
        else
          declare
            new_piv,new_col : Bracket(1..pivots'last-1);
            new_start : Standard_Complex_Matrices.Matrix
                          (1..integer32(n),start'range(2));  
          begin
            for i in nd.children'range(1) loop
              for j in nd.children'range(2) loop
                if ((nd.children(i,j) /= null)
                   and then (nd.children(i,j).roco > 0)) then
                  if nd.children(i,j).tp = top
                   then new_piv := Remove(pivots,nd.children(i,j).top);
                   else new_piv := Remove(pivots,nd.children(i,j).bottom);
                  end if;
                  new_col := Remove(columns,pivots,new_piv);
                  new_start := Special_Plane(integer32(n),integer32(m),
                                             integer32(k),new_col,special);
                  Solve_along_One_Chain
                    (file,nd.children(i,j).all,n,l+1,k,ind,cod,poset,
                     new_piv,new_col,bs,special,new_start,start,
                     planes,report,outlog,npaths,timings);
                end if;
              end loop;
            end loop;
            One_General_Deform
              (file,n,ind,poset,nd,start,target,planes,bs,
               report,outlog,npaths,timings);
          end;
        end if;
      end if;
    end if;
  end Solve_along_One_Chain;

  procedure Solve_along_Two_Chains
                  ( file : in file_type; nd : in Node;
                    n,l_top,k_top,l_bot,k_bot,ind : in natural32;
                    cod : in Bracket; poset : in out Array_of_Array_of_VecMats;
                    top_pivots,top_columns,bot_pivots,bot_columns : in Bracket;
                    top_bs,bot_bs : in Bracket_System;
                    top_special,top_start,top_target,bot_special,bot_start,
                    bot_target : in Standard_Complex_Matrices.Matrix;
                    planes : in VecMat; report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );

  -- DESCRIPTION :
  --   Applies the general solver along the nodes in all chains that end at
  --   the current node.  This is the analogue to the hypersurface solver
  --   where the type of the node may be anything.

  -- ON ENTRY :
  --   file         to write intermediate results on;
  --   nd           current node in the localization poset;
  --   n            working dimension, equation m+p;
  --   l_top        runs from 0 to k_top-1;
  --   k_top        co-dimension condition satisfied incrementing top pivots;
  --   l_bot        runs from 0 to k_bot-1;
  --   k_bot        co-dimension condition satisfied decrementing bottom pivots;
  --   poset        structure with all solution p-planes;
  --   ind          ind-1 planes are already folded in;
  --   top_pivots   pivot elements used for the special top m-plane;
  --   top_columns  which columns of the special top m-plane are used;
  --   bot_pivots   pivot elements used for the special bottom m-plane;
  --   bot_columns  which columns of the special bottom m-plane are used;
  --   top_bs       Laplace expansion of the polynomial equations;
  --   bot_bs       Laplace expansion of the polynomial equations;
  --   top_special  special m-plane for top pivots;
  --   top_start    (m+1-k)-plane used at the start of the deformation;
  --   top_target   (m+1-k)-plane used as target satisfied with top pivots;
  --   bot_special  special m-plane for top pivots;
  --   bot_start    (m+1-k)-plane used at the start of the deformation;
  --   bot_target   (m+1-k)-plane used as target satisfied with bottom pivots;
  --   planes       sequence of (m+1-k)-planes;
  --   report       indicates whether intermediate output during continuation;
  --   outlog       flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   poset        updated structure of solution p-planes;
  --   npaths       updated number of paths traced at each level;
  --   timings      updated CPU user timings at each level.

  procedure Solve_along_Two_Chains_Deforming_Top_and_Bottom
              ( file : in file_type; nd : in Node;
                n,l_top,k_top,l_bot,k_bot,ind : in natural32;
                cod : in Bracket; poset : in out Array_of_Array_of_VecMats;
                top_pivots,top_columns,bot_pivots,bot_columns : in Bracket;
                top_bs,bot_bs : in Bracket_System;
                top_special,top_start,top_target,bot_special,bot_start,
                bot_target : in Standard_Complex_Matrices.Matrix;
                planes : in VecMat; report,outlog : in boolean;
                npaths : in out Standard_Natural_Vectors.Vector;
                timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Assumes that k_top > l_top and k_bot > l_bot so that the deformations
  --   involve both incrementing top and decrementing bottom pivots.

    m : constant natural32 := n - natural32(nd.p);
    top_piv,top_col : Bracket(1..top_pivots'last-1);
    bot_piv,bot_col : Bracket(1..bot_pivots'last-1);
    new_top_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),top_start'range(2));
    new_bot_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),bot_start'range(2));

  begin
    for i in nd.children'range(1) loop
      for j in nd.children'range(2) loop
        if ((nd.children(i,j) /= null)
          and then (nd.children(i,j).roco > 0)) then
          if nd.children(i,j).tp = top or nd.children(i,j).tp = mixed then
            top_piv := Remove(top_pivots,nd.children(i,j).top);
            top_col := Remove(top_columns,top_pivots,top_piv);
            new_top_start
              := Special_Plane(integer32(n),integer32(m),integer32(k_top),
                               top_col,top_special);
          end if;
          if nd.children(i,j).tp = bottom or nd.children(i,j).tp = mixed then
            bot_piv := Remove(bot_pivots,nd.children(i,j).bottom);
            bot_col := Remove(bot_columns,bot_pivots,bot_piv);
            new_bot_start
              := Special_Plane(integer32(n),integer32(m),integer32(k_bot),
                               bot_col,bot_special);
          end if;
          Solve_along_Two_Chains
            (file,nd.children(i,j).all,
             n,l_top+1,k_top,l_bot+1,k_bot,ind,cod,poset,
             top_piv,top_col,bot_piv,bot_col,top_bs,bot_bs,
             top_special,new_top_start,top_start,
             bot_special,new_bot_start,bot_start,
             planes,report,outlog,npaths,timings);
        end if;
      end loop;
    end loop;
    Two_General_Deform
      (file,integer32(n),integer32(ind),poset,nd,top_start,top_target,
       bot_start,bot_target,planes,top_bs,bot_bs,report,outlog,npaths,timings);
  end Solve_along_Two_Chains_Deforming_Top_and_Bottom;
 
  procedure Switch_Top_and_Solve_along_Two_Chains
              ( file : in file_type; nd : in Node;
                n,l_bot,k_bot,ind : in natural32; cod : in Bracket;
                poset : in out Array_of_Array_of_VecMats;
                bot_pivots,bot_columns : in Bracket;
                bot_bs : in Bracket_System;
                bot_special,bot_start,bot_target
                  : in Standard_Complex_Matrices.Matrix;
                planes : in VecMat; report,outlog : in boolean;
                npaths : in out Standard_Natural_Vectors.Vector;
                timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Assumes that l_top = k_top, l_bot < k_bot, and ind > cod'first so
  --   that first a new top chain can be started, which is then solved
  --   along with the existing chain for decrementing bottom pivots.

    m : constant natural32 := n - natural32(nd.p);
    new_k_top : constant natural32 := cod(integer32(ind));
    kd : constant natural32 := n+1-new_k_top;
    new_top_pivots : constant Bracket(1..integer32(m)) := Complement(n,nd.top);
    new_top_special : constant Standard_Complex_Matrices.Matrix
                                 (1..integer32(n),1..integer32(m))
                    := Special_Top_Plane(integer32(m),nd.top);
    new_top_target : constant Standard_Complex_Matrices.Matrix
                                (1..integer32(n),1..integer32(m+1-new_k_top))
      := planes(integer32(ind)).all;
    new_top_columns : Bracket(1..integer32(m));
    new_top_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),1..integer32(m+1-new_k_top));
    new_top_bm : Bracket_Monomial := Maximal_Minors(n,kd);
    new_top_bs : Bracket_System(0..integer32(Number_of_Brackets(new_top_bm)))
               := Minor_Equations(kd,kd-natural32(nd.p),new_top_bm);

  begin
    for i in new_top_columns'range loop
      new_top_columns(i) := natural32(i);
    end loop;
    new_top_start
      := Special_Plane(integer32(n),integer32(m),integer32(new_k_top),
                       new_top_columns,new_top_special);
    Solve_along_Two_Chains_Deforming_Top_and_Bottom
      (file,nd,n,0,new_k_top,l_bot,k_bot,ind,cod,poset,
       new_top_pivots,new_top_columns,bot_pivots,bot_columns,
       new_top_bs,bot_bs,new_top_special,new_top_start,new_top_target,
       bot_special,bot_start,bot_target,planes,
       report,outlog,npaths,timings);
    Clear(new_top_bm); Clear(new_top_bs);
  end Switch_Top_and_Solve_along_Two_Chains;
 
  procedure Switch_Top_and_Solve_along_One_Chain
              ( file : in file_type; nd : in Node; n,ind : in natural32;
                 cod : in Bracket; poset : in out Array_of_Array_of_VecMats;
                 planes : in VecMat; report,outlog : in boolean;
                 npaths : in out Standard_Natural_Vectors.Vector;
                 timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Assumes that l_top = k_top, l_bot = k_bot, with nd.tp = top
  --   and ind > cod'first so that a new top chain can be started.

    m : constant natural32 := n - natural32(nd.p);
    new_k_top : constant natural32 := cod(integer32(ind));
    kd : constant natural32 := n+1-new_k_top;
    new_top_pivots : constant Bracket(1..integer32(m)) := Complement(n,nd.top);
    new_top_special : constant Standard_Complex_Matrices.Matrix
                                 (1..integer32(n),1..integer32(m))
                    := Special_Top_Plane(integer32(m),nd.top);
    new_top_target : constant Standard_Complex_Matrices.Matrix
                                (1..integer32(n),1..integer32(m+1-new_k_top))
      := planes(integer32(ind)).all;
    new_top_columns : Bracket(1..integer32(m));
    new_top_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),1..integer32(m+1-new_k_top));
    new_top_bm : Bracket_Monomial := Maximal_Minors(n,kd);
    new_top_bs : Bracket_System(0..integer32(Number_of_Brackets(new_top_bm)))
               := Minor_Equations(kd,kd-natural32(nd.p),new_top_bm);

  begin
    for i in new_top_columns'range loop
      new_top_columns(i) := natural32(i);
    end loop;
    new_top_start
      := Special_Plane(integer32(n),integer32(m),integer32(new_k_top),
                       new_top_columns,new_top_special);
    Solve_along_One_Chain
      (file,nd,n,0,new_k_top,ind,cod(cod'first..integer32(ind)),poset,
       new_top_pivots,new_top_columns,new_top_bs,
       new_top_special,new_top_start,new_top_target,
       planes,report,outlog,npaths,timings);
    Clear(new_top_bm); Clear(new_top_bs);
  end Switch_Top_and_Solve_along_One_Chain;
 
  procedure Switch_Bottom_and_Solve_along_Two_Chains
              ( file : in file_type; nd : in Node;
                n,l_top,k_top,ind : in natural32; cod : in Bracket;
                poset : in out Array_of_Array_of_VecMats;
                top_pivots,top_columns : in Bracket;
                top_bs : in Bracket_System;
                top_special,top_start,top_target
                  : in Standard_Complex_Matrices.Matrix;
                planes : in VecMat; report,outlog : in boolean;
                npaths : in out Standard_Natural_Vectors.Vector;
                timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Assumes that l_top < k_top, l_bot = k_bot, and ind > cod'first so
  --   that first a new bottom chain can be launched, which is then solved
  --   along with the existing chain for incrementing top pivots.

    m : constant natural32 := n - natural32(nd.p);
    new_k_bot : constant natural32 := cod(integer32(ind));
    kd : constant natural32 := n+1-new_k_bot;
    new_bot_pivots : constant Bracket(1..integer32(m))
                   := Complement(n,nd.bottom);
    new_bot_special : constant Standard_Complex_Matrices.Matrix
                                 (1..integer32(n),1..integer32(m))
                    := Special_Bottom_Plane(integer32(m),nd.bottom);
    new_bot_target : constant Standard_Complex_Matrices.Matrix
                                (1..integer32(n),1..integer32(m+1-new_k_bot))
                   := planes(integer32(ind)).all;
    new_bot_columns : Bracket(1..integer32(m));
    new_bot_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),1..integer32(m+1-new_k_bot));
    new_bot_bm : Bracket_Monomial := Maximal_Minors(n,kd);
    new_bot_bs : Bracket_System(0..integer32(Number_of_Brackets(new_bot_bm)))
               := Minor_Equations(kd,kd-natural32(nd.p),new_bot_bm);

  begin
    for i in new_bot_columns'range loop
      new_bot_columns(i) := natural32(i);
    end loop;
    new_bot_start
      := Special_Plane(integer32(n),integer32(m),integer32(new_k_bot),
                       new_bot_columns,new_bot_special);
    Solve_along_Two_Chains_Deforming_Top_and_Bottom
      (file,nd,n,l_top,k_top,0,new_k_bot,ind,cod,poset,
       top_pivots,top_columns,new_bot_pivots,new_bot_columns,
       top_bs,new_bot_bs,top_special,top_start,top_target,
       new_bot_special,new_bot_start,new_bot_target,
       planes,report,outlog,npaths,timings); 
    Clear(new_bot_bm); Clear(new_bot_bs);
  end Switch_Bottom_and_Solve_along_Two_Chains;
 
  procedure Switch_Bottom_and_Solve_along_One_Chain
              ( file : in file_type; nd : in Node; n,ind : in natural32;
                cod : in Bracket; poset : in out Array_of_Array_of_VecMats;
                planes : in VecMat; report,outlog : in boolean;
                npaths : in out Standard_Natural_Vectors.Vector;
                timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Assumes that l_top < k_top, l_bot = k_bot, and ind > cod'first so
  --   that first a new bottom chain can be launched, which is then solved
  --   along with the existing chain for incrementing top pivots.

    m : constant natural32 := n - natural32(nd.p);
    new_k_bot : constant natural32 := cod(integer32(ind));
    kd : constant natural32 := n+1-new_k_bot;
    new_bot_pivots : constant Bracket(1..integer32(m))
                   := Complement(n,nd.bottom);
    new_bot_special : constant Standard_Complex_Matrices.Matrix
                                 (1..integer32(n),1..integer32(m))
                    := Special_Bottom_Plane(integer32(m),nd.bottom);
    new_bot_target
      : constant Standard_Complex_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-new_k_bot))
      := planes(integer32(ind)).all;
    new_bot_columns : Bracket(1..integer32(m));
    new_bot_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),1..integer32(m+1-new_k_bot));
    new_bot_bm : Bracket_Monomial := Maximal_Minors(n,kd);
    new_bot_bs : Bracket_System(0..integer32(Number_of_Brackets(new_bot_bm)))
               := Minor_Equations(kd,kd-natural32(nd.p),new_bot_bm);

  begin
    for i in new_bot_columns'range loop
      new_bot_columns(i) := natural32(i);
    end loop;
    new_bot_start
      := Special_Plane(integer32(n),integer32(m),integer32(new_k_bot),
                       new_bot_columns,new_bot_special);
    Solve_along_One_Chain
      (file,nd,n,0,new_k_bot,ind,cod(cod'first..integer32(ind)),poset,
       new_bot_pivots,new_bot_columns,new_bot_bs,
       new_bot_special,new_bot_start,new_bot_target,
       planes,report,outlog,npaths,timings); 
    Clear(new_bot_bm); Clear(new_bot_bs);
  end Switch_Bottom_and_Solve_along_One_Chain;

  procedure Switch_Top_Bottom_and_Solve_along_Two_Chains
              ( file : in file_type; nd : in Node; n,ind : in natural32;
                cod : in Bracket; poset : in out Array_of_Array_of_VecMats;
                planes : in VecMat; report,outlog : in boolean;
                npaths : in out Standard_Natural_Vectors.Vector;
                timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Assumes that l_top = k_top, l_bot = k_bot, and ind > cod'first+1
  --   so that first new top and bottom chains can be started which are
  --   then solved along.

    m : constant natural32 := n - natural32(nd.p);
    new_k_top : constant natural32 := cod(integer32(ind));
    new_k_bot : constant natural32 := cod(integer32(ind)+1);
    kd_top : constant natural32 := n+1-new_k_top;
    kd_bot : constant natural32 := n+1-new_k_bot;
    new_top_pivots : constant Bracket(1..integer32(m))
                   := Complement(n,nd.top);
    new_bot_pivots : constant Bracket(1..integer32(m))
                   := Complement(n,nd.bottom);
    new_top_special : constant Standard_Complex_Matrices.Matrix
                                 (1..integer32(n),1..integer32(m))
                    := Special_Top_Plane(integer32(m),nd.top);
    new_bot_special : constant Standard_Complex_Matrices.Matrix
                                 (1..integer32(n),1..integer32(m))
                    := Special_Bottom_Plane(integer32(m),nd.bottom);
    new_top_target
      : constant Standard_Complex_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-new_k_top))
      := planes(integer32(ind)).all;
    new_bot_target
      : constant Standard_Complex_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-new_k_bot))
      := planes(integer32(ind)+1).all;
    new_top_columns,new_bot_columns : Bracket(1..integer32(m));
    new_top_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),1..integer32(m+1-new_k_top));
    new_bot_start : Standard_Complex_Matrices.Matrix
                      (1..integer32(n),1..integer32(m+1-new_k_bot));
    new_top_bm : Bracket_Monomial := Maximal_Minors(n,kd_top);
    new_bot_bm : Bracket_Monomial := Maximal_Minors(n,kd_bot);
    new_top_bs : Bracket_System(0..integer32(Number_of_Brackets(new_top_bm)))
               := Minor_Equations(kd_top,kd_top-natural32(nd.p),new_top_bm);
    new_bot_bs : Bracket_System(0..integer32(Number_of_Brackets(new_bot_bm)))
               := Minor_Equations(kd_bot,kd_bot-natural32(nd.p),new_bot_bm);

  begin
    for i in new_bot_columns'range loop
      new_top_columns(i) := natural32(i);
      new_bot_columns(i) := natural32(i);
    end loop;
    new_top_start := Special_Plane
                       (integer32(n),integer32(m),integer32(new_k_top),
                        new_top_columns,new_top_special);
    new_bot_start := Special_Plane
                       (integer32(n),integer32(m),integer32(new_k_bot),
                        new_bot_columns,new_bot_special);
    Solve_along_Two_Chains_Deforming_Top_and_Bottom
      (file,nd,n,0,new_k_top,0,new_k_bot,ind,cod,poset,
       new_top_pivots,new_top_columns,new_bot_pivots,new_bot_columns,
       new_top_bs,new_bot_bs,new_top_special,new_top_start,new_top_target,
       new_bot_special,new_bot_start,new_bot_target,
       planes,report,outlog,npaths,timings);
    Clear(new_top_bm); Clear(new_top_bs);
    Clear(new_bot_bm); Clear(new_bot_bs);
  end Switch_Top_Bottom_and_Solve_along_Two_Chains;

  procedure Solve_along_Two_Chains
                  ( file : in file_type; nd : in Node;
                    n,l_top,k_top,l_bot,k_bot,ind : in natural32;
                    cod : in Bracket; poset : in out Array_of_Array_of_VecMats;
                    top_pivots,top_columns,bot_pivots,bot_columns : in Bracket;
                    top_bs,bot_bs : in Bracket_System;
                    top_special,top_start,top_target,bot_special,bot_start,
                    bot_target : in Standard_Complex_Matrices.Matrix;
                    planes : in VecMat; report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   Does the dispatching according to the case analysis.

  -- IMPORTANT :
  --   The control structure in the case analysis matches the structure
  --   in Localization_Posets.Recursive_Top_Bottom_Create.

  begin
    if empty(poset,nd.level,nd.label) then
      if nd.level = 0 then
        poset(nd.level)(nd.label)(1)
          := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n,nd));
      elsif nd.roco > 0 then
        if ((l_top < k_top) and (l_bot < k_bot))then
          Solve_along_Two_Chains_Deforming_Top_and_Bottom
            (file,nd,n,l_top,k_top,l_bot,k_bot,ind,cod,poset,
             top_pivots,top_columns,bot_pivots,bot_columns,
             top_bs,bot_bs,top_special,top_start,top_target,
             bot_special,bot_start,bot_target,planes,
             report,outlog,npaths,timings);
        elsif ((l_top = k_top) and (l_bot < k_bot)) then
          if integer32(ind) = cod'first then
            Solve_along_One_Chain
              (file,nd,n,l_bot,k_bot,ind,cod,poset,
               bot_pivots,bot_columns,bot_bs,bot_special,
               bot_start,bot_target,planes,report,outlog,npaths,timings);
          else
            Switch_Top_and_Solve_along_Two_Chains
              (file,nd,n,l_bot,k_bot,ind-1,cod,poset,
               bot_pivots,bot_columns,bot_bs,bot_special,
               bot_start,bot_target,planes,report,outlog,npaths,timings);
          end if;
        elsif ((l_top < k_top) and (l_bot = k_bot)) then
          if integer32(ind) = cod'first then
            Solve_along_One_Chain
              (file,nd,n,l_top,k_top,ind,cod,poset,
               top_pivots,top_columns,top_bs,top_special,
               top_start,top_target,planes,report,outlog,npaths,timings);
          else
            Switch_Bottom_and_Solve_along_Two_Chains
              (file,nd,n,l_top,k_top,ind-1,cod,poset,
               top_pivots,top_columns,top_bs,top_special,
               top_start,top_target,planes,report,outlog,npaths,timings);
          end if;
        else -- ((l_top = k_top) and (l_bot = k_bot))
          if integer32(ind) > cod'first+1 then
            Switch_Top_Bottom_and_Solve_along_Two_Chains
              (file,nd,n,ind-2,cod,poset,planes,report,outlog,npaths,timings);
          elsif integer32(ind) > cod'first then
            if nd.tp = bottom then
              Switch_Bottom_and_Solve_along_One_Chain
                (file,nd,n,ind-1,cod,poset,planes,
                 report,outlog,npaths,timings);
            else
              Switch_Top_and_Solve_along_One_Chain
                (file,nd,n,ind-1,cod,poset,planes,
                 report,outlog,npaths,timings);
            end if;
          end if;
        end if;
      end if;
    end if;
  end Solve_along_Two_Chains;

  procedure One_Solve
                  ( file : in file_type; n : in natural32; cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats;
                    nd : in Node; planes : in VecMat;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

  -- NOTE :
  --   We assume that we are always folding in the last condition, with
  --   codimension cod(cod'last).  In passing recursively to higher levels
  --   in the deformation poset, we take a slice of k, omitting the last one.
  --   Applies the solver first to all grandchildren of the current node,
  --   which is the additional layer compared to the hypersurface case.
  --   This implementation will only work in the non-mixed case.

    m : constant natural32 := n - natural32(nd.p);
    kk : constant natural32 := cod(cod'last);
    kd : constant natural32 := n+1-kk;

    procedure Solve_Grand_Child
                ( lnd : in Link_to_Node; continue : out boolean ) is

    -- DESCRIPTION :
    --   This node lnd is a grandchild of the current node.

    begin
      if Empty(poset,lnd.level,lnd.label) then
        if lnd.level = 0 then
          poset(lnd.level)(lnd.label)(1)
            := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n,lnd.all)); 
        elsif lnd.roco > 0 then
          One_Solve(file,n,cod(cod'first..cod'last-1),poset,
                    lnd.all,planes,report,outlog,npaths,timings);
        end if;
      end if;
      continue := true;
    end Solve_Grand_Child;

    procedure Solve_Grand_Children is
      new Enumerate_Grand_Children(Solve_Grand_Child);

  begin
    if (Empty(poset,nd.level,nd.label) and (nd.roco > 0))
     then
       if cod'last >= cod'first
        then Solve_Grand_Children(nd,kk);
       end if;
       declare
         pivots,columns : Bracket(1..integer32(m));
         special : Standard_Complex_Matrices.Matrix
                     (1..integer32(n),1..integer32(m));
         target : constant Standard_Complex_Matrices.Matrix
                := planes(cod'last).all;
         start : Standard_Complex_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-kk));
         bm : Bracket_Monomial := Maximal_Minors(n,kd);
         bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
            := Minor_Equations(kd,kd-natural32(nd.p),bm); 
       begin
         for i in columns'range loop
           columns(i) := natural32(i);
         end loop;
         if nd.tp = top then
           pivots := Complement(n,nd.top);
           special := Special_Top_Plane(integer32(m),nd.top);
         else
           pivots := Complement(n,nd.bottom);
           special := Special_Bottom_Plane(integer32(m),nd.bottom);
         end if;
         start := Special_Plane(integer32(n),integer32(m),integer32(kk),
                                columns,special);
         One_Solve_along_Chains
           (file,nd,n,0,kk,natural32(cod'last),poset,pivots,columns,bs,
            special,start,target,planes,report,outlog,npaths,timings);
         Clear(bm); Clear(bs);
       end;
    end if;
  end One_Solve;

--  procedure Chain_Solve
--                  ( file : in file_type; n : in natural; cod : in Bracket;
--                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
--                    planes : in VecMat; report,outlog : in boolean;
--                    npaths : in out Standard_Natural_Vectors.Vector;
--                    timings : in out Duration_Array ) is
--
--  -- NOTE :
--  --   The convention is that the last co-dimension condition is treated
--  --   when the type of the node is not mixed, otherwise the last two entries
--  --   of the vector of co-dimension conditions are sliced off when moving
--  --   to the upper levels.
--  --   This is another organization of One_Solve and only works when the
--  --   type of the nodes are not mixed.
--
--  begin
--    if Empty(poset,nd.level,nd.label) then
--      if nd.level = 0 then
--        poset(nd.level)(nd.label)(1)
--          := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n,nd));
--      elsif nd.roco > 0 then
--        declare
--          m : constant natural := n - nd.p;
--          pivots,columns : Bracket(1..m);
--          special : Standard_Complex_Matrices.Matrix(1..n,1..m);
--          kk : constant natural := cod(cod'last);
--          kd : constant natural := n+1-kk;
--          start : Standard_Complex_Matrices.Matrix(1..n,1..m+1-kk);
--          target : constant Standard_Complex_Matrices.Matrix
--                 := planes(planes'last).all;
--          bm : Bracket_Monomial := Maximal_Minors(n,kd);
--          bs : Bracket_System(0..Number_of_Brackets(bm))
--             := Minor_Equations(kd,kd-nd.p,bm);
--        begin
--          for i in columns'range loop
--            columns(i) := i;
--          end loop;
--          if nd.tp = top
--           then pivots := Complement(n,nd.top);
--                special := Special_Top_Plane(m,nd.top);
--           else pivots := Complement(n,nd.bottom);
--                special := Special_Bottom_Plane(m,nd.bottom);
--          end if;
--          start := Special_Plane(n,m,kk,columns,special);
--          Solve_along_One_Chain
--            (file,nd,n,0,cod(cod'last),cod'last,cod,poset,pivots,columns,
--             bs,special,start,target,planes,report,outlog,npaths,timings);
--          Clear(bm); Clear(bs);
--        end;
--      end if;
--    end if;
--  end Chain_Solve;

  procedure Solve ( file : in file_type; n : in natural32; cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

  -- REQUIREMENT ONE :
  --   The convention is that the last co-dimension condition is treated
  --   when the type of the node is not mixed, otherwise the last two entries
  --   of the vector of co-dimension conditions are sliced off when moving
  --   to the upper levels.
  -- REQUIREMENT TWO :
  --   The nodes that are not mixed appear at the top of the poset.

  begin
    if Empty(poset,nd.level,nd.label) then
      if nd.level = 0
       then poset(nd.level)(nd.label)(1)
              := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n,nd));
       elsif nd.roco > 0 then
         if nd.tp /= mixed then
           declare
             m : constant natural32 := n - natural32(nd.p);
             pivots,columns : Bracket(1..integer32(m));
             special : Standard_Complex_Matrices.Matrix
                         (1..integer32(n),1..integer32(m));
             kk : constant natural32 := cod(cod'last);
             kd : constant natural32 := n+1-kk;
             start : Standard_Complex_Matrices.Matrix
                       (1..integer32(n),1..integer32(m+1-kk));
             target : constant Standard_Complex_Matrices.Matrix
                    := planes(planes'last).all;
             bm : Bracket_Monomial := Maximal_Minors(n,kd);
             bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
                := Minor_Equations(kd,kd-natural32(nd.p),bm);
           begin
             for i in columns'range loop
               columns(i) := natural32(i);
             end loop;
             if nd.tp = top then
               pivots := Complement(n,nd.top);
               special := Special_Top_Plane(integer32(m),nd.top);
             else
               pivots := Complement(n,nd.bottom);
               special := Special_Bottom_Plane(integer32(m),nd.bottom);
             end if;
             start := Special_Plane(integer32(n),integer32(m),integer32(kk),
                                    columns,special);
             Solve_along_One_Chain
               (file,nd,n,0,cod(cod'last),natural32(cod'last),
                cod,poset,pivots,columns,
                bs,special,start,target,planes,report,outlog,npaths,timings);
             Clear(bm); Clear(bs);
           end;
         else
           declare
             m : constant natural32 := n - natural32(nd.p);
             top_col,bot_col : Bracket(1..integer32(m));
             kk_top : constant natural32 := cod(cod'last-1);
             kk_bot : constant natural32 := cod(cod'last);
             kd_top : constant natural32 := n+1-kk_top;
             kd_bot : constant natural32 := n+1-kk_bot;
             top_bm : Bracket_Monomial := Maximal_Minors(n,kd_top);
             bot_bm : Bracket_Monomial := Maximal_Minors(n,kd_bot);
             top_bs : Bracket_System(0..integer32(Number_of_Brackets(top_bm)))
                    := Minor_Equations(kd_top,kd_top-natural32(nd.p),top_bm);
             bot_bs : Bracket_System(0..integer32(Number_of_Brackets(bot_bm)))
                    := Minor_Equations(kd_bot,kd_bot-natural32(nd.p),bot_bm);
             top_special : constant Standard_Complex_Matrices.Matrix
                                      (1..integer32(n),1..integer32(m))
                         := Special_Top_Plane(integer32(m),nd.top);
             bot_special : constant Standard_Complex_Matrices.Matrix
                                      (1..integer32(n),1..integer32(m))
                         := Special_Bottom_Plane(integer32(m),nd.bottom);
             top_piv : constant Bracket(1..integer32(m))
                     := Complement(n,nd.top);
             bot_piv : constant Bracket(1..integer32(m))
                     := Complement(n,nd.bottom);
             top_start,top_target 
               : Standard_Complex_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-kk_top));
             bot_start,bot_target
               : Standard_Complex_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-kk_bot));
           begin
             for i in top_col'range loop
               top_col(i) := natural32(i);
             end loop;
             top_start := Special_Plane(integer32(n),integer32(m),
                                        integer32(kk_top),top_col,top_special);
             top_target := planes(planes'last-1).all;
             for i in bot_col'range loop
               bot_col(i) := natural32(i);
             end loop;
             bot_start := Special_Plane(integer32(n),integer32(m),
                                        integer32(kk_bot),bot_col,bot_special);
             bot_target := planes(planes'last).all;
             Solve_along_Two_Chains
               (file,nd,n,0,kk_top,0,kk_bot,natural32(cod'last-1),cod,poset,
                top_piv,top_col,bot_piv,bot_col,top_bs,bot_bs,
                top_special,top_start,top_target,bot_special,bot_start,
                bot_target,planes,report,outlog,npaths,timings);
             Clear(top_bm); Clear(top_bs); Clear(bot_bm); Clear(bot_bs);
           end;
         end if;
       end if;
    end if;
  end Solve;

  procedure Recursive_Quantum_Solve
               ( n,q,nb : in natural32; cnt : in out natural32;
                 nd : in Node; expbp : in Bracket_Polynomial;
                 poset : in out Array_of_Array_of_VecMats;
                 planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                 npaths : in out Standard_Natural_Vectors.Vector;
                 timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   This additional layer is added to avoid the repeated construction
  --   of the structure of the equations, that is now in expbp.
  --   This is the q-analogue to the Recursive Hypersurface Solver.

    back_cnt : constant natural32 := cnt;

  begin
    if Empty(poset,nd.level,nd.label)
     then
       if nd.level = 0
        then
          poset(nd.level)(nd.label)(1)
            := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n*(q+1),nd));
        else
          for i in nd.children'range(1) loop
            for j in nd.children'range(2) loop
              if ((cnt < nb) and then (nd.children(i,j) /= null))
               then Recursive_Quantum_Solve
                      (n,q,nb,cnt,nd.children(i,j).all,expbp,poset,
                       planes,s,npaths,timings);
              end if;
              exit when (cnt = nb);
            end loop;
          end loop;
          cnt := back_cnt;
          Quantum_Deform(n,q,nb,cnt,poset,nd,expbp,planes,s,npaths,timings);
       end if;
    end if;
  end Recursive_Quantum_Solve;

  procedure Recursive_Quantum_Solve
               ( file : in file_type;
                 n,q,nb : in natural32; cnt : in out natural32;
                 nd : in Node; expbp : in Bracket_Polynomial;
                 poset : in out Array_of_Array_of_VecMats;
                 planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                 report,outlog : in boolean;
                 npaths : in out Standard_Natural_Vectors.Vector;
                 timings : in out Duration_Array ) is

  -- DESCRIPTION :
  --   This additional layer is added to avoid the repeated construction
  --   of the structure of the equations, that is now in expbp.
  --   This is the q-analogue to the Recursive Hypersurface Solver.

    back_cnt : constant natural32 := cnt;

  begin
    if Empty(poset,nd.level,nd.label) then
      if nd.level = 0 then
        poset(nd.level)(nd.label)(1)
          := new Standard_Complex_Matrices.Matrix'(Leaf_Plane(n*(q+1),nd));
      else
        for i in nd.children'range(1) loop
          for j in nd.children'range(2) loop
            if ((cnt < nb) and then (nd.children(i,j) /= null)) then
              Recursive_Quantum_Solve
                (file,n,q,nb,cnt,nd.children(i,j).all,expbp,poset,
                 planes,s,report,outlog,npaths,timings);
            end if;
            exit when (cnt = nb);
          end loop;
        end loop;
        cnt := back_cnt;
        Quantum_Deform(file,n,q,nb,cnt,poset,nd,expbp,planes,s,
                       report,outlog,npaths,timings);
       end if;
    end if;
  end Recursive_Quantum_Solve;

  procedure Solve ( n,q,nb : in natural32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

    bm : Bracket_Monomial := Maximal_Minors(n,n);
    bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(n,n-natural32(nd.p),bm);
    cnt : natural32 := 0;

  begin
    Recursive_Quantum_Solve
      (n,q,nb,cnt,nd,bs(1),poset,planes,s,npaths,timings);
    Clear(bm); Clear(bs);
  end Solve;

  procedure Solve ( file : in file_type; n,q,nb : in natural32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

    bm : Bracket_Monomial := Maximal_Minors(n,n);
    bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(n,n-natural32(nd.p),bm);
    cnt : natural32 := 0;

  begin
    Recursive_Quantum_Solve
      (file,n,q,nb,cnt,nd,bs(1),poset,planes,s,report,outlog,npaths,timings);
    Clear(bm); Clear(bs);
  end Solve;

  procedure One_Solve
                  ( file : in file_type;
                    n,q,nb : in natural32; cnt : in out natural32;
                    cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

    m : constant natural32 := n - natural32(nd.p);
    kk : constant natural32 := cod(cod'last);
    kd : constant natural32 := n+1-kk;

    procedure Solve_Grand_Child
                ( lnd : in Link_to_Node; continue : out boolean ) is

    -- DESCRIPTION :
    --   This node lnd is a grandchild of the current node.

    begin
      if Empty(poset,lnd.level,lnd.label) then
        if lnd.level = 0 then
          poset(lnd.level)(lnd.label)(1)
            := new Standard_Complex_Matrices.Matrix'
                     (Leaf_Plane(n*(q+1),lnd.all)); 
        elsif lnd.roco > 0 then
          One_Solve
            (file,n,q,nb,cnt,cod(cod'first..cod'last-1),poset,
             lnd.all,planes,s,report,outlog,npaths,timings);
        end if;
      end if;
      continue := true;
    end Solve_Grand_Child;

    procedure Solve_Grand_Children is
      new Enumerate_Grand_Children(Solve_Grand_Child);

  begin
    if (Empty(poset,nd.level,nd.label) and (nd.roco > 0)) then
      if cod'last >= cod'first
       then Solve_Grand_Children(nd,kk);
      end if;
      declare
        pivots,columns : Bracket(1..integer32(m));
        mod_piv : Bracket(1..nd.p);
        special : Standard_Complex_Matrices.Matrix
                    (1..integer32(n),1..integer32(m));
        target : constant Standard_Complex_Matrices.Matrix
               := planes(cod'last).all;
        start : Standard_Complex_Matrices.Matrix
                  (1..integer32(n),1..integer32(m+1-kk));
        bm : Bracket_Monomial := Maximal_Minors(n,kd);
        bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
           := Minor_Equations(kd,kd-natural32(nd.p),bm); 
      begin
        for i in columns'range loop
          columns(i) := natural32(i);
        end loop;
        if nd.tp = top then
          mod_piv := Modulo(nd.top,n);
          pivots := Complement(n,mod_piv);
          special := Special_Top_Plane(integer32(m),mod_piv);
        else
          mod_piv := Modulo(nd.bottom,n);
          pivots := Complement(n,mod_piv);
          special := Special_Bottom_Plane(integer32(m),mod_piv);
        end if;
        start := Special_Plane(integer32(n),integer32(m),integer32(kk),
                               columns,special);
        One_Quantum_Solve_along_Chains
          (file,nd,n,q,nb,0,kk,natural32(cod'last),cnt,poset,pivots,columns,
           bs,special,start,target,planes,s,report,outlog,npaths,timings);
        Clear(bm); Clear(bs);
      end;
    end if;
  end One_Solve;

  procedure Solve ( file : in file_type;
                    n,q,nb : in natural32; cnt : in out natural32;
                    cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array ) is

  begin
    null;
  end Solve;

-- DESTRUCTORS :

  procedure Clear ( avm : in out Array_of_VecMats ) is
  begin
    for i in avm'range loop
      Deep_Clear(avm(i));
    end loop;
  end Clear;

  procedure Clear ( avm : in out Link_to_Array_of_VecMats ) is

    procedure free is
      new unchecked_deallocation(Array_of_VecMats,Link_to_Array_of_VecMats);

  begin
    if avm /= null
     then Clear(avm.all);
          free(avm);
    end if;
  end Clear;

  procedure Clear ( avm : in out Array_of_Array_of_VecMats ) is
  begin
    for i in avm'range loop
      if avm(i) /= null
       then Clear(avm(i).all);
      end if;
    end loop;
  end Clear;

end Deformation_Posets;
