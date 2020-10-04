with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;   
with Set_Structure,Set_Structure_io;
with Main_Vertex_Points;                 use Main_Vertex_Points;
with Trees_of_Vectors;                   use Trees_of_Vectors;
with Trees_of_Vectors_io;                use Trees_of_Vectors_io;
with Set_Structures_and_Volumes;         use Set_Structures_and_Volumes;
with Generic_Position;
with Drivers_for_Coefficient_Systems;    use Drivers_for_Coefficient_Systems;

package body Drivers_for_Implicit_Lifting is

  procedure Implicit_Lifting_Info is

    i : array(1..20) of string(1..65);

  begin
    i( 1):="  The mixed volume of the tuple  of  Newton  polytopes  gives  an";
    i( 2):="upper  bound  for  the  number of complex isolated solutions with";
    i( 3):="nonzero coefficients.  This BKK bound (called  after  Bernshtein,";
    i( 4):="Koushnirenko  and  Khovanskii) is exact when the coefficients are";
    i( 5):="sufficiently chosen at  random  or  when  the  polytopes  are  in";
    i( 6):="generic position w.r.t. each other.                              ";
    i( 7):="  By implicit lifting, the mixed volume is computed by lifting up";
    i( 8):="the origin of the first polytope and by computing those facets of";
    i( 9):="the lower hull of the Minkowski sum that are spanned by edge-edge";
    i(10):="combinations.   This  procedure  to compute the mixed volume uses";
    i(11):="recursion on the dimension.                                      ";
    i(12):="  A random coefficient start system has the same Newton polytopes";
    i(13):="as  the  target  system, but all coefficients are randomly chosen";
    i(14):="complex numbers.  To solve a  random  coefficient  start  system,";
    i(15):="homotopy  continuation  is invoked in a recursive way, similar to";
    i(16):="the computation of the mixed volume.                             ";
    i(17):="  By constructing a set structure for the first k  equations  and";
    i(18):="computing  the  mixed  volumes  of  the  remaining  n-k equations";
    i(19):="obtained  after  elimination,  a  combined  Bezout-BKK  bound  is";
    i(20):="computed.                                                        ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Implicit_Lifting_Info;

  procedure Driver_for_Mixture_Bezout_BKK
                 ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                   q : out Poly_Sys; qsols : out Solution_List;
                   b : out natural32 ) is

    welcome : constant string := "Mixed-Volume Computation by Implicit Lifting";

   -- VARIABLES :

    n,bkk,k : natural32;
    ans : character;
    timer : timing_widget;
    tdft,gft,solsft : file_type;
    td : Tree_of_Vectors;
    pp,qq : Poly_Sys(p'range);
    qqsols : Solution_List;

 -- SWITCHES TO BE SET :

    verpts   : boolean;   -- to extract the vertex points
    havetd   : boolean;   -- already a tree of directions available
    wanttd   : boolean;   -- put a tree of directions on separate file
    tosolve  : boolean;   -- a system has to be solved
    ranstart : boolean;   -- for random coefficient start system
    contrep  : boolean;   -- reporting during continuation

  begin
    new_line; put_line(welcome);
    n := natural32(p'length);
    new_line;
    put("Do you first want to extract the vertex points ? (y/n) ");
    Ask_Yes_or_No(ans);
    verpts := (ans = 'y');
    if verpts then
      put_line("Computing the vertex points...");
      Copy(p,pp);
      Vertex_Points(file,pp);
    else
      pp := p;
    end if;
    new_line;
    put_line("MENU for Mixture between Bezout and BKK Bound :");
    put_line("  k = 0 : for computation of the BKK bound;");
    put("  k = "); put(n,1); put_line(" : for a generalized Bezout number;");
    put("  0 < k < "); put(n,1); put_line(" : gives a mixed Bezout BKK bound.");
    put("Give the number k between 0 and "); put(n,1); put(" : ");
    Read_Natural(k);
    if k > 0 then
      new_line;
      put_line("Do you have a good set structure");
      put("for the first "); put(k,1); put(" polynomials ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        declare
          ns : Standard_Natural_Vectors.Vector(p'range);
        begin
          for i in 1..integer32(k) loop
	    put("  Give the number of sets for polynomial ");
	    put(i,1); put(" : "); Read_Natural(ns(i));
          end loop;
	  ns(integer32(k)+1..integer32(n))
            := (integer32(k)+1..integer32(n) => 0);
	  Set_Structure.Init(ns);
	  put_line("Give the set structure :");
	  Set_Structure_io.get;
        end;
      else
        Build_RPS(k,n,pp);
      end if;
      new_line(file);
      put_line(file,
               "******* MIXTURE BETWEEN BEZOUT AND BKK BOUND *******");
      new_line(file);
      put(file,"Set structure for the first "); put(file,k,1);
      put_line(file," polynomials : ");
      Set_Structure_io.put(file);
    end if;
    new_line;
    put("Do you have a tree of useful directions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Reading the name of the file where the tree is.");
      declare
        tdf : file_type;
      begin
        Read_Name_and_Open_File(tdf);
        get(tdf,n-k,td);
        Close(tdf);
        havetd := true;
      exception 
        when others =>
          put_line("There is something wrong with your tree ...");
          havetd := false;
      end;
    else 
      havetd := false;
    end if;
    if not havetd then
      put("Do you want the useful directions on separate file ? (y/n) ");
      Ask_Yes_or_No(ans);
      wanttd := (ans = 'y');
      if wanttd then
        put_line("Reading a name of a file to write tree on.");
        Read_Name_and_Create_File(tdft);
      end if;
    else
      wanttd := false;
    end if;
    Driver_for_Coefficient_System
      (file,pp,k,byebye,qq,gft,solsft,tosolve,ranstart,contrep);
    tstart(timer);
    if tosolve then
      new_line(file);
      put_line(file,"INTERMEDIATE RESULTS :");
      new_line(file);
      if havetd then
        Mixed_Solve(file,k,n,pp,td,bkk,qq,qqsols);
      elsif wanttd then
        Bezout_BKK(k,n,pp,td,bkk);
        put(tdft,td); Close(tdft);
        Mixed_Solve(file,k,n,pp,td,bkk,qq,qqsols);
      else
        Bezout_BKK(k,n,pp,td,bkk);
        Mixed_Solve(file,k,n,pp,td,bkk,qq,qqsols);
      end if;
      if k > 0 then
        put(gft,natural32(qq'last),qq);
        new_line(file);
        put_line(file,"THE START SYSTEM : ");
        new_line(file);
        put_line(file,qq);
      end if;
      if ranstart then
        new_line(gft); put_line(gft,"THE SOLUTIONS :"); new_line(gft);
        put(gft,Length_Of(qqsols),n,qqsols);
        Close(gft);
      else
        put(solsft,Length_Of(qqsols),n,qqsols); Close(solsft);
      end if;
      q := qq; qsols := qqsols;
    elsif havetd then
      bkk := Bezout_BKK(k,n,pp,td);
    elsif wanttd then
      Bezout_BKK(k,n,pp,td,bkk);
      put(tdft,td);
      Close(tdft);
    else
      Bezout_Bkk(k,n,pp,td,bkk);
    end if;
    tstop(timer);
    b := bkk;
    new_line(file);
    if k > 0
     then put(file,"The mixed Bezout BKK bound equals ");
     else put(file,"The BKK bound equals ");
    end if;
    put(file,bkk,1); put_line(file,".");
    new_line(file);
    if not Is_Null(td) then
      put_line(file,"The tree of useful directions : ");
      put(file,td);
      new_line(file);
      if k = 0 then
        if Generic_Position(pp,td)
         then put_line(file,"The polytopes may be in generic position.");
         else put_line(file,"The polytopes are not in generic position.");
        end if;
        new_line(file);
      end if;
    end if;
    print_times(file,timer,"mixed homotopy method");
    new_line(file);
    if ans = 'y' then
      new_line(file);
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      new_line(file);
      put(file,natural32(qq'last),qq);
      new_line(file);
      put_line(file,"THE SOLUTIONS : ");
      put(file,Length_Of(qqsols),n,qqsols);
    end if;
    if verpts
     then Clear(pp);
    end if;
  end Driver_for_Mixture_Bezout_BKK;

  procedure Driver_for_Mixture_Bezout_BKK
                ( file : in file_type; p : in Laur_Sys; byebye : in boolean;
                  q : out Laur_Sys; qsols : out Solution_List;
                  b : out natural32 ) is
 
    pp,pq : Poly_Sys(p'range);

  begin
    pp := Laurent_to_Polynomial_System(p);
    Driver_for_Mixture_Bezout_BKK(file,pp,byebye,pq,qsols,b);
    q := Polynomial_to_Laurent_System(pq);
  end Driver_for_Mixture_Bezout_BKK;

end Drivers_for_Implicit_Lifting;
