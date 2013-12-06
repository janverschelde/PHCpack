with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Linear_Product_System;
with Set_Structure;

package body Random_Product_Start_Systems is

  type Boolean_Array is array ( natural32 range <> ) of boolean;

  function leq ( d1,d2 : Degrees ) return boolean is

  -- DESCRIPTION :
  --   returns true if for all k: d1(k) <= d2(k)

  begin
    for k in d1'range loop
      if d1(k) > d2(k) 
       then return false;
      end if;
    end loop;
    return true;
  end leq;

  function Dominated ( d : Degrees; p : Poly ) return boolean is

   -- DESCRIPTION :
   --   Returns true if the term indicated by the degrees d
   --   is already in the set structure;
   --   the polynomial p contains all terms that are already
   --   belonging to the set structure.

    res : boolean := false;

    procedure Scan ( t : in Term; continue : out boolean ) is
    begin
      if leq(d,t.dg)
       then continue := false; res := true;
       else continue := true;
      end if;
    end Scan;
    procedure Scan_Terms is new Visiting_Iterator (Scan);

  begin
    Scan_Terms(p);
    return res;
  end Dominated;

  procedure Divide ( dg : in out Standard_Natural_Vectors.Vector;
                     i,d : in natural32; occupied : in out Boolean_Array ) is

  -- DESCRIPTION :
  --   Divides the monomial by every unknown that occurs already
  --   in one set of the ith component of the set structure.
  --   Every set that contains already one unknown of the monomial is
  --   marked by setting the corresponding entry in occupied on true.

    j : natural32;

  begin
    for k in dg'range loop
      j := 1;
      while (dg(k) /= 0) and (j <= d) loop
        if not occupied(j) and Set_Structure.Is_In(i,j,natural32(k))
         then occupied(j) := true; dg(k) := dg(k)-1;
        end if;
        j := j+1;
      end loop;
    end loop;
  end Divide;

  procedure Augment ( ba : in Boolean_Array; index : in out natural32 ) is

  -- DESCRIPTION :
  --   Augments the index till ba(index) = false.

  begin
    while ba(index) loop
      index := index + 1;
      exit when index >= ba'last;
    end loop;
  end Augment;

  procedure Build_Set_Structure ( i,d : in natural32; p : in Poly ) is

    temp : Poly := Null_Poly;

    procedure Process ( t : in Term; continue : out boolean ) is

      ind : natural32 := 1; -- number of the set
      processed : Boolean_Array(1..d) := (1..d => false);
      tdg : Standard_Natural_Vectors.Vector(t.dg'range) := t.dg.all;

    begin
      if not Dominated(t.dg,temp) then
        Divide(tdg,i,d,processed);
        Augment(processed,ind);
        for k in tdg'range loop
          for j in 1..tdg(k) loop 
            if not Set_Structure.Is_In(i,ind,natural32(k)) then
              Set_Structure.Add(i,ind,natural32(k));
              processed(ind) := true;
            end if;
            Augment(processed,ind);
          end loop;
        end loop;
        Add(temp,t);
      end if;
      continue := true;
    end Process;
    procedure Process_Terms is new Visiting_Iterator (Process);

  begin
    Process_Terms(p);
    Clear(temp);
  end Build_Set_Structure;

  procedure Build_Set_Structure ( p : in Poly_Sys ) is

    n : constant integer32 := integer32(p'length);
    ns : Standard_Natural_Vectors.Vector(1..n);

  begin
    for i in p'range loop
      ns(i) := natural32(Degree(p(i)));
    end loop;
    Set_Structure.Init(ns);
    for i in p'range loop
      Build_Set_Structure(natural32(i),natural32(Degree(p(i))),p(i));
    end loop;
  end Build_Set_Structure;

  procedure Build_Random_Product_System ( n : in natural32 ) is

  -- DESCRIPTION :
  --   Based upon the set structure of p, a random product system
  --   will be constructed.

    h : Vector(0..integer32(n));
    first : boolean;

  begin
    for i in 1..n loop
      for j in 1..Set_Structure.Number_Of_Sets(i) loop
        h := (0..integer32(n) => Create(0.0));
        first := true;
	for k in 1..n loop
          if Set_Structure.Is_In(i,j,k) then
            if first
             then h(integer32(k)) := Create(1.0); first := false;
             else h(integer32(k)) := Random1;
            end if;
          end if;
        end loop;
        h(0) := Random1;
        Standard_Linear_Product_System.Add_Hyperplane(i,h);
      end loop;
    end loop;
  end Build_Random_Product_System;

  procedure Construct ( p : in Poly_Sys; q : in out Poly_Sys;
		        sols : in out Solution_List ) is

    nl : natural32 := 0;
    n : constant natural32 := natural32(p'length);

  begin
    Build_Set_Structure(p);
    Standard_Linear_Product_System.Init(n);
    Build_Random_Product_System(n);
    q := Standard_Linear_Product_System.Polynomial_System;
    Standard_Linear_Product_System.Solve(sols,nl);
    Standard_Linear_Product_System.Clear;
  end Construct;

end Random_Product_Start_Systems;
