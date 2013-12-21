with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Symbolic_Symmetry_Group_io;

package body Drivers_for_Symmetry_Group_io is

  procedure Read_Permutation_Group
               ( n : in natural32; g,v : in out List_of_Permutations;
                 allperms : out boolean ) is

    ans : character;
    nb : natural32;

  begin
    new_line;
    put("Is the group the full permutation group ? (y/n) ");
    Ask_Yes_or_No(ans);  allperms := (ans = 'y');
    if ans = 'y'
     then
       g := SymGrp(integer32(n));
     else
       put("The neutral element of the group is represented as ");
       for i in 1..n loop
         declare
           sb : Symbol;
         begin
           sb := (sb'range => ' ');
           sb := Symbol_Table.Get(i);
           Symbol_Table_io.put(sb);
           put(" ");
         end;
       end loop;
       new_line;
       put("Give the number of generating elements in the group : ");
       Read_Natural(nb);
       put("Give "); put(nb,1);
       put_line(" vector representations of the generating elements :");
       Symbolic_Symmetry_Group_io.Get(g,n,nb);
    end if;
    put("Do you want the generation of the group ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then v := Generate(g);
     else v := g;
    end if;  -- v = either full group or group supplied by user
  end Read_Permutation_Group;

  procedure Read_Symmetry_Group
               ( n : in natural32; g,v : in out List_of_Permutations;
                 allperms,signsym,allsigns : out boolean ) is
   
    ans : character;
    nb : natural32;
    pv,fv,fg : List_of_Permutations;

  begin
    Read_Permutation_Group(n,pv,fv,allperms);
    put("Is there any sign symmetry to take into account ? (y/n) ");
    Ask_Yes_or_No(ans);
    signsym := (ans = 'y');
    if ans = 'y' then
      put("Contains the group all sign permutations ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        allsigns := true;
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
        Symbolic_Symmetry_Group_io.Get(fg,n,nb);
      end if;
      put("Do you want the generation of the group ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then fv := Generate(Union(fg,pv));
       else fv := Union(fg,pv);
      end if;
    else
      allsigns := false;
    end if;
    g := Union(fg,pv); v := fv;
  end Read_Symmetry_Group;

end Drivers_for_Symmetry_Group_io;
