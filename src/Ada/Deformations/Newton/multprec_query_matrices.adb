with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;

package body Multprec_Query_Matrices is

  procedure Show_Element ( A : in character; i,j : in integer32;
                           aij : in Complex_Number ) is

  begin
    put(A);
    put("("); put(i,1); put(","); put(j,1); put(") = ");
    put(aij); new_line;
  end Show_Element;

  procedure Write_Matrix ( a : in Matrix ) is
  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        Show_Element('A',i,j,a(i,j));
      end loop;
    end loop;
  end Write_Matrix;

  procedure Show_Matrix ( a : in Matrix ) is

    ans : character;
    cnt : integer32 := 0;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        cnt := cnt + 1;
        if cnt mod 20 = 0 then
          put("Do you want to see more ? (y/n) "); 
          Ask_Yes_or_No(ans);
          if ans = 'y'
           then cnt := 0;
           else cnt := 21;
          end if;
        end if;
        exit when (cnt > 20);
        Show_Element('A',i,j,a(i,j));
      end loop;
      exit when (cnt > 20);
    end loop;
  end Show_Matrix;

  procedure Query_Matrix ( a : in Matrix ) is

    row,col : integer32 := 0;

  begin
    put("Querying a "); put(a'last(1),1); put("-by-");
    put(a'last(2),1); put_line("-matrix...");
    loop
      put("Give row index (0 to exit) : "); get(row);
      exit when (row = 0);
      if row > a'last(1) then
        put("row index "); put(row,1); put(" > "); put(a'last(1),1);
        put_line("  please try again...");
      else
        put("Give column index (0 to exit) : "); get(col);
        exit when (col = 0);
        if col > a'last(2) then
          put("column index "); put(col,1); put(" > "); put(a'last(2),1);
          put_line("  please try again...");
        else
          Show_Element('A',row,col,a(row,col));
        end if;
      end if;
    end loop;
  end Query_Matrix;

  procedure Query_Matrices ( a,b : in Matrix ) is

    row,col : integer32 := 0;

  begin
    put("Querying a pair of "); put(a'last(1),1); put("-by-");
    put(a'last(2),1); put_line("-matrices...");
    loop
      put("Give row index (0 to exit) : "); get(row);
      exit when (row = 0);
      if row > a'last(1) then
        put("row index "); put(row,1); put(" > "); put(a'last(1),1);
        put_line("  please try again...");
      else
        put("Give column index (0 to exit) : "); get(col);
        exit when (col = 0);
        if col > a'last(2) then
          put("column index "); put(col,1); put(" > "); put(a'last(2),1);
          put_line("  please try again...");
        else
          Show_Element('A',row,col,a(row,col));
          Show_Element('B',row,col,b(row,col));
          put("difference = ");
          put(a(row,col)-b(row,col)); new_line;
        end if;
      end if;
    end loop;
  end Query_Matrices;

  function Difference_of_Matrices ( a,b : Matrix ) return Floating_Number is

    dif : Complex_Number;
    abs_dif : Floating_Number;
    sum : Floating_Number := Create(integer(0));

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        dif := a(i,j) - b(i,j);
        abs_dif := AbsVal(dif);
        Add(sum,abs_dif);
        Clear(dif); Clear(abs_dif);
      end loop;
    end loop;
    return sum;
  end Difference_of_Matrices;

  procedure Show_Difference_of_Matrices ( a,b : in Matrix ) is

    ans : character;
    cnt : integer32 := 0;
    dif : Complex_Number;
    abs_dif : Floating_Number;
    sum : Floating_Number := Create(integer(0));
    tol : Floating_Number := Create(1.0E-8);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        cnt := cnt + 1;
        if cnt mod 10 = 0 then
          put("Do you want to see more ? (y/n) "); 
          Ask_Yes_or_No(ans);
          if ans = 'y'
           then cnt := 0;
           else cnt := 11;
          end if;
        end if;
        exit when (cnt > 10);
        Show_Element('A',i,j,a(i,j));
        Show_Element('B',i,j,b(i,j));
        dif := a(i,j) - b(i,j);
        put("difference = "); put(dif); 
        abs_dif := AbsVal(dif);
        if abs_dif < tol then
          new_line;
        else
          put(" > "); put(tol,2); put_line("!!!");
        end if;
        Add(sum,abs_dif);
        Clear(dif); Clear(abs_dif);
      end loop;
      exit when (cnt > 10);
    end loop;
    put("Sum of all differences : "); put(sum); new_line;
    Clear(tol);
  end Show_Difference_of_Matrices;

  procedure Show_Differences_in_Matrices ( a,b : in Matrix ) is

    ans : character;
    cnt : integer32 := 1;
    dif : Complex_Number;
    abs_dif : Floating_Number;
    sum : Floating_Number := Create(integer(0));
    tol : Floating_Number := Create(1.0E-8);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        if cnt mod 10 = 0 then
          put("Do you want to see more ? (y/n) "); 
          Ask_Yes_or_No(ans);
          if ans = 'y'
           then cnt := 0;
           else cnt := 11;
          end if;
        end if;
        exit when (cnt > 10);
        dif := a(i,j) - b(i,j);
        abs_dif := AbsVal(dif);
        if abs_dif > tol then
          cnt := cnt + 1;
          Show_Element('A',i,j,a(i,j));
          Show_Element('B',i,j,b(i,j));
          put("difference = "); put(dif); 
          put(" > "); put(tol,2); put_line("!!!");
        end if;
        Add(sum,abs_dif);
        Clear(dif); Clear(abs_dif);
      end loop;
      exit when (cnt > 10);
    end loop;
    put("Sum of all differences : "); put(sum); new_line;
    Clear(tol);
  end Show_Differences_in_Matrices;

end Multprec_Query_Matrices;
