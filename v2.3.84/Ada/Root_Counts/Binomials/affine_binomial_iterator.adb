with Standard_Integer_VecVecs;

package body Affine_Binomial_Iterator is

-- STATE VARIABLES :

  incidence_matrix : Standard_Integer_Matrices.Link_to_Matrix;
  number_of_variables,number_of_monomials : integer32;
  current_level : integer32; -- top of the stack
  selections : Standard_Integer_VecVecs.Link_to_VecVec;
  nonzerovar : Standard_Integer_VecVecs.Link_to_VecVec;
  number_of_selections : Standard_Integer_Vectors.Link_to_Vector;
  number_of_equations : Standard_Integer_Vectors.Link_to_Vector;
  maximum_selections : integer32;
  -- generated,generated_last : List; -- to avoid duplicate reporting
  -- but this is too inefficient to maintain

-- OPERATIONS :

  procedure Initialize_Subsets ( n,max : in integer32 ) is

  -- The stack is allocated: we allocate all levels of the stack 
  -- and set the current level to 1, with no selections made.

  begin
    number_of_variables := n;
    maximum_selections := max;
    selections := new Standard_Integer_VecVecs.VecVec(1..n+1);
    for i in selections'range loop
      selections(i) := new Standard_Integer_Vectors.Vector(1..n);
    end loop;
    number_of_selections := new Standard_Integer_Vectors.Vector(1..n+1);
    current_level := 1;
    number_of_selections(current_level) := 0;
  end Initialize_Subsets;

  procedure Initialize_Iterator
              ( A : in Standard_Integer_Matrices.Matrix;
                max : in integer32 ) is

  -- The stack size can never become larger than the number of variables
  -- because each element that is put on the stack has exactly one more
  -- variable set to zero.  The number_of_selections(i) stores the index
  -- of the current monomial minus one because this index is increased
  -- each time we consider the top of the stack.

    stack_size : constant integer32 := A'last(2) + 1;

  begin
    number_of_variables := A'last(2);
    number_of_monomials := A'last(1);
    maximum_selections := max;
    incidence_matrix := new Standard_Integer_Matrices.Matrix'(A);
    selections := new Standard_Integer_VecVecs.VecVec(1..stack_size);
    nonzerovar := new Standard_Integer_VecVecs.VecVec(1..stack_size);
    for i in selections'range loop
      selections(i) := new Standard_Integer_Vectors.Vector'(A'range(2) => 0);
      nonzerovar(i) := new Standard_Integer_Vectors.Vector'(A'range(2) => 0);
    end loop;
    number_of_selections
      := new Standard_Integer_Vectors.Vector'(1..stack_size => 0);
    number_of_equations
      := new Standard_Integer_Vectors.Vector'(1..stack_size => 0);
    current_level := 1;
    number_of_selections(current_level) := 0;
  end Initialize_Iterator;

  function Set_to_Zero
             ( A : Standard_Integer_Matrices.Matrix; i : integer32;
               s : Standard_Integer_Vectors.Vector ) return boolean is
  begin
    for k in A'range(2) loop
      if A(i,k) > 0 and s(k) = 1
       then return true;
      end if;
    end loop;
    return false;
  end Set_to_Zero;

  function All_Set_to_Zero 
             ( A : Standard_Integer_Matrices.Matrix;
               s : Standard_Integer_Vectors.Vector ) return boolean is
  begin
    for i in A'range(1) loop
      if not Set_to_Zero(A,i,s)
       then return false;
      end if;
    end loop;
    return true;
  end All_Set_to_Zero;

  function Set_Nonzero
             ( A : Standard_Integer_Matrices.Matrix; i : integer32;
               s : Standard_Integer_Vectors.Vector ) return boolean is
  begin
    for k in A'range(2) loop
      if A(i,k) > 0 and s(k) = 0
       then return false;
      end if;
    end loop;
    return true; -- for all A(i,k) > 0 we had s(k) = 1
  end Set_Nonzero;

  procedure Update_Present_Variables
              ( s : in out Standard_Integer_Vectors.Vector;
                A : in Standard_Integer_Matrices.Matrix; i : in integer32 ) is
  begin
    for j in A'range(2) loop
      if A(i,j) > 0
       then s(j) := 1;
      end if;
    end loop;
  end Update_Present_Variables;

  procedure Enumerate ( A : in Standard_Integer_Matrices.Matrix;
                        s0_max : in integer32 ) is

    s0 : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    -- s0 contains the selected variables set to be zero
    s0_cnt : integer32 := 0; -- counts number of ones in s0
    s1 : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    -- s1 contains the variables which should not be set to zero
    --nq : constant integer32 := A'last(1)/2;
   -- eq : Standard_Integer_Vectors.Vector(1..nq) := (1..nq => 0);
    -- eq(k) = 1 means that k-th equation is skipped and remains
    eq_cnt : integer32 := 0; -- counts number of ones in eq
    skip : boolean;
    continue : boolean := true;

    procedure Enum ( k : in integer32 ) is

    -- DESCRIPTION :
    --   Examines the k-th monomial.

    begin
     -- put("k = "); put(k,1); 
     -- put("  s0 = "); put(s0);
     -- put("  s1 = "); put(s1); new_line;
      if k > A'last(1) then
       -- Report(s0,s0_cnt,s1,eq,eq_cnt,continue);
       -- if not Is_In(generated,s0) then
       --   Append(generated,generated_last,s0);
          Report(s0,s0_cnt,s1,eq_cnt,continue);
       -- end if;
      elsif Set_to_Zero(A,k,s0) then
        Enum(k+1);
      else
        skip := false;
        if k mod 2 = 1 then -- can we skip binomial ?
          if not Set_to_Zero(A,k+1,s0) 
           then skip := true; -- eq((k+1)/2) := 1;
                eq_cnt := eq_cnt + 1;
          end if;
        end if;
        if skip then
          declare
            s1back : constant Standard_Integer_Vectors.Vector(s1'range) := s1;
          begin
            Update_Present_Variables(s1,A,k);
            Update_Present_Variables(s1,A,k+1);
            Enum(k+2);
            s1 := s1back; -- eq((k+1)/2) := 0;
            eq_cnt := eq_cnt - 1;
          end;
        end if;
        if continue then
         -- if not All_Set_to_Zero(A,s0) then
            for j in A'range(2) loop
              if A(k,j) > 0 and s1(j) = 0 then
                if s0_cnt < s0_max then
                  s0(j) := 1; s0_cnt := s0_cnt + 1;
                  Enum(k+1);
                  s0(j) := 0; s0_cnt := s0_cnt - 1;
                end if;
              end if;
              exit when not continue;
            end loop;
         -- end if;
        end if;
      end if;
    end Enum;

  begin
    Enum(A'first(1));
   -- Clear(generated);
  end Enumerate;

  procedure Next_Subset
              ( setzero : out Standard_Integer_Vectors.Vector;
                cntzero : out integer32 ) is

    n : constant integer32 := number_of_variables;
    nbsl : integer32;
    accu : Standard_Integer_Vectors.Vector(1..n);

  begin
    if current_level < 1
     then cntzero := -1; return;
    end if;
    while current_level > 0 loop
      nbsl := number_of_selections(current_level);
      if nbsl = n then
        cntzero := 0;
        for i in selections(current_level)'range loop
          setzero(i) := selections(current_level)(i);
          if setzero(i) = 1
           then cntzero := cntzero + 1;
          end if;
        end loop;
        current_level := current_level - 1;  -- pop from the stack
        return;
      else
        cntzero := 0;
        for i in 1..nbsl loop
          accu(i) := selections(current_level)(i);
          if accu(i) = 1
           then cntzero := cntzero + 1;
          end if;
        end loop;
        nbsl := nbsl + 1;  -- consider next element
        -- replace current element on top of the stack
        number_of_selections(current_level) := nbsl;
        accu(nbsl) := 0;   -- do not select the variable
        for i in 1..nbsl loop
          selections(current_level)(i) := accu(i);
        end loop;
        -- add one element on top of the stack
        if cntzero < maximum_selections then
          current_level := current_level + 1;
          accu(nbsl) := 1;   -- do select the variable
          for i in 1..nbsl loop
            selections(current_level)(i) := accu(i);
          end loop;
          number_of_selections(current_level) := nbsl;
        end if;
      end if;
    end loop;
  end Next_Subset;

--  procedure Write_Stack is
--
--  -- DESCRIPTION :
--  --   Writes the current content of the stack,
--  --   useful to debug the Next_Selection procedure.
--
--  begin
--    for i in 1..current_level loop
--      put("s0("); put(i,1); put(") : ");
--      put(selections(i).all);
--      put(" | "); put(number_of_selections(i),1);
--      put(" | ");
--      put(nonzerovar(i).all); new_line;
--    end loop;
--  end Write_Stack;

  procedure Next_Selection
              ( setzero : out Standard_Integer_Vectors.Vector;
                cntzero : out integer32;
                nonzero : out Standard_Integer_Vectors.Vector;
                cntrest : out integer32 ) is

    nm : constant integer32 := number_of_monomials;
    nv : constant integer32 := number_of_variables;
    s0_accu : Standard_Integer_Vectors.Vector(1..nv);
    s1_accu : Standard_Integer_Vectors.Vector(1..nv);
    current_monomial,k,eq : integer32;
    skip : boolean;

    procedure Assign_Selection is

    -- DESCRIPTION :
    --   Assigns the top of the stack to the result parameters.

    begin
      cntzero := 0;
      for i in selections(current_level)'range loop
        setzero(i) := selections(current_level)(i);
        nonzero(i) := nonzerovar(current_level)(i);
        if setzero(i) = 1
         then cntzero := cntzero + 1;
        end if;
      end loop;
      cntrest := number_of_equations(current_level);
    end Assign_Selection;

    procedure Top_of_Stack is

    -- DESCRIPTION :
    --   Copies the top of the stack onto the work variables.

    begin
      cntzero := 0;
      for i in s0_accu'range loop
        s0_accu(i) := selections(current_level)(i);
        s1_accu(i) := nonzerovar(current_level)(i);
        if s0_accu(i) = 1
         then cntzero := cntzero + 1;
        end if;
      end loop;
      eq := number_of_equations(current_level);
    end Top_of_Stack;

    procedure Push_Combinations is

    -- DESCRIPTION :
    --   Pushes onto the stack the next combinations.

    begin
      if cntzero < maximum_selections then
        for j in reverse 1..nv loop -- reverse is same as recursive
          if incidence_matrix(current_monomial,j) > 0 then
            if s1_accu(j) = 0 then
              s0_accu(j) := 1;
              for i in s0_accu'range loop
                selections(current_level)(i) := s0_accu(i);
                nonzerovar(current_level)(i) := s1_accu(i);
              end loop;
              number_of_selections(current_level) := current_monomial;
              number_of_equations(current_level) := eq;
              current_level := current_level + 1;
              s0_accu(j) := 0;
            end if;
          end if;
        end loop;
      end if;
    end Push_Combinations;

    procedure Skip_Binomial is

    -- DESCRIPTION :
    --   Handles the case where we skip the current binomial equation.

    begin
      skip := false;
      if current_monomial mod 2 = 1 then -- can we skip binomial ?
        if not Set_to_Zero(incidence_matrix.all,current_monomial+1,s0_accu)
         then skip := true;
        end if;
      end if;
      if skip then
        k := current_monomial;
        Update_Present_Variables(s1_accu,incidence_matrix.all,k);
        Update_Present_Variables(s1_accu,incidence_matrix.all,k+1);
        current_level := current_level + 1;
        number_of_equations(current_level) := eq + 1;
        for i in s0_accu'range loop
          selections(current_level)(i) := s0_accu(i);
          nonzerovar(current_level)(i) := s1_accu(i);
        end loop;
        number_of_selections(current_level) := k + 1;
       -- put_line("The stack after skipping :"); Write_Stack;
      end if;
    end Skip_Binomial;

  begin
    if current_level < 1
     then cntzero := -1; cntrest := nm/2; return;
    end if;
    while current_level > 0 loop
      current_monomial := number_of_selections(current_level) + 1;
     -- put("current level : "); put(current_level,1);
     -- put("  current monomial : "); put(current_monomial,1);
     -- put_line("  the stack :"); Write_Stack;
      if current_monomial > nm then
       -- if not Is_In(generated,selections(current_level).all) then
       --   Append(generated,generated_last,selections(current_level).all);
          Assign_Selection;
          current_level := current_level - 1; -- pop stack
          return;
       -- else
       --   current_level := current_level - 1; -- pop stack and continue
       -- end if;
      else
        Top_of_Stack;
        if Set_to_Zero(incidence_matrix.all,current_monomial,s0_accu) then
          -- skipping the monomial, adjusting top of the stack
          -- notice that the current_monomial was already increased by 1
          number_of_selections(current_level) := current_monomial;
        elsif Set_Nonzero(incidence_matrix.all,current_monomial,s1_accu) then
          -- the recursive enumeration does not proceed once a monomial
          -- is prevented from being zero by the selection is s1
          current_level := current_level - 1; -- pop the stack
          Skip_Binomial;
        else
          Push_Combinations;
          current_level := current_level - 1; -- for the last + 1
          Skip_Binomial;
        end if;
      end if;
    end loop;
    cntzero := -1; -- if current_level became zero in while loop
  end Next_Selection;

 -- function Generated_Selections return List is
 -- begin
 --   return generated;
 -- end Generated_Selections;

  procedure Clear_Subsets is
  begin
    Standard_Integer_VecVecs.Deep_Clear(selections);
    Standard_Integer_Vectors.Clear(number_of_selections);
  end Clear_Subsets;

  procedure Clear_Iterator is
  begin
    Standard_Integer_Matrices.Clear(incidence_matrix);
    Standard_Integer_VecVecs.Deep_Clear(selections);
    Standard_Integer_VecVecs.Deep_Clear(nonzerovar);
    Standard_Integer_Vectors.Clear(number_of_selections);
    Standard_Integer_Vectors.Clear(number_of_equations);
   -- Clear(generated);
  end Clear_Iterator;

end Affine_Binomial_Iterator;
