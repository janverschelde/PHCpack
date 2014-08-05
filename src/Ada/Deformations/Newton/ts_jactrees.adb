with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with Standard_Derivative_Trees; 
with Standard_Jacobian_Trees;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Derivative_Trees;
with DoblDobl_Jacobian_Trees;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Derivative_Trees;
with QuadDobl_Jacobian_Trees;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Derivative_Trees; 
with Multprec_Jacobian_Trees;

procedure ts_jactrees is

-- DESCRIPTION :
--   Computes all derivatives of a polynomial system.

-- plain recursive enumeration without creation of a tree:

  procedure Standard_Show_all_Derivatives
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure DoblDobl_Show_all_Derivatives
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure QuadDobl_Show_all_Derivatives
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Multprec_Show_all_Derivatives
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys );

  procedure Standard_Show_all_Derivatives
              ( p : in Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    dp : Poly_Sys(1..n);

  begin
    for i in 1..n loop
      put("Derivative w.r.t. "); put(i,1); put(" : ");
      dp(i) := Diff(p,i);
      put(dp(i)); new_line;
    end loop;
    Standard_Show_all_Derivatives(dp);
  end Standard_Show_all_Derivatives;

  procedure DoblDobl_Show_all_Derivatives
              ( p : in DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    dp : Poly_Sys(1..n);

  begin
    for i in 1..n loop
      put("Derivative w.r.t. "); put(i,1); put(" : ");
      dp(i) := Diff(p,i);
      put(dp(i)); new_line;
    end loop;
    DoblDobl_Show_all_Derivatives(dp);
  end DoblDobl_Show_all_Derivatives;

  procedure QuadDobl_Show_all_Derivatives
              ( p : in QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    dp : Poly_Sys(1..n);

  begin
    for i in 1..n loop
      put("Derivative w.r.t. "); put(i,1); put(" : ");
      dp(i) := Diff(p,i);
      put(dp(i)); new_line;
    end loop;
    QuadDobl_Show_all_Derivatives(dp);
  end QuadDobl_Show_all_Derivatives;

  procedure Multprec_Show_all_Derivatives
              ( p : in Multprec_Complex_Polynomials.Poly ) is

    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    dp : Poly_Sys(1..n);

  begin
    for i in 1..n loop
      put("Derivative w.r.t. "); put(i,1); put(" : ");
      dp(i) := Diff(p,i);
      put(dp(i)); new_line;
    end loop;
    Multprec_Show_all_Derivatives(dp);
  end Multprec_Show_all_Derivatives;

  procedure Standard_Show_all_Derivatives
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        put("All derivatives of polynomial ");
        put(i,1); put_line(" : ");
        Standard_Show_all_Derivatives(p(i));
      end if;
    end loop;
  end Standard_Show_all_Derivatives;

  procedure DoblDobl_Show_all_Derivatives
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        put("All derivatives of polynomial ");
        put(i,1); put_line(" : ");
        DoblDobl_Show_all_Derivatives(p(i));
      end if;
    end loop;
  end DoblDobl_Show_all_Derivatives;

  procedure QuadDobl_Show_all_Derivatives
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        put("All derivatives of polynomial ");
        put(i,1); put_line(" : ");
        QuadDobl_Show_all_Derivatives(p(i));
      end if;
    end loop;
  end QuadDobl_Show_all_Derivatives;

  procedure Multprec_Show_all_Derivatives
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    use Multprec_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        put("All derivatives of polynomial ");
        put(i,1); put_line(" : ");
        Multprec_Show_all_Derivatives(p(i));
      end if;
    end loop;
  end Multprec_Show_all_Derivatives;

-- output of trees :

  procedure Write ( nd : in Standard_Derivative_Trees.Node ) is

    use Standard_Derivative_Trees;

  begin
    put("Polynomial "); put(nd.p); put_line(" has derivatives : ");
    for i in nd.d'range loop
      put("  w.r.t. "); put(i,1); put(" : ");
      if nd.d(i) = null
       then put_line("0");
       else put(nd.d(i).p); new_line;
      end if;
    end loop;
  end Write;

  procedure Write ( nd : in DoblDobl_Derivative_Trees.Node ) is

    use DoblDobl_Derivative_Trees;

  begin
    put("Polynomial "); put(nd.p); put_line(" has derivatives : ");
    for i in nd.d'range loop
      put("  w.r.t. "); put(i,1); put(" : ");
      if nd.d(i) = null
       then put_line("0");
       else put(nd.d(i).p); new_line;
      end if;
    end loop;
  end Write;

  procedure Write ( nd : in QuadDobl_Derivative_Trees.Node ) is

    use QuadDobl_Derivative_Trees;

  begin
    put("Polynomial "); put(nd.p); put_line(" has derivatives : ");
    for i in nd.d'range loop
      put("  w.r.t. "); put(i,1); put(" : ");
      if nd.d(i) = null
       then put_line("0");
       else put(nd.d(i).p); new_line;
      end if;
    end loop;
  end Write;

  procedure Write ( nd : in Multprec_Derivative_Trees.Node ) is

    use Multprec_Derivative_Trees;

  begin
    put("Polynomial "); put(nd.p); put_line(" has derivatives : ");
    for i in nd.d'range loop
      put("  w.r.t. "); put(i,1); put(" : ");
      if nd.d(i) = null
       then put_line("0");
       else put(nd.d(i).p); new_line;
      end if;
    end loop;
  end Write;

  procedure Write_Derivative_Nodes
              ( nd : in Standard_Derivative_Trees.Node ) is

    use Standard_Derivative_Trees;

    procedure Write ( nd : in Node; cont : out boolean ) is
    begin
      Write(nd);
      cont := true;
    end Write;
    procedure Write_all_Derivatives is new Enumerate(Write);

  begin
    Write_all_Derivatives(nd);
  end Write_Derivative_Nodes;

  procedure Write_Derivative_Nodes
              ( nd : in DoblDobl_Derivative_Trees.Node ) is

    use DoblDobl_Derivative_Trees;

    procedure Write ( nd : in Node; cont : out boolean ) is
    begin
      Write(nd);
      cont := true;
    end Write;
    procedure Write_all_Derivatives is new Enumerate(Write);

  begin
    Write_all_Derivatives(nd);
  end Write_Derivative_Nodes;

  procedure Write_Derivative_Nodes
              ( nd : in QuadDobl_Derivative_Trees.Node ) is

    use QuadDobl_Derivative_Trees;

    procedure Write ( nd : in Node; cont : out boolean ) is
    begin
      Write(nd);
      cont := true;
    end Write;
    procedure Write_all_Derivatives is new Enumerate(Write);

  begin
    Write_all_Derivatives(nd);
  end Write_Derivative_Nodes;

 -- procedure Write_Derivative_Nodes
 --             ( nd : in Multprec_Derivative_Trees.Node ) is

 --   use Multprec_Derivative_Trees;

 --   procedure Write ( nd : in Node; cont : out boolean ) is
 --   begin
 --     Write(nd);
 --     cont := true;
 --   end Write;
 --   procedure Write_all_Derivatives is new Enumerate(Write);

 -- begin
 --   Write_all_Derivatives(nd);
 -- end Write_Derivative_Nodes;

  procedure Write_Derivative_Nodes
              ( rnd : in Standard_Derivative_Trees.Array_of_Nodes ) is

    use Standard_Derivative_Trees;

    procedure Write ( nd : in Node; cont : out boolean ) is
    begin
      Write(nd);
      cont := true;
    end Write;
    procedure Write_all_Derivatives is new Enumerate_Nodes(Write);

  begin
    Write_all_Derivatives(rnd);
  end Write_Derivative_Nodes;

  procedure Write_Derivative_Nodes
              ( rnd : in DoblDobl_Derivative_Trees.Array_of_Nodes ) is

    use DoblDobl_Derivative_Trees;

    procedure Write ( nd : in Node; cont : out boolean ) is
    begin
      Write(nd);
      cont := true;
    end Write;
    procedure Write_all_Derivatives is new Enumerate_Nodes(Write);

  begin
    Write_all_Derivatives(rnd);
  end Write_Derivative_Nodes;

  procedure Write_Derivative_Nodes
              ( rnd : in QuadDobl_Derivative_Trees.Array_of_Nodes ) is

    use QuadDobl_Derivative_Trees;

    procedure Write ( nd : in Node; cont : out boolean ) is
    begin
      Write(nd);
      cont := true;
    end Write;
    procedure Write_all_Derivatives is new Enumerate_Nodes(Write);

  begin
    Write_all_Derivatives(rnd);
  end Write_Derivative_Nodes;

  procedure Write_Derivative_Nodes
              ( rnd : in Multprec_Derivative_Trees.Array_of_Nodes ) is

    use Multprec_Derivative_Trees;

    procedure Write ( nd : in Node; cont : out boolean ) is
    begin
      Write(nd);
      cont := true;
    end Write;
    procedure Write_all_Derivatives is new Enumerate_Nodes(Write);

  begin
    Write_all_Derivatives(rnd);
  end Write_Derivative_Nodes;

  procedure Write ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        put("  "); put(jm(i,j));
      end loop;
      new_line;
    end loop;
  end Write;

  procedure Write ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        put("  "); put(jm(i,j));
      end loop;
      new_line;
    end loop;
  end Write;

  procedure Write ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        put("  "); put(jm(i,j));
      end loop;
      new_line;
    end loop;
  end Write;

  procedure Write ( jm : in Multprec_Complex_Jaco_Matrices.Jaco_Mat ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        put("  "); put(jm(i,j));
      end loop;
      new_line;
    end loop;
  end Write;

  procedure Write_Jacobian_Nodes
              ( nd : in Standard_Jacobian_Trees.Node ) is

    use Standard_Jacobian_Trees;

    procedure Write ( nd : Node; cont : out boolean ) is
    begin
      put_line("Jacobian matrix"); Write(nd.a.all); 
      put_line("has derivatives");
      for i in nd.d'range loop
        put("  w.r.t. "); put(i,1); put(" : ");
        if nd.d(i) = null
         then put_line("0");
         else new_line; Write(nd.d(i).a.all);
        end if;
      end loop;
      cont := true;
    end Write;
    procedure Write_all_Jacobian_Matrices is new Enumerate_Nodes(Write);

  begin
    Write_all_Jacobian_Matrices(nd);
  end Write_Jacobian_Nodes;

  procedure Write_Jacobian_Nodes
              ( nd : in DoblDobl_Jacobian_Trees.Node ) is

    use DoblDobl_Jacobian_Trees;

    procedure Write ( nd : Node; cont : out boolean ) is
    begin
      put_line("Jacobian matrix"); Write(nd.a.all); 
      put_line("has derivatives");
      for i in nd.d'range loop
        put("  w.r.t. "); put(i,1); put(" : ");
        if nd.d(i) = null
         then put_line("0");
         else new_line; Write(nd.d(i).a.all);
        end if;
      end loop;
      cont := true;
    end Write;
    procedure Write_all_Jacobian_Matrices is new Enumerate_Nodes(Write);

  begin
    Write_all_Jacobian_Matrices(nd);
  end Write_Jacobian_Nodes;

  procedure Write_Jacobian_Nodes
              ( nd : in QuadDobl_Jacobian_Trees.Node ) is

    use QuadDobl_Jacobian_Trees;

    procedure Write ( nd : Node; cont : out boolean ) is
    begin
      put_line("Jacobian matrix"); Write(nd.a.all); 
      put_line("has derivatives");
      for i in nd.d'range loop
        put("  w.r.t. "); put(i,1); put(" : ");
        if nd.d(i) = null
         then put_line("0");
         else new_line; Write(nd.d(i).a.all);
        end if;
      end loop;
      cont := true;
    end Write;
    procedure Write_all_Jacobian_Matrices is new Enumerate_Nodes(Write);

  begin
    Write_all_Jacobian_Matrices(nd);
  end Write_Jacobian_Nodes;

  procedure Write_Jacobian_Nodes
              ( nd : in Multprec_Jacobian_Trees.Node ) is

    use Multprec_Jacobian_Trees;

    procedure Write ( nd : Node; cont : out boolean ) is
    begin
      put_line("Jacobian matrix"); Write(nd.a.all); 
      put_line("has derivatives");
      for i in nd.d'range loop
        put("  w.r.t. "); put(i,1); put(" : ");
        if nd.d(i) = null
         then put_line("0");
         else new_line; Write(nd.d(i).a.all);
        end if;
      end loop;
      cont := true;
    end Write;
    procedure Write_all_Jacobian_Matrices is new Enumerate_Nodes(Write);

  begin
    Write_all_Jacobian_Matrices(nd);
  end Write_Jacobian_Nodes;

  procedure Standard_Enumerate_all_Derivatives
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Derivative_Trees;

    dp : constant Array_of_Nodes(p'range) := Create(p);

  begin
    Write_Derivative_Nodes(dp);
  end Standard_Enumerate_all_Derivatives;

  procedure DoblDobl_Enumerate_all_Derivatives
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Derivative_Trees;

    dp : constant Array_of_Nodes(p'range) := Create(p);

  begin
    Write_Derivative_Nodes(dp);
  end DoblDobl_Enumerate_all_Derivatives;

  procedure QuadDobl_Enumerate_all_Derivatives
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Derivative_Trees;

    dp : constant Array_of_Nodes(p'range) := Create(p);

  begin
    Write_Derivative_Nodes(dp);
  end QuadDobl_Enumerate_all_Derivatives;

  procedure Multprec_Enumerate_all_Derivatives
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    use Multprec_Derivative_Trees;

    dp : constant Array_of_Nodes(p'range) := Create(p);

  begin
    Write_Derivative_Nodes(dp);
  end Multprec_Enumerate_all_Derivatives;

  procedure Standard_Enumerate_all_Jacobian_Matrices
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range(1),1..n) := Create(p);
    dp : constant Node(n) := Create(jm);

  begin
    Write_Jacobian_Nodes(dp);
  end Standard_Enumerate_all_Jacobian_Matrices;

  procedure DoblDobl_Enumerate_all_Jacobian_Matrices
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range(1),1..n) := Create(p);
    dp : constant Node(n) := Create(jm);

  begin
    Write_Jacobian_Nodes(dp);
  end DoblDobl_Enumerate_all_Jacobian_Matrices;

  procedure QuadDobl_Enumerate_all_Jacobian_Matrices
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range(1),1..n) := Create(p);
    dp : constant Node(n) := Create(jm);

  begin
    Write_Jacobian_Nodes(dp);
  end QuadDobl_Enumerate_all_Jacobian_Matrices;

  procedure Multprec_Enumerate_all_Jacobian_Matrices
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range(1),1..n) := Create(p);
    dp : constant Node(n) := Create(jm);

  begin
    Write_Jacobian_Nodes(dp);
  end Multprec_Enumerate_all_Jacobian_Matrices;

-- selective derivative calculation :

  procedure Read_Vector ( v : in out Vector )is
  begin
    for i in v'range loop
      put("How many times to differentiate w.r.t. ");
      put(i,1); put(" ? ");
      get(v(i));
    end loop;
  end Read_Vector;

  procedure Interactive_Derivatives
              ( p : in Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;
    use Standard_Derivative_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    dp : Node(n) := Initialize(p);
    v : Vector(1..n);
    d : Poly;

  begin
    Read_Vector(v);
    Derivative(dp,v,d);
    put("The derivative of "); put(dp.p);
    put(" w.r.t."); put(v); put(" is "); put(d); new_line;
    Clear(d);
    put_line("The current tree of derivatives : ");
    Write_Derivative_Nodes(dp);
  end Interactive_Derivatives;

  procedure Standard_Interactive_Jacobians
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range,1..n) := Create(p);
    nd : Node(n) := Initialize(jm);
    v : Vector(1..n);
    da : Link_to_Jaco_Mat;
    ans : character;

  begin
    loop
      Read_Vector(v);
      Derivative(nd,v,da);
      put_line("The derivative of "); Write(nd.a.all);
      put(" w.r.t. "); put(v); put(" is ");
      if da = null
       then put_line("0");
       else new_line; Write(da.all);
      end if;
      put_line("The current tree of derivatives : ");
      Write_Jacobian_Nodes(nd);
      put("Do you wish to continue ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end Standard_Interactive_Jacobians;

  procedure DoblDobl_Interactive_Jacobians
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range,1..n) := Create(p);
    nd : Node(n) := Initialize(jm);
    v : Vector(1..n);
    da : Link_to_Jaco_Mat;
    ans : character;

  begin
    loop
      Read_Vector(v);
      Derivative(nd,v,da);
      put_line("The derivative of "); Write(nd.a.all);
      put(" w.r.t. "); put(v); put(" is ");
      if da = null
       then put_line("0");
       else new_line; Write(da.all);
      end if;
      put_line("The current tree of derivatives : ");
      Write_Jacobian_Nodes(nd);
      put("Do you wish to continue ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end DoblDobl_Interactive_Jacobians;

  procedure QuadDobl_Interactive_Jacobians
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range,1..n) := Create(p);
    nd : Node(n) := Initialize(jm);
    v : Vector(1..n);
    da : Link_to_Jaco_Mat;
    ans : character;

  begin
    loop
      Read_Vector(v);
      Derivative(nd,v,da);
      put_line("The derivative of "); Write(nd.a.all);
      put(" w.r.t. "); put(v); put(" is ");
      if da = null
       then put_line("0");
       else new_line; Write(da.all);
      end if;
      put_line("The current tree of derivatives : ");
      Write_Jacobian_Nodes(nd);
      put("Do you wish to continue ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end QuadDobl_Interactive_Jacobians;

  procedure Multprec_Interactive_Jacobians
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range,1..n) := Create(p);
    nd : Node(n) := Initialize(jm);
    v : Vector(1..n);
    da : Link_to_Jaco_Mat;
    ans : character;

  begin
    loop
      Read_Vector(v);
      Derivative(nd,v,da);
      put_line("The derivative of "); Write(nd.a.all);
      put(" w.r.t. "); put(v); put(" is ");
      if da = null
       then put_line("0");
       else new_line; Write(da.all);
      end if;
      put_line("The current tree of derivatives : ");
      Write_Jacobian_Nodes(nd);
      put("Do you wish to continue ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end Multprec_Interactive_Jacobians;

  procedure Main is

    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    stlp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    ddlp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    qdlp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    mplp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    ans : character;

  begin
    new_line;
    put_line("Testing the calculation of trees of derivatives ...");
    new_line;
    put_line("MENU to derive one single polynomial :");
    put_line("  0. show all derivatives in recursive calculation;");
    put_line("  1. interactive selective derivative calculation;");
    put_line("MENU to derive Jacobian Matrices, standard coefficients :");
    put_line("  2. show all derivatives in recursive calculation;");
    put_line("  3. complete tree calculation of system and printing;");
    put_line("  4. complete tree calculation of all Jacobian matrices;");
    put_line("  5. selective calculation of Jacobian matrices.");
    put_line("MENU to derive Jacobian Matrices, dobldobl coefficients :");
    put_line("  6. show all derivatives in recursive calculation;");
    put_line("  7. complete tree calculation of system and printing;");
    put_line("  8. complete tree calculation of all Jacobian matrices;");
    put_line("  9. selective calculation of Jacobian matrices.");
    put_line("MENU to derive Jacobian Matrices, quaddobl coefficients :");
    put_line("  A. show all derivatives in recursive calculation;");
    put_line("  B. complete tree calculation of system and printing;");
    put_line("  C. complete tree calculation of all Jacobian matrices;");
    put_line("  D. selective calculation of Jacobian matrices.");
    put_line("MENU to derive Jacobian Matrices, multiprecision coefficients :");
    put_line("  E. show all derivatives in recursive calculation;");
    put_line("  F. complete tree calculation of system and printing;");
    put_line("  G. complete tree calculation of all Jacobian matrices;");
    put_line("  H. selective calculation of Jacobian matrices.");
    put("Type 0, 1, .., 9, or A, B, .., H to select : ");
    Ask_Alternative(ans,"0123456789ABCDEFGH");
    new_line;
    if ans = '0' or ans = '1' then
      put("Give the number of variables : "); get(n);
      Symbol_Table.Init(n);
      put("Give p : "); get(p);
      put("Your p : "); put(p);
      new_line;
      case ans is
        when '0' => Standard_Show_All_Derivatives(p);
        when '1' => Interactive_Derivatives(p);
        when others => null;
      end case;
    else
      case ans is
        when '2' | '3' | '4' | '5' => get(stlp);
        when '6' | '7' | '8' | '9' => get(ddlp);
        when 'A' | 'B' | 'C' | 'D' => get(qdlp);
        when 'E' | 'F' | 'G' | 'H' => get(mplp);
        when others => null;
      end case;
      case ans is
        when '2' => Standard_Show_all_Derivatives(stlp.all);
        when '3' => Standard_Enumerate_all_Derivatives(stlp.all);
        when '4' => Standard_Enumerate_all_Jacobian_Matrices(stlp.all);
        when '5' => Standard_Interactive_Jacobians(stlp.all);
        when '6' => DoblDobl_Show_all_Derivatives(ddlp.all);
        when '7' => DoblDobl_Enumerate_all_Derivatives(ddlp.all);
        when '8' => DoblDobl_Enumerate_all_Jacobian_Matrices(ddlp.all);
        when '9' => DoblDobl_Interactive_Jacobians(ddlp.all);
        when 'A' => QuadDobl_Show_all_Derivatives(qdlp.all);
        when 'B' => QuadDobl_Enumerate_all_Derivatives(qdlp.all);
        when 'C' => QuadDobl_Enumerate_all_Jacobian_Matrices(qdlp.all);
        when 'D' => QuadDobl_Interactive_Jacobians(qdlp.all);
        when 'E' => Multprec_Show_all_Derivatives(mplp.all);
        when 'F' => Multprec_Enumerate_all_Derivatives(mplp.all);
        when 'G' => Multprec_Enumerate_all_Jacobian_Matrices(mplp.all);
        when 'H' => Multprec_Interactive_Jacobians(mplp.all);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_jactrees;
