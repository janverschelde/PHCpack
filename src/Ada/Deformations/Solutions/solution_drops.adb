with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Solution_Drops is

  function Drop ( s : Standard_Complex_Solutions.Solution; k : natural32 )
                return Standard_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : Standard_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : DoblDobl_Complex_Solutions.Solution; k : natural32 )
                return DoblDobl_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : DoblDobl_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : TripDobl_Complex_Solutions.Solution; k : natural32 )
                return TripDobl_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : TripDobl_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : QuadDobl_Complex_Solutions.Solution; k : natural32 )
                return QuadDobl_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : QuadDobl_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : PentDobl_Complex_Solutions.Solution; k : natural32 )
                return PentDobl_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : PentDobl_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : OctoDobl_Complex_Solutions.Solution; k : natural32 )
                return OctoDobl_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : OctoDobl_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : DecaDobl_Complex_Solutions.Solution; k : natural32 )
                return DecaDobl_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : DecaDobl_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : HexaDobl_Complex_Solutions.Solution; k : natural32 )
                return HexaDobl_Complex_Solutions.Solution is
  begin
    if s.n <= 1 or k < 1 or k > natural32(s.n) then
      return s;
    else
      declare
        res : HexaDobl_Complex_Solutions.Solution(s.n-1);
        ind : constant integer32 := integer32(k);
      begin
        res.t := s.t;
        res.m := s.m;
        for i in 1..(ind-1) loop
          res.v(i) := s.v(i);
        end loop;
        for i in (ind+1)..s.n loop
          res.v(i-1) := s.v(i);
        end loop;
        res.err := s.err;
        res.rco := s.rco;
        res.res := s.res;
        return res;
      end;
    end if;
  end Drop;

  function Drop ( s : Standard_Complex_Solutions.Solution_List;
                  k : natural32 )
                return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

  function Drop ( s : DoblDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

  function Drop ( s : TripDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return TripDobl_Complex_Solutions.Solution_List is

    use TripDobl_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

  function Drop ( s : QuadDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

  function Drop ( s : PentDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return PentDobl_Complex_Solutions.Solution_List is

    use PentDobl_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

  function Drop ( s : OctoDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return OctoDobl_Complex_Solutions.Solution_List is

    use OctoDobl_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

  function Drop ( s : DecaDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return DecaDobl_Complex_Solutions.Solution_List is

    use DecaDobl_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

  function Drop ( s : HexaDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return HexaDobl_Complex_Solutions.Solution_List is

    use HexaDobl_Complex_Solutions;

    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Drop(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Drop;

end Solution_Drops;
