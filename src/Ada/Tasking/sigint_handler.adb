package body Sigint_Handler is
 
  protected body Handler is
 
    entry Wait when call_count > 0 is
    begin
      call_count := call_count - 1;
    end Wait;
 
    procedure Handle is
    begin
      call_count := call_count + 1;
    end Handle;
 
  end Handler;
 
end Sigint_Handler;
