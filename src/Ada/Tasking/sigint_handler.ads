with Ada.Interrupts.Names;
 
package Sigint_Handler is

  -- defines a simple handler for the SIGINT signal

  protected Handler is

    entry Wait; -- decreases the call_count

    procedure Handle; -- increases the call_count

    pragma Interrupt_Handler(Handle);
    pragma Attach_Handler(Handle, Ada.Interrupts.Names.Sigint);

  private
    call_count : natural := 0;
  end Handler;
 
end Sigint_Handler;
