////////////////////////////////////////////////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////s           www.testbench.in           s////
////s                                      s////
////s        SystemVerilog Tutorial        s////
////s           gopi@testbench.in          s////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////////////////////////////////////////////////
`ifndef GUARD_TESTCASE
`define GUARD_TESTCASE

`include "environment.sv"

class small_packet extends packet;

	constraint small_c { data.size > 200 ; }

endclass

program testcase(mem_interface.MEM mem_intf,input_interface.IP input_intf,output_interface.OP output_intf, gloabl_interface.IG global_intf);

environment env;

small_packet spkt;


initial
begin
$display(" ******************* Start of testcase ****************");
//env = new(mem_intf,input_intf,output_intf, global_intf);
//env.run();
spkt = new();
env = new(mem_intf,input_intf,output_intf, global_intf);
env.build();
env.drvr.gpkt = spkt;
env.reset();
env.cfg_dut();
env.start();
env.wait_for_end();
env.report();



#1000;
end

final
$display(" ******************** End of testcase *****************");


endprogram 

`endif
