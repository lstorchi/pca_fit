////////////////////////////////////////////////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////s           www.testbench.in           s////
////s                                      s////
////s        SystemVerilog Tutorial        s////
////s           gopi@testbench.in          s////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////////////////////////////////////////////////
`ifndef GUARD_TOP
`define GUARD_TOP

`include "interface.sv"

module top();

/////////////////////////////////////////////////////
// Clock Declaration and Generation                //
/////////////////////////////////////////////////////
bit Clock;

initial
  forever #10 Clock = ~Clock;
  
///////////////////////////////
// Interface global signal	///
///////////////////////////////
gloabl_interface global_intf(Clock);
  
/////////////////////////////////////////////////////
//  Memory interface instance                      //
/////////////////////////////////////////////////////

mem_interface mem_intf(Clock);

/////////////////////////////////////////////////////
//  Input interface instance                       //
/////////////////////////////////////////////////////

input_interface input_intf(Clock);

/////////////////////////////////////////////////////
//  output interface instance                      //
/////////////////////////////////////////////////////
//Define the output interface to the RTL module
output_interface output_intf(Clock);

/////////////////////////////////////////////////////
//  Program block Testcase instance                //
/////////////////////////////////////////////////////

testcase TC (mem_intf,input_intf,output_intf, global_intf);

/////////////////////////////////////////////////////
//  DUT instance and signal connection             //
/////////////////////////////////////////////////////

	TrackFitterRTL 	DUT    
				(.clk(Clock),
				.reset(global_intf.reset),
				.dv_in_0(input_intf.data_valid),
				.data_in_x_0(input_intf.data_in_x),			
				.data_in_y_0(input_intf.data_in_y),
				.data_in_z_0(input_intf.data_in_z),
				.dv_out(output_intf.ready),
				.data_out(output_intf.data_out),
				.mem_en(mem_intf.mem_en),
				.mem_rd_wr(mem_intf.mem_rd_wr),
				.mem_add(mem_intf.mem_add),
				.mem_data(mem_intf.mem_data)
	);


endmodule


`endif


























