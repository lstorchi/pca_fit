////////////////////////////////////////////////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////s           www.testbench.in           s////
////s                                      s////
////s        SystemVerilog Tutorial        s////
////s           gopi@testbench.in          s////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////////////////////////////////////////////////
`ifndef GUARD_INTERFACE
`define GUARD_INTERFACE

///////////////////////////////
// Interface global signal	///
///////////////////////////////
interface gloabl_interface(input bit clock);
	logic           reset; 
  
	modport IG(output reset,input clock);
  
endinterface

//////////////////////////////////////////
// Interface declaration for the memory///
//////////////////////////////////////////

interface mem_interface(input bit clock);
  logic [7:0] mem_data;
  logic [1:0] mem_add;
  logic       mem_en;
  logic       mem_rd_wr;
  
  clocking cb@(posedge clock);
     default input #1 output #1;		//define the input and output skew
     output mem_data;
     output mem_add;
     output mem_en;
     output mem_rd_wr;
  endclocking
  
  modport MEM(clocking cb,input clock);

endinterface

////////////////////////////////////////////
// Interface for the input side of switch.//
////////////////////////////////////////////
//define the interface with the information to each input layer
interface input_interface(input bit clock);
	logic           data_valid;			//data valid signal
	logic     [7:0] data_in_x;			//information of stub in each layer
	logic     [7:0] data_in_y;
	logic     [7:0] data_in_z;

	clocking cb@(posedge clock);
		default input #1 output #1;
		output    data_valid;
		output    data_in_x;
		output    data_in_y;
		output    data_in_z;
	endclocking
  
	modport IP(clocking cb,input clock);
  
endinterface

/////////////////////////////////////////////////
// Interface for the output side of the switch.//
// output_interface is for only one output port//
/////////////////////////////////////////////////

interface output_interface(input bit clock);
  logic    [7:0] data_out;
  logic    ready;
  
  clocking cb@(posedge clock);
    default input #1 output #1;
    input     data_out;
    input     ready;
  endclocking
  
  modport OP(clocking cb,input clock);

endinterface


//////////////////////////////////////////////////

`endif 
