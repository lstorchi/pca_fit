////////////////////////////////////////////////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////s           www.testbench.in           s////
////s                                      s////
////s        SystemVerilog Tutorial        s////
////s                                      s////
////s           gopi@testbench.in          s////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////////////////////////////////////////////////
`ifndef GUARD_ENV
`define GUARD_ENV

`include "Driver.sv"
`include "Receiver.sv"
`include "Scoreboard.sv"

class environment ;

	//declare virtual interface considering the modport declaration
	virtual mem_interface.MEM    	mem_intf;
	virtual input_interface.IP  	input_intf;
	virtual output_interface.OP 	output_intf;
	virtual gloabl_interface.IG		global_intf;
	
	//define the Driver class and the mailbox
	Driver drvr;
	mailbox drvr2sb;
	
	//define the receiver
	Receiver rcvr;
	mailbox rcvr2sb;
	
	//define the scoreboard
	Scoreboard sb;

function new(virtual mem_interface.MEM    mem_intf_new,
             virtual input_interface.IP  input_intf_new,
             virtual output_interface.OP output_intf_new,
			 virtual gloabl_interface.IG global_intf_new 
			 );

  this.mem_intf      = mem_intf_new;
  this.input_intf    = input_intf_new;
  this.output_intf   = output_intf_new;
  this.global_intf	= global_intf_new;
  
  $display(" %0d : Environemnt : created env object",$time);
endfunction : new

function void build();
	$display(" %0d : Environemnt : start of build() method",$time);
	
	drvr2sb = new();
	drvr= new(input_intf, drvr2sb);
	rcvr2sb = new();
	rcvr= new(output_intf,rcvr2sb);
	
	sb = new(drvr2sb,rcvr2sb);
	
	$display(" %0d : Environemnt : end of build() method",$time);
endfunction :build

//method to reset the signals
task reset();
	$display(" %0d : Environemnt : start of reset() method",$time);
	
	//assign reset to signal
	mem_intf.cb.mem_data <= 0; 
	mem_intf.cb.mem_add <= 0; 
	mem_intf.cb.mem_en <= 0; 
	mem_intf.cb.mem_rd_wr <= 0; 
	input_intf.cb.data_valid <= 0; 
	input_intf.cb.data_in_x <= 0; 
	input_intf.cb.data_in_y <= 0; 
	input_intf.cb.data_in_z <= 0; 
	
	// Reset the DUT 
	global_intf.reset <= 0;
	repeat (10) @ global_intf.clock; 
	global_intf.reset <= 1; 
	repeat (4) @ global_intf.clock; 
	global_intf.reset <= 0;
	
	$display(" %0d : Environemnt : end of reset() method",$time);
endtask : reset

//in the configuration DUT we store the constant in the memory
task cfg_dut();
	$display(" %0d : Environemnt : start of cfg_dut() method",$time);
	
	mem_intf.cb.mem_en <= 1;
  @(posedge mem_intf.clock);
  mem_intf.cb.mem_rd_wr <= 1;

  @(posedge mem_intf.clock);
  mem_intf.cb.mem_add  <= 8'h0;
  mem_intf.cb.mem_data <= `P0;
  $display(" 0 : Environemnt : Port 0 Address %h ",$time,`P0);

  @(posedge mem_intf.clock);
  mem_intf.cb.mem_add  <= 8'h1;
  mem_intf.cb.mem_data <= `P1;
  $display(" 0 : Environemnt : Port 1 Address %h ",$time,`P1);

  @(posedge mem_intf.clock);
  mem_intf.cb.mem_add  <= 8'h2;
  mem_intf.cb.mem_data <= `P2;
  $display(" 0 : Environemnt : Port 2 Address %h ",$time,`P2);

  @(posedge mem_intf.clock);
  mem_intf.cb.mem_add  <= 8'h3;
  mem_intf.cb.mem_data <= `P3;
  $display(" 0 : Environemnt : Port 3 Address %h ",$time,`P3);

  @(posedge mem_intf.clock);
  mem_intf.cb.mem_en    <=0;
  mem_intf.cb.mem_rd_wr <= 0;
  mem_intf.cb.mem_add   <= 0;
  mem_intf.cb.mem_data  <= 0;
 
	$display(" %0d : Environemnt : end of cfg_dut() method",$time);
endtask : cfg_dut

task start();
	$display(" %0d : Environemnt : start of start() method",$time);
	
	fork
    drvr.start();
    rcvr.start();
    sb.start();
	join_any
	
	$display(" %0d : Environemnt : end of start() method",$time);
endtask : start

task wait_for_end();
 $display(" %0d : Environemnt : start of wait_for_end() method",$time);
 $display(" %0d : Environemnt : end of wait_for_end() method",$time);
endtask : wait_for_end

//method use to run all the function in the environment class
task run();
 $display(" %0d : Environemnt : start of run() method",$time);
 build();
 reset();
 cfg_dut();
 start();
 wait_for_end();
 report();
 $display(" %0d : Environemnt : end of run() method",$time);
endtask : run

task report();
$display("\n\n*************************************************");
   if( 0 == error)
       $display("********            TEST PASSED         *********");
   else
       $display("********    TEST Failed with 0 errors *********",error);

   $display("*************************************************\n\n");
endtask : report
       
endclass

`endif
