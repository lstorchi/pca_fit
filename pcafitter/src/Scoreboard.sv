////////////////////////////////////////////////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////s           www.testbench.in           s////
////s                                      s////
////s        SystemVerilog Tutorial        s////
////s           gopi@testbench.in          s////
////s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s////
////////////////////////////////////////////////
`ifndef GUARD_SCOREBOARD
`define GUARD_SCOREBOARD

`include "Packet.sv"
`include "Coverage.sv"

class Scoreboard;

	//define the mailbox for driver and receiver
	mailbox drvr2sb;
	mailbox rcvr2sb;
	
	//instantiate the coverage function
	coverage cov = new();

	function new(mailbox drvr2sb,mailbox rcvr2sb);
	  this.drvr2sb = drvr2sb;
	  this.rcvr2sb = rcvr2sb;
	endfunction:new


	task start();
	  packet pkt_rcv,pkt_exp;
	  forever
	  begin
		rcvr2sb.get(pkt_rcv);
		$display(" 0 : Scorebooard : Scoreboard received a packet from receiver ",$time);
		drvr2sb.get(pkt_exp);
		if(pkt_rcv.compare(pkt_exp))
		begin
			$display(" 0 : Scoreboardd :Packet Matched ",$time);
			cov.sample(pkt_exp);
		end
		else
		  error++;
	  end
	endtask : start

endclass

`endif