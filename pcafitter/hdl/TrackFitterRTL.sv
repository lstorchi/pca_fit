 
 module TrackFitterRTL (
				clk,
				reset,
				dv_in_0,
				data_in_x_0,			
				data_in_y_0,
				data_in_z_0,
				dv_out,
				data_out,
				mem_en, mem_rd_wr, mem_add, mem_data
				);
 input clk; 
 input reset;
 input dv_in_0;
 input [7:0] data_in_x_0;
 input [7:0] data_in_y_0;
 input [7:0] data_in_z_0;


 input mem_en; 
 input mem_rd_wr;
 input [1:0] mem_add;
 input  [7:0] mem_data;

 output         dv_out; 
 output  [7:0]   data_out; 

	 
endmodule //router

