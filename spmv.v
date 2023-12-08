// -----------------------------------------------------------------------------
// Engineer       : Zhang Haoyang
// Design Name    : s_spmv
// Module Name    : s_spmv
// Target Devices : 
// Tool Versions  : 
// Description    : 
//    sparse matrix-vector multiplication PE
//    B=A @ x
//    with custom compression format
//    
//    
//    
//    
// Revision       :
// Version        Date        Author        Descriptin
// 0              2023/12/07  HY Zhang    
// 
// 
// 
// -----------------------------------------------------------------------------
module  spmv
module  spmv
  #(
    parameter       NNZ             =   'd100000  ,       //
    parameter       WIDTH_NNZ       =   'd4       ,       //
    parameter       WIDTH_value     =   'd16       ,       //
    parameter       WIDTH_col       =   'd16       ,       //
    parameter       WIDTH_ADDR      =   'd4        ,
    parameter       WIDTH_ADDR_VEC  =   'd16          
  )   

(   output wire     [WIDTH_value-1:0]   final_result                  , 
    input                               clk                           ,
    input                               rst_n     
);      


wire        [WIDTH_col-1:0]             col_index_fifo                ;        
wire                                    EOR_fifo                      ;
//----------------------------------------------
//                control
//---------------------------------------------
reg         [WIDTH_NNZ:0]               cnt_en                        ;
wire                                    rd_en                         ;
assign rd_en  =  ((cnt_en>=1) & (cnt_en<=NNZ+4))  ?   1'b1   :  1'b0  ;
reg                                     REG_rd_en_val                 ;
wire                                    wire_rd_en_val                ;
assign wire_rd_en_val       =           REG_rd_en_val                 ;
reg                                     eor_1d                        ;
reg                                     eor_2d                        ;
reg                                     eor_3d                        ;
reg                                     eor_4d                        ;
reg                                     next_row_sig                  ;
reg                                     next_row_sig_1d               ;
reg         [WIDTH_ADDR_VEC-1:0]        row_count                     ;//WIDTH: use row_count to fetch vector value
reg         [WIDTH_ADDR-1:0]            ADDR_col                      ;
wire        [WIDTH_ADDR_VEC-1:0]        wire_row_count                ;
assign wire_row_count       =           row_count                     ;
wire        [WIDTH_value-1:0]           x_lower                       ;//-------obtain x
wire        [WIDTH_value-1:0]           x_upper                       ;//-------obtain x
reg         [WIDTH_value-1:0]           partial_upper                 ;
reg                                     diag_signal                   ;
reg                                     diag_signal_1d                ;
reg                                     diag_signal_2d                ;
reg         [WIDTH_value-1:0]           partial_lower                 ;
wire        [WIDTH_value-1:0]           current_lower                 ;
wire        [WIDTH_value-1:0]           write_through_lower            ;
reg         [WIDTH_col-1:0]             partial_sum_up                ;
reg         [WIDTH_col-1:0]             partial_sum_low               ;
reg         [WIDTH_col-1:0]             col_index_1d                  ;
reg         [WIDTH_col-1:0]             col_index_2d                  ;
reg         [WIDTH_col-1:0]             col_index_3d                  ;
reg         [WIDTH_col-1:0]             col_index_4d                  ;
wire        [WIDTH_col-1:0]             ADDR_low_sum_read             ;
assign ADDR_low_sum_read    =           col_index_1d                  ;
wire        [WIDTH_col-1:0]             ADDR_low_sum_write            ;
assign ADDR_low_sum_write   =           col_index_3d                  ;
wire                                    same_row_sig                  ;
assign same_row_sig  =   (col_index_2d==col_index_3d) ?  1'b1 :1'b0   ;//element in upper part have the same col, means they are in the same row of lower part         
wire                                    rd_and_wr_conflict            ;//debug indicator:concurrent rd and wr
assign  rd_and_wr_conflict = (col_index_2d==col_index_4d)             ;
reg         [WIDTH_value-1:0]           final_all_output              ;
assign final_result =                   final_all_output              ;

//--------------------------
always @(posedge clk or negedge rst_n) begin
    if(!rst_n)
        cnt_en          <=          10'd0                 ;
    else
        cnt_en          <=          cnt_en+10'd1          ;  
end

//---------control-------rden delay-----------1 delay-------for val and eor

always @(posedge clk or negedge rst_n) begin
  if(!rst_n)
    REG_rd_en_val           <=        1'b0                   ;
  else
    REG_rd_en_val           <=        rd_en                  ;
end 

//---------control------------EOR delay----------

always @(posedge clk or negedge rst_n) begin
  if (!rst_n ) begin
    eor_1d         <=          'b0            ;
    eor_2d          <=         'b0            ;
    eor_3d         <=          'b0            ;
    eor_4d          <=         'b0            ;
  end
  else begin
    eor_1d          <=          EOR_fifo      ;
    eor_2d          <=          eor_1d        ;
    eor_3d          <=          eor_2d        ;
    eor_4d          <=          eor_3d        ;
  end  
end

always @(posedge clk or negedge rst_n) begin
  if (!rst_n ) 
    next_row_sig          <=          'b0               ;
  else if (rd_en & cnt_en>= 'd4) begin
    if (  eor_2d ==1'b0) 
      next_row_sig        <=          'b1               ;
    else  
      next_row_sig        <=          'b0               ;
  end
  else
    next_row_sig          <=          'b0               ;
end 

always @(posedge clk or negedge rst_n) begin
  if (!rst_n ) 
    next_row_sig_1d          <=          'b0            ;
  else if (rd_en & cnt_en>= 'd5) begin
    if (  eor_3d ==1'b0) 
      next_row_sig_1d        <=          'b1            ;
    else  
      next_row_sig_1d        <=          'b0            ;
  end
  else
    next_row_sig_1d           <=         'b0            ;
end 

always @(posedge clk or negedge rst_n) begin
  if (!rst_n )    begin
    col_index_1d          <=          'b0             ;
    col_index_2d          <=          'b0             ;
    col_index_3d          <=          'b0             ;
    col_index_4d          <=          'b0             ;
  end
  else begin
    col_index_1d          <=          col_index_fifo  ;
    col_index_2d          <=          col_index_1d    ;
    col_index_3d          <=          col_index_2d    ;
    col_index_4d          <=          col_index_3d    ;
  end
end
//-----------------------------------------------------------------------
//          
//------------------------------------------------------------------------

//--------row count

always @(posedge clk or negedge rst_n) begin
    if(!rst_n || ADDR_col ==0)  //due to the delay of rd
        row_count       <=          15'd0                 ;
    else if (rd_en==1'b1 & EOR_fifo==1'd0)
        row_count       <=          row_count+15'd1       ;
end

//-------diag control

always @(posedge clk or negedge rst_n) begin
  if(!rst_n)
    diag_signal       <=      1'b0              ;
  else if (row_count == col_index_fifo && cnt_en >=2)              
    diag_signal       <=      1'b1              ;
  else
    diag_signal       <=      1'b0              ;
end 

always @(posedge clk or negedge rst_n) begin
  if(!rst_n)  begin
    diag_signal_1d       <=      1'b0              ;
    diag_signal_2d       <=      1'b0              ;
  end
  else begin       
    diag_signal_1d       <=      diag_signal              ;
    diag_signal_2d       <=      diag_signal_1d             ;
  end
end 
//-------upper mult

always @(posedge clk or negedge rst_n) begin
    if(!rst_n )
        partial_upper       <=          'd0                             ;
    else if (rd_en)
        partial_upper       <=          value_fifo   *    x_upper       ;
    else
        partial_upper       <=          'd0                             ;

end    

//-------lower mult

always @(posedge clk or negedge rst_n) begin
    if(!rst_n || diag_signal == 1'b1)
        partial_lower       <=          'd0                             ;
    else if (rd_en)
        partial_lower       <=          value_fifo   *    x_lower       ;
    else
        partial_lower       <=          'd0                             ;
end
//-------summation

always @(posedge clk or negedge rst_n) begin
  if (!rst_n)
    partial_sum_up          <=       'd0        ;
  else if (rd_en) begin
    if (next_row_sig)
      partial_sum_up          <=        'd0             + partial_upper     ;
    else if (diag_signal_2d)
      partial_sum_up          <=        partial_sum_up  + partial_upper + partial_sum_low ;
    else 
      partial_sum_up          <=        partial_sum_up  + partial_upper     ;
  end
end


always @(posedge clk or negedge rst_n) begin
  if (!rst_n)
    partial_sum_low          <=       'd0        ;
  else if (rd_en) begin
    if (same_row_sig)
      partial_sum_low         <=       partial_sum_low + partial_lower ;    // conflict:    consequent read
    else if (col_index_2d == col_index_4d)
      partial_sum_low         <=       write_through_lower + partial_lower ; // conflict:    concurrent read and write
    else 
      partial_sum_low         <=       current_lower + partial_lower ;
  end
end


//------col address


always @(posedge clk or negedge rst_n) begin
    if (!rst_n ||  rd_en==1'd0) 
        ADDR_col          <=       'd0             ;
    else 
        ADDR_col          <=       ADDR_col+'d1 ;
end

//-------val addr
reg [3:0]       ADDR_val                          ;

always @(posedge clk or negedge rst_n) begin
  if (!rst_n)
    ADDR_val          <=        'd0               ;
  else
    ADDR_val          <=        ADDR_col          ;
end


//------final final output-----from up ----
always @(posedge clk or negedge rst_n) begin
  if (!rst_n)
    final_all_output      <=  0;
  else  begin 
  if (next_row_sig)
    final_all_output      <=      partial_sum_up  ;
  if (next_row_sig && next_row_sig_1d)
    final_all_output      <=      partial_sum_up  + partial_sum_low;
  end
end
wire single_element_in_one_row      ;
assign single_element_in_one_row = next_row_sig && next_row_sig_1d  ;

blk_mem_gen_val val (
  .clka(clk),    // input wire clka
  .ena(wire_rd_en_val),      // input wire ena
  .addra(ADDR_val),  // input wire [3 : 0] addra
  .douta(value_fifo)  // output wire [15 : 0] douta
);
blk_mem_gen_col col_index (
  .clka(clk),    // input wire clka
  .ena(rd_en),      // input wire ena
  .addra(ADDR_col),  // input wire [3 : 0] addra
  .douta(col_index_fifo)  // output wire [15 : 0] douta
);
blk_mem_gen_eor eor (
  .clka(clk),    // input wire clka
  .ena(rd_en),      // input wire ena
  .addra(ADDR_col),  // input wire [3 : 0] addra
  .douta(EOR_fifo)  // output wire [15 : 0] douta
);

blk_mem_gen_1 vector_value_ram_upper (
  .clka(clk),    // input wire clka
  .ena(rd_en),      // input wire ena
  .addra(col_index_fifo),  // input wire [3 : 0] addra
  .douta(x_upper)  // output wire [15 : 0] douta
);

blk_mem_gen_2 ram_lower (
  .clka(clk),    // input wire clka
  .ena(rd_en),      // input wire ena
  .addra(wire_row_count),  // input wire [3 : 0] addra
  .douta(x_lower)  // output wire [15 : 0] douta
);

blk_output output_ram (
  .clka(clk),    // input wire clka
  .ena(rd_en),      // input wire ena
  .wea(),      // input wire [0 : 0] wea
  .addra(ADDR_low_sum_read),  // input wire [3 : 0] addra
  .dina(),    // input wire [15 : 0] dina
  .douta(current_lower),  // output wire [15 : 0] douta
  .clkb(clk),    // input wire clkb
  .enb(rd_en),      // input wire enb
  .web(rd_en),      // input wire [0 : 0] web
  .addrb(ADDR_low_sum_write),  // input wire [3 : 0] addrb
  .dinb(partial_sum_low),    // input wire [15 : 0] dinb
  .doutb(write_through_lower)  // output wire [15 : 0] doutb
);

  endmodule
