// based on shrtdct.c
// https://www.kurims.kyoto-u.ac.jp/~ooura/fft-j.html


// w4r = sin(pi/4)
// ckr = sqrt(2/8) * cos(pi/2*k/8)
// cki = sqrt(2/8) * sin(pi/2*k/8)

`define w4r $signed(32'b 00000000000000000000000010110101)
`define c1r $signed(32'b 00000000000000000000000001111101)
`define c1i $signed(32'b 00000000000000000000000000011000)
`define c2r $signed(32'b 00000000000000000000000001110110)
`define c2i $signed(32'b 00000000000000000000000000110000)
`define c3r $signed(32'b 00000000000000000000000001101010)
`define c3i $signed(32'b 00000000000000000000000001000111)
`define c4r $signed(32'b 00000000000000000000000001011010)



module transform_quantize(
        input clk,
        input nrst,
        input signed [31:0] din,
        output signed [31:0] dout,
        output dout_en,
        output accept_in
    );
	 
	 
function signed [31:0] abpcd ( // ab+cd
        input signed [31:0] a,
        input signed [31:0] b,
        input signed [31:0] c,
        input signed [31:0] d
    );

    abpcd = (a*b + c*d) >>> 8;
endfunction


function signed [31:0] abmcd ( // ab-cd
        input signed [31:0] a,
        input signed [31:0] b,
        input signed [31:0] c,
        input signed [31:0] d
    );
    abmcd = (a*b - c*d) >>> 8;
endfunction

function signed [31:0] abpc( // a*(b+c)
        input signed [31:0] a,
        input signed [31:0] b,
        input signed [31:0] c
    );
    abpc = ( a * (b + c) ) >>> 8;
endfunction

function signed [31:0] abmc( // a*(b-c)
        input signed [31:0] a,
        input signed [31:0] b,
        input signed [31:0] c
    );
    abmc = ( a * (b - c) ) >>> 8;
endfunction

function signed [31:0] mulq( //  multiply by b and round, b = reciprocal of quantization coeff.
        input signed [31:0] a,
        input signed [31:0] b
    );
    mulq = ( a * b + 32768 ) >>> 16;
endfunction



function signed [31:0] mulq_inv( multiply by b, b = quantization coeff.
        input signed [31:0] a,
        input signed [31:0] b
    );
    mulq_inv = (a * b);
endfunction



    reg [1:0] state; // 0: input, 1: dct->quantize->idct, 2: output
    reg [7:0] stage;
    reg signed [31:0] a[0:63]; // input data
    reg [5:0] idx_in, idx_out;

    assign dout_en = (state == 2'd2);
    assign dout = (state == 2'd2 ? a[idx_out] : 32'bz);
    assign accept_in = (state == 2'd0);

    reg signed [31:0] qtab [0:63]; // quantization table
    reg signed [31:0] qtab_inv [0:63]; // reciprocal of quantization table


    // initialize quantization table
    // https://www.jstage.jst.go.jp/article/itej/67/2/67_136/_pdf
    initial begin
        $readmemb("/home/sysele2303/eeic-exp2-fpga/quantize_table.bin", qtab);
        $readmemb("/home/sysele2303/eeic-exp2-fpga/quantize_table_inv.bin", qtab_inv);
    end

    // auxiliary registers for dct, idct
    reg signed [31:0] x0r [0:7];
    reg signed [31:0] x0i [0:7];
    reg signed [31:0] x1r [0:7];
    reg signed [31:0] x1i [0:7];
    reg signed [31:0] x2r [0:7];
    reg signed [31:0] x2i [0:7];
    reg signed [31:0] x3r [0:7];
    reg signed [31:0] x3i [0:7];


    integer j;
    always @(posedge clk) begin
        if (!nrst) begin
            idx_in <= 0;
            idx_out <= 0;
            state <= 0;
            stage <= 0;
        end else begin

            // input
            if (state == 0) begin
                a[idx_in] <= din;
                idx_in <= idx_in+6'd1;
                if (idx_in == 63) begin
                    state <= 1;
                end
            end

            else if (state == 1) begin

                // dct: x-axis
                if (stage == 0) begin
                    stage <= 1;
                    for (j = 0; j <= 7; j = j+1) begin
                        x1r[j] <= abpcd(`c1r, a[ 8+j], `c1i, a[56+j]);
                        x1i[j] <= abmcd(`c1r, a[56+j], `c1i, a[ 8+j]);
                        x3r[j] <= abpcd(`c3r, a[24+j], `c3i, a[40+j]);
                        x3i[j] <= abmcd(`c3r, a[40+j], `c3i, a[24+j]);
                        x2r[j] <= abpcd(`c2r, a[16+j], `c2i, a[48+j]);
                        x2i[j] <= abmcd(`c2r, a[48+j], `c2i, a[16+j]);
                        x0r[j] <= abpc(`c4r, a[j], a[32+j]);
                        x0i[j] <= abmc(`c4r, a[j], a[32+j]);
                    end
                end
                else if(stage == 1) begin
                    stage <= 2;
                    for (j=0; j<=7; j=j+1) begin
                        x2r[j] <= x0r[j] - x2r[j];
                        x2i[j] <= x0i[j] - x2i[j];
                        x0r[j] <= x0r[j] + x2r[j];
                        x0i[j] <= x0i[j] + x2i[j];
                        x3r[j] <= x1r[j] - x3r[j];
                        x1i[j] <= x1i[j] + x3i[j];
                        x1r[j] <= x1r[j] + x3r[j];
                        x3i[j] <= x3i[j] - x1i[j];
                    end
                end
                else if (stage == 2) begin
                    stage <= 3;
                    for (j = 0; j <= 7; j = j+1) begin
                        x1i[j] <= abpc(`w4r, x3r[j], x1i[j]);
                        x3r[j] <= abmc(`w4r, x3r[j], x1i[j]);
                    end
                end
                else if (stage == 3) begin
                    stage <= 4;
                    for (j = 0; j <= 7; j = j+1) begin
                        a[   j] <= x0r[j] + x1r[j];
                        a[56+j] <= x0r[j] - x1r[j];
                        a[16+j] <= x0i[j] + x1i[j];
                        a[40+j] <= x0i[j] - x1i[j];
                        a[32+j] <= x2r[j] - x3i[j];
                        a[24+j] <= x2r[j] + x3i[j];
                        a[48+j] <= x2i[j] - x3r[j];
                        a[ 8+j] <= x2i[j] + x3r[j];
                    end
                end

                // dct: y-axis
                else if (stage == 4) begin
                    stage <= 5;
                    for (j = 0; j <= 7; j = j+1) begin
                        x1r[j] <= abpcd(`c1r, a[1+(j<<3)], `c1i, a[7+(j<<3)]);
                        x1i[j] <= abmcd(`c1r, a[7+(j<<3)], `c1i, a[1+(j<<3)]);
                        x3r[j] <= abpcd(`c3r, a[3+(j<<3)], `c3i, a[5+(j<<3)]);
                        x3i[j] <= abmcd(`c3r, a[5+(j<<3)], `c3i, a[3+(j<<3)]);
                        x2r[j] <= abpcd(`c2r, a[2+(j<<3)], `c2i, a[6+(j<<3)]);
                        x2i[j] <= abmcd(`c2r, a[6+(j<<3)], `c2i, a[2+(j<<3)]);
                        x0r[j] <= abpc(`c4r, a[j<<3], a[4+(j<<3)]);
                        x0i[j] <= abmc(`c4r, a[j<<3], a[4+(j<<3)]);
                    end
                end
                else if(stage == 5) begin
                    stage <= 6;
                    for (j=0; j<=7; j=j+1) begin
                        x2r[j] <= x0r[j] - x2r[j];
                        x2i[j] <= x0i[j] - x2i[j];
                        x0r[j] <= x0r[j] + x2r[j];
                        x0i[j] <= x0i[j] + x2i[j];
                        x3r[j] <= x1r[j] - x3r[j];
                        x1i[j] <= x1i[j] + x3i[j];
                        x1r[j] <= x1r[j] + x3r[j];
                        x3i[j] <= x3i[j] - x1i[j];
                    end
                end
                else if (stage == 6) begin
                    stage <= 7;
                    for (j = 0; j <= 7; j = j+1) begin
                        x1i[j] <= abpc(`w4r, x3r[j], x1i[j]);
                        x3r[j] <= abmc(`w4r, x3r[j], x1i[j]);
                    end
                end
                else if (stage == 7) begin
                    stage <= 8;
                    for (j = 0; j <= 7; j = j+1) begin
                        a[  (j<<3)] <= x0r[j] + x1r[j];
                        a[7+(j<<3)] <= x0r[j] - x1r[j];
                        a[2+(j<<3)] <= x0i[j] + x1i[j];
                        a[5+(j<<3)] <= x0i[j] - x1i[j];
                        a[4+(j<<3)] <= x2r[j] - x3i[j];
                        a[3+(j<<3)] <= x2r[j] + x3i[j];
                        a[6+(j<<3)] <= x2i[j] - x3r[j];
                        a[1+(j<<3)] <= x2i[j] + x3r[j];
                    end
                end

                // quantize
                 else if (stage == 8) begin
                    stage <= 9;
                    for (j=0; j<=63; j=j+1) begin
                        a[j] <= mulq(a[j], qtab_inv[j]);
                    end
                 end
                  else if (stage == 9) begin
                    stage <= 10;
                    for (j=0; j<=63; j=j+1) begin
                        a[j] <= mulq_inv(a[j], qtab[j]);
                    end
                 end

                // idct: x-axis
                else if (stage == 10) begin
                    stage <= 11;
                    for (j = 0; j <= 7; j = j+1) begin
                        x0r[j] <= a[   j] + a[56+j];
                        x1r[j] <= a[   j] - a[56+j];
                        x0i[j] <= a[16+j] + a[40+j];
                        x1i[j] <= a[16+j] - a[40+j];
                        x2r[j] <= a[32+j] + a[24+j];
                        x3r[j] <= a[32+j] - a[24+j];
                        x2i[j] <= a[48+j] + a[ 8+j];
                        x3i[j] <= a[48+j] - a[ 8+j];
                    end
                end
                else if(stage == 11) begin
                    stage <= 12;
                    for (j=0; j<=7; j=j+1) begin
                        x0r[j] <= x0r[j] + x2r[j];
                        x2i[j] <= x0i[j] + x2i[j];
                        x2r[j] <= x0r[j] - x2r[j];
                        x0i[j] <= x0i[j] - x2i[j];
                        x3i[j] <= abmc(`w4r, x1i[j], x3i[j]);
                        x1i[j] <= abpc(`w4r, x1i[j], x3i[j]);
                    end
                end
                else if (stage == 12) begin
                    stage <= 13;
                    for (j = 0; j <= 7; j = j+1) begin
                        x3i[j] <= x1i[j] - x3r[j];
                        x1i[j] <= x1i[j] + x3r[j];
                        x3r[j] <= x1r[j] - x3i[j];
                        x1r[j] <= x1r[j] + x3i[j];
                    end
                end
                else if (stage == 13) begin
                    stage <= 14;
                    for (j = 0; j <= 7; j = j+1) begin
                        a[   j] <= abpc(`c4r, x0r[j], x2i[j]);
                        a[32+j] <= abmc(`c4r, x0r[j], x2i[j]);
                        a[16+j] <= abmcd(`c2r, x2r[j], `c2i, x0i[j]);
                        a[48+j] <= abpcd(`c2r, x0i[j], `c2i, x2r[j]);
                        a[ 8+j] <= abmcd(`c1r, x1r[j], `c1i, x1i[j]);
                        a[56+j] <= abpcd(`c1r, x1i[j], `c1i, x1r[j]);
                        a[24+j] <= abmcd(`c3r, x3r[j], `c3i, x3i[j]);
                        a[40+j] <= abpcd(`c3r, x3i[j], `c3i, x3r[j]);
                    end
                end

                // idct: y-axis
                else if (stage == 14) begin
                    stage <= 15;
                    for (j = 0; j <= 7; j = j+1) begin
                        x0r[j] <= a[  (j<<3)] + a[7+(j<<3)];
                        x1r[j] <= a[  (j<<3)] - a[7+(j<<3)];
                        x0i[j] <= a[2+(j<<3)] + a[5+(j<<3)];
                        x1i[j] <= a[2+(j<<3)] - a[5+(j<<3)];
                        x2r[j] <= a[4+(j<<3)] + a[3+(j<<3)];
                        x3r[j] <= a[4+(j<<3)] - a[3+(j<<3)];
                        x2i[j] <= a[6+(j<<3)] + a[1+(j<<3)];
                        x3i[j] <= a[6+(j<<3)] - a[1+(j<<3)];
                    end
                end
                else if(stage == 15) begin
                    stage <= 16;
                    for (j=0; j<=7; j=j+1) begin
                        x0r[j] <= x0r[j] + x2r[j];
                        x2i[j] <= x0i[j] + x2i[j];
                        x2r[j] <= x0r[j] - x2r[j];
                        x0i[j] <= x0i[j] - x2i[j];
                        x3i[j] <= abmc(`w4r, x1i[j], x3i[j]);
                        x1i[j] <= abpc(`w4r, x1i[j], x3i[j]);
                    end
                end
                else if (stage == 16) begin
                    stage <= 17;
                    for (j = 0; j <= 7; j = j+1) begin
                        x3i[j] <= x1i[j] - x3r[j];
                        x1i[j] <= x1i[j] + x3r[j];
                        x3r[j] <= x1r[j] - x3i[j];
                        x1r[j] <= x1r[j] + x3i[j];
                    end
                end
                else if (stage == 17) begin
                    stage <= 0;
                    state <= 2;
                    for (j = 0; j <= 7; j = j+1) begin
                        a[  (j<<3)] <= abpc(`c4r, x0r[j], x2i[j]);
                        a[4+(j<<3)] <= abmc(`c4r, x0r[j], x2i[j]);
                        a[2+(j<<3)] <= abmcd(`c2r, x2r[j], `c2i, x0i[j]);
                        a[6+(j<<3)] <= abpcd(`c2r, x0i[j], `c2i, x2r[j]);
                        a[1+(j<<3)] <= abmcd(`c1r, x1r[j], `c1i, x1i[j]);
                        a[7+(j<<3)] <= abpcd(`c1r, x1i[j], `c1i, x1r[j]);
                        a[3+(j<<3)] <= abmcd(`c3r, x3r[j], `c3i, x3i[j]);
                        a[5+(j<<3)] <= abpcd(`c3r, x3i[j], `c3i, x3r[j]);
                    end
                end
            end

            // output
            else if (state == 2) begin
                idx_out <= idx_out+6'd1;
                if (idx_out == 63) begin
                    state <= 0;
                end
            end
        end
    end

endmodule // transform_quantize

