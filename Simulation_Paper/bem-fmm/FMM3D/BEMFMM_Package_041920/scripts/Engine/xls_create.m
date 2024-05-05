function CE_int = xls_create(MNE_out)

CE_int = readcell('16chs_rack_input.xls');
row_dis = ones(1,32);

dI_dt_fin = MNE_out.dI_dt/1e6;
int_MSO_tmp = 100*dI_dt_fin/140;
int_MSO = upsample(int_MSO_tmp,2);

CE_int(3,6:37) = mat2cell(floor(int_MSO),row_dis);
CE_int(2,1:36) = mat2cell(zeros(36,1),ones(36,1));

writecell(CE_int,'16chs_rack_output.xls') 
end