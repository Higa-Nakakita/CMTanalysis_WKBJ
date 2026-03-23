CMP = gfortran
#OPT1 = -O3
OPT2 = -no-pie  #to eliminate error 'R_X86_64_32 against `.rodata' can not be used when making a PIE object' #20251024 This option in invalid in minarnda (unrecognized command line option ‘-no-pie’)
OPT3 = -mcmodel=large #to eliminate error "再配置がオーバーフローしないように切り詰められました"
sac_lib = /opt/sac/lib/sacio.a /opt/sac/lib/libsacio.a /opt/sac/lib/libsac.a #for rsac1, wsac1
EXECUTE_file = SWFI_CMT_WKBJ_260109

BIN = .

all:
	$(CMP) $(OPT2) $(OPT3) -o $(BIN)/$(EXECUTE_file) *.f $(sac_lib)
clean:
	rm -f $(EXECUTE_file)
