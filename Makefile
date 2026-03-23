CMP = gfortran
EXE = mk_divdel.out
all:
	${CMP} -o ${EXE} mk_divdel_info.f90 subsphere_tmp.f subcalc.f reflect_inc.f
clean:
	rm -f ${EXE}
        	
