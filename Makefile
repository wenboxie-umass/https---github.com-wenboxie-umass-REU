all:
	g++ -std=c++11 HLM_KMP.cpp -o exe_file
	rm HL_KMP.txt
	clear
	./exe_file 11 11
compile:
	g++ -std=c++11 HLM_KMP.cpp -o exe_file
run:
	./exe_file 11 11
show:
	cat HL_KMP.txt
para:
	g++ -Wall -Itrng-4.19 -Ltrng-4.19/src/.libs HLM_KMP_1D_Para_Yao.cpp -o exe_parallel -ltrng4 -fopenmp -std=c++11
	./exe_parallel