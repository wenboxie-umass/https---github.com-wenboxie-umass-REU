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
