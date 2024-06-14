clean:
	rm -f lib/param.o lib/param_py.o lib/param.so
	rm -rf firstTest/

lib/param.o: lib/param.F90
	gfortran -c -fPIC -o lib/param.o lib/param.F90

lib/param_py.o: lib/param_py.F90
	gfortran -c -fPIC -o lib/param_py.o lib/param_py.F90

lib/param.so: lib/param.o lib/param_py.o
	gfortran -shared -g -o lib/param.so lib/param_py.o lib/param.o

run: lib/param.so
	./comp_test.py
