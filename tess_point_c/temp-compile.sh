gcc -fPIC -Iinclude src/vec.c src/mat_ra3.c src/tess_stars2px.c -c -lm
gcc -shared -Wl,-soname,libtess_stars2px.so.1 -o libtess_stars2px.so.1.0.1  tess_stars2px.o mat_ra3.o vec.o