.PHONY: unix32 unix64 win32 win64

all: unix32 unix64 win32 win64

unix32:
	env PREFIX=`pwd`/$@ CXX="g++ -m32" make -C ../core all
	env PREFIX=`pwd`/$@ CXX="g++ -m32" LIBS="-lML32i3 -lrt" make -C ../mlink all

unix64:
	env PREFIX=`pwd`/$@ CXX="g++ -m64" make -C ../core all
	env PREFIX=`pwd`/$@ CXX="g++ -m64" LIBS="-lML64i3 -lrt" make -C ../mlink all

win32:
	env PREFIX=`pwd`/$@ CXX="i686-w64-mingw32-g++" make -C ../core all
	env PREFIX=`pwd`/$@ CXX="i686-w64-mingw32-g++" LIBS="-lml32i3m" make -C ../mlink all

win64:
	env PREFIX=`pwd`/$@ CXX="x86_64-w64-mingw32-g++" make -C ../core all
	env PREFIX=`pwd`/$@ CXX="x86_64-w64-mingw32-g++" LIBS="-lml64i3m" make -C ../mlink all

clean_unix32:
	env PREFIX=`pwd`/unix32 make -C ../core clean
	env PREFIX=`pwd`/unix32 make -C ../mlink clean

clean_unix64:
	env PREFIX=`pwd`/unix64 make -C ../core clean
	env PREFIX=`pwd`/unix64 make -C ../mlink clean

clean_win32:
	env PREFIX=`pwd`/win32 make -C ../core clean
	env PREFIX=`pwd`/win32 make -C ../mlink clean

clean_win64:
	env PREFIX=`pwd`/win64 make -C ../core clean
	env PREFIX=`pwd`/win64 make -C ../mlink clean

clean: clean_unix32 clean_unix64 clean_win32 clean_win64
