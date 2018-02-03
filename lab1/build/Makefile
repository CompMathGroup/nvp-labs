.PHONY: unix32 unix64 win32 win64 osx

# all: unix32 unix64 win32 win64 osx
all: unix32 unix64 win32 win64

OSXCROSS=/home/uranix/osxcross/target/bin/

unix32:
	env PREFIX=`pwd`/$@ CXX="g++ -m32" make -C ../core all
	env PREFIX=`pwd`/$@ CXX="g++ -m32" LIBS="-Wl,-Bdynamic -lML32i3 -lrt" ARCH=unix make -C ../mlink all

unix64:
	env PREFIX=`pwd`/$@ CXX="g++ -m64" make -C ../core all
	env PREFIX=`pwd`/$@ CXX="g++ -m64" LIBS="-Wl,-Bdynamic -lML64i3 -lrt" ARCH=unix make -C ../mlink all

osx:
	env PREFIX=`pwd`/$@ CXX="${OSXCROSS}x86_64-apple-darwin15-clang++ -stdlib=libstdc++" AR="${OSXCROSS}x86_64-apple-darwin15-ar" make -C ../core all
	env PREFIX=`pwd`/$@ CXX="${OSXCROSS}x86_64-apple-darwin15-clang++ -stdlib=libstdc++" LIBS="-lMLi3 -framework Foundation" ARCH=osx make -C ../mlink all

win32:
	env PREFIX=`pwd`/$@ CC="i686-w64-mingw32-gcc -mwindows" CXX="i686-w64-mingw32-g++ -mwindows" make -C ../core all
	env PREFIX=`pwd`/$@ CC="i686-w64-mingw32-gcc -mwindows" CXX="i686-w64-mingw32-g++ -mwindows" LIBS="-static-libstdc++ -static-libgcc -lml32i3m" ARCH=win make -C ../mlink all

win64:
	env PREFIX=`pwd`/$@ CC="x86_64-w64-mingw32-gcc -mwindows" CXX="x86_64-w64-mingw32-g++ -mwindows" make -C ../core all
	env PREFIX=`pwd`/$@ CC="x86_64-w64-mingw32-gcc -mwindows" CXX="x86_64-w64-mingw32-g++ -mwindows" LIBS="-static-libstdc++ -static-libgcc -lml64i3m" ARCH=win make -C ../mlink all

clean_unix32:
	env PREFIX=`pwd`/unix32 make -C ../core clean
	env PREFIX=`pwd`/unix32 make -C ../mlink clean

clean_unix64:
	env PREFIX=`pwd`/unix64 make -C ../core clean
	env PREFIX=`pwd`/unix64 make -C ../mlink clean

clean_osx:
	env PREFIX=`pwd`/osx make -C ../core clean
	env PREFIX=`pwd`/osx make -C ../mlink clean

clean_win32:
	env PREFIX=`pwd`/win32 make -C ../core clean
	env PREFIX=`pwd`/win32 make -C ../mlink clean

clean_win64:
	env PREFIX=`pwd`/win64 make -C ../core clean
	env PREFIX=`pwd`/win64 make -C ../mlink clean

clean: clean_unix32 clean_unix64 clean_win32 clean_win64 clean_osx
