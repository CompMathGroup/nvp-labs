#!/bin/bash

cd ..

rm -rf build-{l,w}in{32,64}

mkdir -p build-{l,w}in{32,64}

(cd build-lin32 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cmake/linux-i686.cmake .. && make -j8 && cp ellipt ../release/ellipt32)
(cd build-lin64 && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j8 && cp ellipt ../release/ellipt64) &
(cd build-win32 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cmake/windows-i686.cmake .. && make -j8 && cp ellipt.exe ../release/ellipt32.exe)
(cd build-win64 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cmake/windows-x86_64.cmake .. && make -j8 && cp ellipt.exe ../release/ellipt64.exe)

cp *.geom release

cd release

zip Lab3.zip ellipt32* ellipt64* *.geom
