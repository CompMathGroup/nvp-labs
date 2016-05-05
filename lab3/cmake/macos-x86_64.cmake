set(CMAKE_SYSTEM_NAME Darwin)

set(OSXCROSS "/home/uranix/osxcross/target/bin/")

set(CMAKE_C_COMPILER "${OSXCROSS}x86_64-apple-darwin15-clang")
set(CMAKE_CXX_COMPILER "${OSXCROSS}x86_64-apple-darwin15-clang++")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++" CACHE STRING "c++ flags")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -stdlib=libc++" CACHE STRING "c flags")
