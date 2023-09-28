cmake_minimum_required(VERSION 3.5...3.26)
project(example)

#find_package(pybind11 REQUIRED)  # or 

# 設定編譯選項
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DNDEBUG -DBAMG_LONG_LONG -DNCHECKPTR -fPIC")
find_package(ZLIB REQUIRED)

# set(PYTHON_EXECUTABLE "/Users/andylee/Library/Caches/pypoetry/virtualenvs/freefem-webplot-PdjqP-PO-py3.10/bin/python")
# set(PYTHON_INCLUDE_DIRS "/Users/andylee/Library/Caches/pypoetry/virtualenvs/freefem-webplot-PdjqP-PO-py3.10")
# find_package(Python3 COMPONENTS Interpreter Development)
# find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
# include_directories(${PYTHON_INCLUDE_DIRS})

add_subdirectory(extern/pybind11)

# 設定頭文件目錄
include_directories(
    /usr/local/lib/ff++/4.13/include
    /Users/andylee/Work/FreeFem-sources/3rdparty/include
)

# 添加源文件
add_library(pyfem MODULE pyfem.cpp)

# 設定輸出目錄和動態連接選項
set_target_properties(pyfem PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  LIBRARY_OUTPUT_NAME pyfem
  PREFIX ""
)

# 在Mac上設定動態連接選項
if(APPLE)
  set_target_properties(pyfem PROPERTIES
    LINK_FLAGS "-bundle -undefined dynamic_lookup"
  )
endif()

target_link_libraries(pyfem PRIVATE pybind11::embed)
target_link_libraries(pyfem PRIVATE ZLIB::ZLIB)
link_libraries(pyfem python3.10)