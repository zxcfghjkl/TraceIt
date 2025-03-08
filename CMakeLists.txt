cmake_minimum_required(VERSION 3.10)
project(TraceIt)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OS-specific configuration
if(WIN32)
    # Windows
    add_executable(TraceIt
        src/main.cpp
    )

    # PDCurses для Windows
    target_link_libraries(TraceIt pdcurses)

    # Эквивалент usleep() для Windows
    target_link_libraries(TraceIt winmm)  # Для Sleep()

else()
    # Linux/Mac/Unix
    add_executable(TraceIt
        src/main.cpp
    )

    # Поиск и линкинг ncurses
    find_package(Curses REQUIRED)
    target_include_directories(TraceIt PRIVATE ${CURSES_INCLUDE_DIR})
    target_link_libraries(TraceIt ${CURSES_LIBRARIES})

    # Линкинг cmath
    target_link_libraries(TraceIt m)
endif()

# fstream, iostream, vector
target_link_libraries(TraceIt stdc++fs)
