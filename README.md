# TraceIt
A simple console 3d physics simulator that I made with C++, should be compatible with any unix-like system (only tested on Linux, Android), Windows.

# Dependencies
## Runtime:
Install ncurses (if not yet installed)
On Arch Linux run:
```bash
$ sudo pacman -Syu ncurses
```
## Build:
Install clang and cmake (if not yet installed)
On Arch Linux run:
```bash
$ sudo pacman -Syu clang cmake
```
# How to build
Clone this repo, create build directory in it and run (in the build/):
```bash
$ cmake .. && make
```
Or for Windows:
```bash
$ cmake .. && cmake --build
```
Then, for the binary to run properly you'll need to put the .obj file in src/ to build/ (not yet implemented in the CMakeLists.txt)
