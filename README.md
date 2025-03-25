# TraceIt
A simple console 3d physics simulator that I made with C++, should be compatible with any unix-like system (only tested on Linux, Android) and Windows. Heavilly inspired by OneLoneCoder's ConsoleGameEngine

# Dependencies
## Runtime:
Install ncurses (if not yet installed)
On Arch Linux run:
```bash
sudo pacman -Syu ncurses
```
On Debian (based) run:
```bash
sudo apt install ncurses
```
## Build:
Install clang and cmake (if not yet installed)
On Arch Linux run:
```bash
sudo pacman -Syu cmake
```
On Debian (based) run:
```bash
sudo apt install cmake
```
# How to build
Clone this repo and create "build" directory:
```bash
git clone https://github.com/zxcfghjkl/TraceIt.git && cd TraceIt && mkdir build
```
Initialize build and compile TraceIt:
```bash
cmake .. && make
```
Or for Windows:
```bash
cmake .. && cmake --build
```
Then, for the binary to run properly you'll need to put the .obj file in src/ (and PDCurses.dll on Windows) to build/
