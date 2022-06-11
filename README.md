# 22praktikum

## Installing GTest:
```bash
git clone https://github.com/google/googletest 
cd googletest 
mkdir build 
cmake .. -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/path/to/your/gtest/installation/dir 
```
## Installing Eigen
### (Ubuntu)
`sudo apt install libeigen3-dev`  <br>
### CentOS / REHL / Fedora
`sudo dnf install eigen3-devel`

## Setting up:
```bash
git clone https://github.com/guidokanschat/22praktikum/ 
cd 22praktikum 
mkdir build 
cd build 
mkdir YourBuildName
cd YourBuildName
cmake ../.. -DGTEST_ROOT=/path/to/your/gtest/installation/dir
```

## Using doxygen comments:
```c++
/// -> before function / class / module declaration
/// @param a: int value a
/// @param b: int value b
/// @returns : sum of a + b
double foobar(int a, int b){
  ...
  return a+b
}
```
run doxygen in root dir of the project to compile. Output will be at doc/html/index.html
