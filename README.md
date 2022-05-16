# 22praktikum

## Installing GTest:

git clone https://github.com/google/googletest
cd googletest
mkdir build
cmake .. -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/path/to/your/gtest/installation/dir

## installing Eigen (Ubuntu)
sudo apt install libeigen3-dev

## Setting up:
git clone https://github.com/guidokanschat/22praktikum/
cd 22praktikum
mkdir build
cd build
mkdir YourBuildName
cd YourBuildName
cmake ../.. -DGTEST_ROOT=/path/to/your/gtest/installation/dir
