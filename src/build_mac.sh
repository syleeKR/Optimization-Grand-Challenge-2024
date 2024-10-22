export CC=/opt/homebrew/bin/gcc-14
export CXX=/opt/homebrew/bin/g++-14

rm -r build

python setup.py build_ext --inplace

