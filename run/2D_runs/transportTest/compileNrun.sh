pushd ../../../Debug
cmake .. -DCMAKE_BUILD_TYPE=Debug && make
popd
rm -r post* *.pvd
python3 main.py
