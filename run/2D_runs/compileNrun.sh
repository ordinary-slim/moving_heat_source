pushd ../../Debug
cmake .. -DCMAKE_BUILD_TYPE=Debug && make
popd
pushd ../../Release
cmake .. -DCMAKE_BUILD_TYPE=Release && make
popd
rm -r post* *.pvd
python3 quadRun.py
