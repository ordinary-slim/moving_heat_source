pushd ../../Release
cmake .. -DCMAKE_BUILD_TYPE=Release && make
popd
rm -r post* *.pvd
python3 quadRun.py
