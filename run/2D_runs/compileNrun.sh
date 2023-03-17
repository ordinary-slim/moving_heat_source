set -e
pushd ../../Release
cmake .. -DCMAKE_BUILD_TYPE=Release && make
popd
time python3 quadRun.py
