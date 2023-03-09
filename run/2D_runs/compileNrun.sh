set -e
pushd ../../Release
cmake .. -DCMAKE_BUILD_TYPE=Release && make
popd
python3 quadRun.py
