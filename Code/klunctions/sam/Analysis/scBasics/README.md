## Mac OSX issues
### Python package installation error related to gcc
Try the following (this example is for installing the python package `annoy`): 
1. Install homebrew (https://brew.sh/)
2. Install gcc: `brew install gcc`
3. `env CXX=/usr/local/Cellar/gcc/8.2.0/bin/g++-8 CC=/usr/local/Cellar/gcc/8.2.0/bin/gcc-8 pip install annoy`

