# obme 
A fork of the Intel-NTU OpenDVC codec developed by Prof. Shao-Yi Chien and Dr. Chia-Han Lee. Official website: http://media.ee.ntu.edu.tw/research/opendvc/

######################################################################################################################################################

Install Instructions (release):
mkdir b
cd b
cmake ..
make

Install Instructions (debug):
mkdir b
cd b
cmake .. -DCMAKE_BUILD_TYPE=DEBUG
make

Install generates bin files, unzips ldpca ladder files, and sets local paths in config.h file.

######################################################################################################################################################

Changes:
1) Changed build tools (make to cmake).

2) Modified SI generation steps: Motion vectors now predicted using Chroma-plane motion estimation and applied to Luma using MC. This enables
   unidirection decoding hierarchy.

3) Experimenting with additional modules:
    a) OBMC -  (h.263 Advance Prediction mode: a form of overlapped block motion estimation to remove blocking artifacts),
    b) Advanced Precision - higher precision motion vectors,
    c) Block matching algorithm (testing exhaustive search vs. modified spiral search).
