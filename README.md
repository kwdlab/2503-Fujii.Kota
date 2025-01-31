# Overview
This project uses the quasi-isomorphic cryptographic library Microsoft SEAL to accelerate the computation. The program will compare the change in computation time and maximum memory usage using the matrix product algorithm and the Strassen algorithm, respectively.

# Description
- The bit_matrix.cpp is in the seal_project file and can be calculated using the matrix product algorithm.
- The bit_strassen.cpp is in the seal_project file and can be calculated using Strassen's algorithm.
- The two cpp programs can change the matrix size, quasi-isomorphic cryptographic parameters, and polynomial degree values in addition to the numerical input.
- When the program is executed, the computation time and maximum memory usage are displayed.

# Requirements
- macOS version 14.4.1
- Microsoft SEAL 4.1
- Apple Clang 15.0.0(clang-1500.3.9.4）

# Install/Usage
- Install of Microsoft SEAL
'git clone https://github.com/microsoft/SEAL.git
cd SEAL'

- Install of Homebrew
'brew install gnu-time'

- Usage
Step1.Compile
'g++ -std=c++17 bit_strassen.cpp -o seal_example -I/usr/local/include/SEAL-4.1 -L/usr/local/lib -lseal-4.1 -mmacosx-version-min=14.4'

Step2.Run
'/opt/homebrew/bin/gtime --verbose ./seal_example'

# Author
[Kouta Fujii](https://github.com/kwdlab/2503-Fujii.Kota)

# References
- [Microsoft SEAL](https://www.microsoft.com/en-us/research/project/microsoft-seal/)
- [Strassen](https://ja.wikipedia.org/wiki/シュトラッセンのアルゴリズム)

# License
[MIT](https://opensource.org/license/mit) 
