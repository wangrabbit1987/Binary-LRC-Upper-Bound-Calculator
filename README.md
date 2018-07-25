# Binary-LRC-Upper-Bound-Calculator
###############################
        Introduction
###############################
For given codelength n, distance d, and locality parameter r, this calculator output the upper bound on the dimesnion k for an [n,k,d] binary linear LRC with locality r, based on the results of the follwing paper:
A. Wang, Z. Zhang and D. Lin. "Bounds for Binary Linear Locally Repairable Codes via a Sphere-Packing Approach".
The full version of the paper is available at https://arxiv.org/abs/1701.05989 
###############################
        Installation
###############################
Download 'main.cpp' and 'makefile', then run 'make' in the same directory.
###############################
           Usages
###############################
1.Computing the upper bound on the dimesnion k for an [n,k,d] binary linear LRC with locality r:
  ./upperbound n d r
2.Computing the upper bound in the fast mode. This mode can accelerate the computation for large code parameters, say n>255. However, the output may be not accurate when the code length is too large, say n>400.
  ./upperbound -f n d r
###############################
         Contact
###############################
If you have any questions, please contact Anyu Wang at wanganyu@iie.ac.cn
