# TSPN

## Method Description
TSPN is an effective method for predicting potential SM-miRNA associations using truncated Schatten-p norm. It uses Gaussian kernels to preprocess SM similarity and miRNA similarity and uses truncated schatten-p norm instead of rank functions to significantly improve the prediction accuracy. <a href="https://www.example.com" title="论文链接(待添加)">论文链接</a>

## Operating Environmention
PyCharm == 2021.2.1\
Python == 3.11.2\
Windows == 10\
Processor == Intel(R) Core(TM) i5 8300H 2.30GHz CPU with 16G RAM

## Required Packages
numpy == 1.24.1\
matplotlib == 3.6.2\
scipy == 1.10.0\
scikit-learn == 1.2.0

## File Description

### Dataset:
This folder contains three datasets (dataset 1, dataset 2, and newTest). Taking dataset 1 as an example, it includes the following files:\
**SM number.xlsx:** It contains the CID of 831 SMs and the row number of each SM in the association matrix. For example, the first row of the association matrix represents SM 'CID:137'.\
**miRNA number.xlsx:** It contains the name of 541 human-related miRNAs and the column number of each miRNA in the association matrix. For example, the first column of the association matrix represents miRNA 'hsa-let-7a-1'.\
**miRNA similarity matrix.txt:** It contains the integrated miRNA similarity (the dimension is 541 $\times$ 541), where each row and its corresponding column represent a specific miRNA. The element in the i-th row and j-th column denotes the similarity score between 
$miRNA_{i}$ and $miRNA_{j}$.\
**SM-miRNA-confirmed associations.txt:** It contains 664 SM-miRNA known associations. For example, the first row in the file, (75,177), denotes that the SMs represented in the 76-th row of the association matrix, 'CID:2578', is associated with the miRNA represented in the 178-th column of the association matrix, 'hsa-mir-200c'.\
**SM-miRNA association matrix.txt:** It is the constructed association matrix (the dimension is 831 $\times$ 541), each row of which represents a specific SM, and each column represents a specific miRNA. The (i,j)-th element of the association matrix, 
$m_{ij}$, is set to 1 if $SM_{i}$ is associated with $miRNA_{j}$, otherwise it is set to 0. The matrix has a total of 664 ones and 448907 zeros.\
**Adjacency matrix M.txt:** It is the adjacency matrix of the constructed SM-miRNA heterogeneous network (the dimension is 1342 
 $\times$ 1342), which is the input of TSPN.py.
### Calculate GIPK similarity.py:
This file contains the Python code to calculate the GIPK similarity of the SM-miRNA association matrix. The input is the SM-miRNA association matrix, the SM similarity matrix and the miRNA similarity matrix. The output is the GIPK similarity matrix of SM and miRNA.
### TSPN.PY:
This file contains the Python code for our algorithm. The inputs are the adjacency matrix (the dimension of the adjacency matrix in Dataset 1 is 1342 $\times$ 1342) and the SM similarity matrix (the dimension of the SM similarity matrix in Dataset 1 is 831 $\times$ 831. The output is the prediction score matrix (831 $\times$ 541).
### Running steps:
Step 1: we run “Calculate GIPK similarity.py”. Its input is the SM-miRNA association matrix, the SM similarity matrix and the miRNA similarity matrix. Its output is is the GIPK similarity matrix of SM and miRNA.\
Step 2: we execute “TSPN.py”. Its input is the adjacency matrix of the heterogeneous SM-miRNA network and the integrated SM similarity matrix (used for matrix division). Its output is the SM-miRNA predictive score matrix.
### Contact
If you have any problems or find mistakes, please feel free to contact me: z22070050@s.upc.edu.cn
