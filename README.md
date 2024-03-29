# **Support Vector Laboratory** <img src="man/figures/Logo.png" align="right" width="150" />
             
![GitHub RCMDCHECK](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)
![GitHub ISSUES](https://img.shields.io/github/issues/define957/SupportVectorLab)
![GitHub STARS](https://img.shields.io/github/stars/define957/SupportVectorLab)
![GitHub FORKS](https://img.shields.io/github/forks/define957/SupportVectorLab)
![GitHub LICENSE](https://img.shields.io/github/license/define957/SupportVectorLab)
***


## Introduction

There are many different types of SVMs in this repository. 

## Why I created this project ?

In order to learn `SVMs` better, I built this repository to implement support vector machines. The package is under active development.

## How to install Support Vector Laboratory

Make sure you have installed `devtools`, if you don't have `devtools` installed, please run the following command first 
```{r}
install.packages("devtools")
```
Then you need to install `Rtools` to compile the `C++` code. 

After you installed `devtools` and `Rtools`, please run the following command :
```{r}
devtools::install_github("define957/SupportVectorLab")
```
Then you can have `SupportVectorLab` package on your PC。

## Optimization
+ The dual form is solved by CLipDCD solver (Cliping Dual Coordinate Descent Method, implemented by RcppArmadillo).


## SVMs for classification

+ Hinge Loss Support Vector Classification
+ Pinball Loss Support Vector Classification
+ Least Squares Loss Support Vector Classification
+ Rescaled Quantile Loss Support Vector Classification
+ Bounded Quantile Loss Support Vector Classification
+ Asymmetric Least Squares Loss Support Vector Classification
+ Square Hinge Loss Support Vector Classification

## SVMs for regression
+ Epsilon Insensitive Loss Support Vector Regression
+ Least Squares Loss Support Vector Regression
+ Bounded Quantile Loss Support Vector Regression
+ Asymmetric Least Squares Support Vector Regression

## Twin support vector machines for classification (TSVM)
+ Hinge-TSVM
+ LS-TSVM

## Twin support vector machines for regression (TSVR)
+ Hinge-TSVR

## Multiple birth support vector machine (MBSVM)
+ Hinge-MBSVM
+ Pin-MBSVM
+ Ramp-MBSVM
+ LS-MBSVM


## Kernel options

+ Linear kernel
+ RBF kernel
+ Polynomial kernel

## Development environment and dependency

My enviroment: R 4.3.0, windows 11 x64 &#x2705;

Other test environment detail: 
+ Windows 10/11 x64 &#x2705;
+ Mac osx (ARM platform) &#x2705; 
+ Linux : We haven't tested it yet &#x2753;

Dependency: 

+ Rcpp
+ RcppArmadillo

## Bug report

If you find bug in this package, please post an issue on the [issue](https://github.com/define957/SupportVectorLab/issues) website.


We welcome all suggestions and will improve them in the future!

## Contact us

&#x2709; Email : zhangjiaqi957957@outlook.com

## Licenses

GNU GENERAL PUBLIC LICENSE Version 3 (GPL-3.0)

# Citing Support Vector Laboratory
```
@article{ZHANG2024122759,
title = {Bounded quantile loss for robust support vector machines-based classification and regression},
journal = {Expert Systems with Applications},
volume = {242},
pages = {122759},
year = {2024},
issn = {0957-4174},
doi = {https://doi.org/10.1016/j.eswa.2023.122759},
url = {https://www.sciencedirect.com/science/article/pii/S095741742303261X},
author = {Jiaqi Zhang and Hu Yang},
}
```