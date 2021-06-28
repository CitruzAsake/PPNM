Practical Programming and Numerical Methods 2021 - C
In this repository the solutions for the Exercises, Homeworks and Exam-project for the Spring semester of 2021 for Akgash Sundaralingam are presented.

Student information:
Name: Akgash Sundaralingam
Student ID: 201706696
AU ID: au593102
Student mail: 201706696@post.au.dk

My student number ends with 96 I have the exam question 8: One-sided Jacobi algorithm for Singular Value Decomposition. This exam problem has the following problemformulation:
One-sided Jacobi algorithm for Singular Value Decomposition

Introduction
The singular value decomposition (SVD) of a (real square, for simplicity) matrix A is a representation of the matrix in the form
A = U D VT ,

where matrix D is diagonal with non-negative elements and matrices U and V are orghogonal. The diagonal elements of matrix D can always be chosen non-negative by multiplying the relevant columns of matrix U with (-1).
SVD can be used to solve a number of problems in linear algebra.

Problem
Implement the one-sided Jacobi SVD algorithm.
Algorithm
In this method the elementary iteration is given as
A → A J(θ,p,q)

where the indices (p,q) are swept cyclicly (p=1..n, q=p+1..n) and where the angle θ is chosen such that the columns number p and q of the matrix AJ(θ,p,q) are orthogonal. One can show that the angle should be taken from the following equation (you should use atan2 function),
tan(2θ)=2apTaq /(aqTaq - apTap)

where ai is the i-th column of matrix A (check this).
After the iterations converge and the matrix A'=AJ (where J is the accumulation of the individual rotations) has orthogonal columns, the SVD is simply given as

A=UDVT

where
V=J, Dii=||a'i||, ui=a'i/||a'i||,

where a'i is the i-th column of matrix A' and ui us the i-th column of matrix U.

