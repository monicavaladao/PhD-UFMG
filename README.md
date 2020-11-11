# Proposta de Um Algoritmo Evolutivo Assistido por um Modelo de Aproximação Kriging para Problemas de Otimização de Alto Custo Computacional

> **Contributors:** Mônica A. C. Valadão<sup>1,2,4</sup>, Lucas S. Batista<sup>3,4</sup>  
> <sup>1</sup> *Graduate Program in Electrical Engineering - Universidade Federal de Minas Gerais ([url](https://www.ppgee.ufmg.br/))*  
> <sup>2</sup> *Science and Technology Institute - Universidade Federal dos Vales do Jequitinhonha e Mucuri ([url](http://ufvjm.edu.br/))*   
> <sup>3</sup> *Dept. Electrical Engineering - Universidade Federal de Minas Gerais ([url](http://www.dee.ufmg.br/))*  
> <sup>4</sup> *Operations Research and Complex Systems Lab. - Universidade Federal de Minas Gerais ([url](http://orcslab.ppgee.ufmg.br/))*


# About this repository

This repository contains the source code of the PhD dissertation entitled "Proposta de Um Algoritmo Evolutivo Assistido por um Modelo de Aproximação Kriging para Problemas de Otimização de Alto Custo Computacional", written by Mônica A. C. Valadão and Lucas S. Batista, submitted to the *Graduate Program in Electrical Engineering - Universidade Federal de Minas Gerais* ([url](https://www.ppgee.ufmg.br/)).

# About this work 

Optimization problems that require the evaluation of functions with high computational cost are frequently solved through metamodel-based strategies. Examples of strategies based on metamodels are the Surrogate Model Assisted Evolutionary Algorithms (SAEAs) that are usually employed to solve optimization problems that are computationally expensive to be evaluated and require several function evaluations, such as the ones with a large number of variables.

Currently, SAEAs have been applied in problems involving up to 200 variables. In such methods, the metamodel is used to guide the evolutionary algorithm towards promising regions of the search space and to reduce the number of function evaluations required. However, the cost associated with the update of the metamodel cannot be prohibitive. This work investigates the construction of global approximation models by Kriging and RBF metamodels. Besides presenting the main SAEAs found in the literature and different approaches to embed the metamodel on them, a new SAEA is proposed, named SAEAa, which couples in the same framework a parameter self-adaptation and a  mechanism that allows the choice between different mutation operators. More precisely, it couples mutation operators with distinct features, adding in the SAEAa maintenance of population diversity and selective pressure in the evolutive process. Another feature of SAEAa is that it employs a unidimensional Ordinary Kriging metamodel. Thus, it reduces the computational cost of training this kind of metamodel compared to the standard form Kriging. 

The description of the proposed strategy also addresses aspects that directly influence the quality of metamodel and population diversity, which are not treated in SAEAs existing in the literature. The proposed approach was employed to solve a set of analytical functions of single-objective optimization problems. The results obtained suggest that the SAEAa presents a better performance, in terms of solution quality and computational cost, when compared to recent strategies. Besides, the proposed approach was employed on the solution of a ground-penetrating radar (GPR) antenna design. The solution returned by SAEAa was validated, and it showed to be suitable for application in GPR. Furthermore, it was possible to show a considerable reduction in computational resources (time) from using the proposed approach.
