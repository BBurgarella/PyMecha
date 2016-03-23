# -*- coding: utf-8 -*-
"""
Pymecha was written by Boris Burgarella
the objective of this module is to simplify the micro-mechanics calculations using python

  Copyright (C) <2014>  <Boris Burgarella>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
   
   
"""

import numpy as np
from math import *
import ast

#Variables redondantes
matriceC =  np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],dtype = 'f')
I2 = [[1,0,0],[0,1,0],[0,0,1]]
MatriceC3x3x3x3 = [[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]],[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]],[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]]]


#Fonction utiles à l'interieur d'autres fonctions:

def testdelabase(El,Jt,F,tF,Kt,Kl):
	"""
	This function is made to test that the base used for projecting the tensors is right.

	Arguments: 3*3*3*3 Tensors El, Jt, F,tF,Kt,Kl
	Returns: Resultat, list of which length if 9. 1 means the equation is verified 0 means it isn't.

	If you want to use a base, remember that the whole Resultat list must be at 1 for the base to be usable.

	"""
	

	#Testing the different vectors relations that defines the base

	Resultat = [0,1,2,3,4,5,6,7,8,9,10]
 
	#Test 1 
 
	if checkequal0(produitContracte44(F,F)):
		Resultat[0] = 1
	else:
		Resultat[0] = 0
  
	#Test 2
  
	if checkequal44(produitContracte44(F,tF),Jt):
		Resultat[1] = 1
	else:
		Resultat[1] = 0  
  
	#Test 3
  
	if checkequal44(produitContracte44(tF,F),El):
		Resultat[2] = 1
	else:
		Resultat[2] = 0  
  
	#Test 4
  
	if checkequal0(produitContracte44(F,Kl)):
		Resultat[3] = 1
	else:
		Resultat[3] = 0  
  
	#Test 5
  
	if checkequal0(produitContracte44(F,Kt)):
		Resultat[4] = 1
	else:
		Resultat[4] = 0  
  
	#Test 6
  
	if checkequal0(produitContracte44(Kl,F)):
		Resultat[5] = 1
	else:
		Resultat[5] = 0
  
	#Test 7
  
	if checkequal0(produitContracte44(Kt,F)):
		Resultat[6] = 1
	else:
		Resultat[6] = 0
  
	#Test 8
  
	if checkequal44(produitContracte44(F,El),F):
		Resultat[7] = 1
	else:
		Resultat[7] = 0 
  
	#Test 9
  
	if checkequal44(produitContracte44(Jt,F),F):
		Resultat[8] = 1
	else:
		Resultat[8] = 0 
  
	#Test 10
  
	if checkequal0(produitContracte44(El,F)):
		Resultat[9] = 1
	else:
		Resultat[9] = 0 
  
	#Test 11
  
	if checkequal0(produitContracte44(F,Jt)):
		Resultat[10] = 1
	else:
		Resultat[10] = 0
	print Resultat
    

def delta(i,j): 

#symbole de kroeneker

	"""
	Kroeneker symbol.
	arguments: i and j
	type of arguments: int i and int j
	Returns: 1 if i == j and 0 elsewhat

	"""
	

	if j==i:
		return 1
	else:
		return 0

def matriceordreN(n):

	"""
	takes int n as input, and returns a 3*3*3*3*3 [...] n times tensor.

	example:
	matrice(2) = 
	[0,0,0]
	[0,0,0]
	[0,0,0]
					
	the returned tensor is full of 0

	"""
	
	if n == 1:
		matrice = [0,0,0]
	else:
		matrice = [matriceordreN(n-1),matriceordreN(n-1),matriceordreN(n-1)]

	return matrice
    
def passage24(a):

	"""
	this function is written to help the function that transform the 6x6 stiffness matrix back to it's normal form 3*3*3*3
	arguments: int a
    returns: int values: [i,j] vector: 	
	if a <= 2: i=j=a 
	if a = 3: i = 1 and j = 2
	if a = 4: i = 0 and j = 2
	if a = 5: i = 0 and j = 1
	This function is the inverse fuction of passage42(i,j)

	"""
	
	if a <= 2:
		i = a
		j = a
	elif a == 3:
		i = 1
		j = 2
	elif a == 4:
		i = 0
		j = 2
	elif a == 5:
		i = 0
		j = 1
	return [i,j]
    
    
def passage42(i,j):

	"""
	this function is written to help the function that transform the 6x6 stiffness matrix back to it's normal form 3*3*3*3
	arguments: int [i,j]
	returns: int valuea: [i,j] vector: 	
	if i=j: a = i = j 
	if i = 1 and j = 2: a = 3
	if i = 0 and j = 2: a = 4
	if i = 0 and j = 1: a = 5
	This function is the inverse fuction of passage24(a)

	"""
	
	if i == j:
		a = i
	elif i+j == 3:
		a = 3
	elif i+j == 2:
		a = 4
	elif i+j == 1:
		a = 5
	return a
    
########################################    
# Definition des opérateurs tensoriels #
########################################

#--------------- Ordre 2 ----------------------
def produittensoriel12(A,B):

	"""
	takes two vectors A and B (dimensions: 3x1)
	and returns the tensorial product of A and B

	Type: int or float
	returns: Matrix of float, dimension NxN

	Functions called by this one: matriceordreN(n)

	"""
	
	C = matriceordreN(2)
	for i in range(len(A)):
		for j in range(len(A)):
			ValA = A[i]
			ValB = B[j]
			C[i][j] = ValA*ValB
	return C
    
    
def produitMatScal2(Scal,Mat):
	"""
	simple product between a scalar and a matrix of dimension NxP
	type: float or int
	returns: Matrix of float, dimention NxP

	"""
	
	for i in range(len(Mat)):
		for j in range(len(Mat[i])):
			Mat[i][j] = Scal*Mat[i][j]

	return Mat

def soustrMatordre2(A,B):
	"""
	I wrote this function when I wasn't using numpy, and I needed a 
	function to substrat two matrix, it takes two 3x3 matrix as arguments
	and returns a 3x3 matrix equals to A-B

	Functions called by this one: matriceordreN(n)

	"""
	
	K = matriceordreN(2)
	for i in range(3):
		for j in range(3):
			K[i][j] = A[i][j]-B[i][j]
	return K
    
def addMat2(I,J):
	"""
	I wrote this function when I wasn't using numpy, and I needed a 
	function to add two matrix, it takes two 3x3 matrix as arguments
	and returns a 3x3 matrix equals to A+B

	Functions called by this one: matriceordreN(n)

	"""
	
	K = matriceordreN(2)
	for i in range(3):
		for j in range(3):
			K[i][j] = I[i][j]+J[i][j]
	return K
    
#--------------- Ordre 4 ---------------------- 
    
def transpose4(mat):
	"""
	Takes a 3x3x3x3 tensor as arguments

	returns the transposed tensor

	Functions called by this one: matriceordreN(n)

	"""
	
	tmat = matriceordreN(4)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					tmat[k][l][i][j] = mat[i][j][k][l]
	return tmat

def definitionI():

	"""
	this function determines the identity tensor for the symmetrical 4rth order tensor vectorial space
	arguments: null
	returns: 3x3x3x3 tensor I2

	Functions called by this one: matriceordreN(n)

	"""
	
	I = matriceordreN(4)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					I[i][j][k][l] = 0.5*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))
	return I    

def produittensoriel24(A,B):
	"""
	takes two vectors A and B (dimensions: 3x3)
	and returns the tensorial product of A and B

	Type: int or float
	returns: Matrix of float, dimension 3x3x3x3

	Functions called by this one: matriceordreN(n)
	"""
	
	C = matriceordreN(4)
	for i in range(len(A)):
		for j in range(len(A)):
			for k in range(len(A)):
				for l in range(len(A)):
					ValA = A[i][j]
					ValB = B[k][l]
					C[i][j][k][l] = ValA*ValB
	return C

    
def soustrMat(I,J):
	"""
	I wrote this function when I wasn't using numpy, and I needed a 
	function to substrat two matrix, it takes two 3x3x3x3 matrix as arguments
	and returns a 3x3x3x3 matrix equals to A-B

	Functions called by this one: matriceordreN(n)

	"""
	
	K = matriceordreN(4)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					K[i][j][k][l] = I[i][j][k][l]-J[i][j][k][l]
	return K
    
def addMat(I,J):
	"""
	I wrote this function when I wasn't using numpy, and I needed a 
	function to add two matrix, it takes two 3x3x3x3 matrix as arguments
	and returns a 3x3x3x3 matrix equals to A+B

	Functions called by this one: matriceordreN(n)

	"""
	K = matriceordreN(4)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					K[i][j][k][l] = I[i][j][k][l]+J[i][j][k][l]
	return K
    

def produitMatScal(Scal,Mat):
	"""
	simple product between a scalar and a matrix of dimension 3x3x3x3
	type: float or int
	returns: Matrix of float, dimention 3x3x3x3

	"""
	
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					Mat[i][j][k][l] = Scal*Mat[i][j][k][l]

	return Mat

def produitContracte42(A,a):
	"""
	This functions calculates the "contracted product" (name directly translated from french as I don't know the english name.
	this product is defined by this equations: Aa[i][j] = A[i][i][k][l]*a[k][l] (Einstein notation)

	arguments: a matrix of dimension 3x3 and A tensor of dimension 3x3x3x3
	returns Aa matrix of dimension 3x3

	"""
	
	Aa = matriceordreN(2)
	for i in range(3):
		for j in range(3):
			temp = 0
			for k in range(3):
				for l in range(3):
					temp = temp+(A[i][j][k][l]*a[k][l])
			Aa[i][j] = temp
	return Aa
    
def produitContracte44(A,B):
	"""
	This function calculates the "contracted product" (name directly translated from french as I don't know the english name.
	this product is defined by this equations: AB[i][j][k][l] = A[i][j][m][n]*B[m][n][k][l] (Einstein notation)

	arguments: a matrix of dimension 3x3x3x3 and A tensor of dimension 3x3x3x3
	returns AB matrix of dimension 3x3x3x3

	"""
	AB = matriceordreN(4)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					temp = 0
					for m in range(3):
						for n in range(3):
							temp = temp+(A[i][j][m][n]*B[m][n][k][l])
					AB[i][j][k][l] = temp
	return AB
	
def produitDoublementContracte44(A,B):
	"""
	This function calculates the "contracted product" (name directly translated from french as I don't know the english name.
	this product is defined by this equations: A::B = A[i][j][k][l]*B[i][j][k][l] (Einstein notation)

	arguments: a matrix of dimension 3x3x3x3 and A tensor of dimension 3x3x3x3
	returns AB real

	"""
	AB = 0
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					print AB
					AB = AB+(A[i][j][k][k]*B[i][j][k][l])
	return AB

# Changement d'ordre pour matrices symmétriques


def Voigt(matrice4):
	"""
	This function takes a 4rth order symmetrical tensor and returns it in Voigt notation
	it is the inverse function of UnVoigt(matrix)
	Arguments: 4rth order tensor of dimension 3x3x3x3
	returns: 6x6 matrix

	Functions called by this one: passage24(a)	

	"""
	
	matrice2 = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],dtype = 'f')
	for a in range(6):
		for b in range(6):
			ij = passage24(a)
			kl = passage24(b)
			i = ij[0]
			j = ij[1]
			k = kl[0]
			l = kl[1]
			matrice2[a][b] = matrice4[i][j][k][l]
	return matrice2
    
def UnVoigt(matrice2):
	"""
	This function takes a 6x6 matrix and returns it's equivalent 4rth order symmetrical tensor
	it is the inverse function of Voigt(Tensor)
	Arguments: 6x6 matrix
	Returns: 4rth order tensor of dimension 3x3x3x3


	Functions called by this one: passage42(i,j)	

	"""
	
	matrice4 = matriceordreN(4)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					a = passage42(i,j)
					b = passage42(k,l)
					matrice4[i][j][k][l] = matrice2[a][b]
	return matrice4

def SixXSix(Matrice4):
	"""
	This function takes a 4rth order symmetrical tensor and returns it in modified Voigt notation (with sqrt(2) and 2 multipled to some terms, which simplifies the strain and stress vector notation) 
	it is the inverse function of Un6x6(matrix)
	Arguments: 4rth order tensor of dimension 3x3x3x3
	returns: 6x6 matrix

	Functions called by this one: passage24(a)	

	"""
	
	matrice2 = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],dtype = 'f')
	for a in range(6):
	   for b in range(6):
			coefa = 1
			coefb = 1
			ij = passage24(a)
			kl = passage24(b)
			i = ij[0]
			j = ij[1]
			k = kl[0]
			l = kl[1]
			if a > 2:
				coefa = sqrt(2)
			if b > 2:
				coefb = sqrt(2)
			matrice2[a][b] = coefa*coefb*Matrice4[i][j][k][l]
	return matrice2
   
def Un6x6(Matrice2):
	"""
	This function takes a 6x6 matrix and returns it's equivalent 4rth order symmetrical tensor
	it is the inverse function of SixXSix(Tensor)
	Arguments: 6x6 matrix
	Returns: 4rth order tensor of dimension 3x3x3x3


	Functions called by this one: passage42(i,j)	

	"""
	
	matrice4 = matriceordreN(4)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					a = passage42(i,j)
					b = passage42(k,l)
					if a >= 3:
						facta = 1/sqrt(2)
					else:
						facta = 1
					if b >= 3:
						factb = 1/sqrt(2)
					else:
						factb = 1
					matrice4[i][j][k][l] = facta*factb*Matrice2[a][b]
	return matrice4
   

# Projection espace vectoriel J,K
def CalculJ():
	"""
	this function calculates J vector, this vector is part of the base of the symmetrical 4rth order isotropic tensor vectorial space
	J is a 3x3x3x3 tensor

	"""

	i = np.array([[1,0,0],[0,1,0],[0,0,1]])
	J = produittensoriel24(i,i)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					J[i][j][k][l] = (1.0/3.0)*J[i][j][k][l]
	return J


def CalculK(J,I):
	"""
	this function calculates K vector, this vector is part of the base of the symmetrical 4rth order isotropic tensor vectorial space
	k is a 3x3x3x3 tensor

	it takes as arguments J (calculated by CalculJ() function and I the Identity (of the 4rth order symmetrical tensors)

	"""
	
	K = matriceordreN(4)
	K = soustrMat(I,J)
	return K
    
# Projection espace vectoriel El,Jt,Kt,Kl,F, isotropie transverse

def CalculKe(n,I):
	"""
	this function calculates Ke vector, this vector is part of the base of the symmetrical 4rth order transverse isotropic tensor vectorial space
	ke is a 3x3x3x3 tensor

	it takes as arguments n, which is the Transverse isotropy axis (3x1 vector) and I the Identity (of the 4rth order symmetrical tensors)

	"""
	It2 = matriceordreN(2)
	nxn = produittensoriel12(n,n)
	for i in range(3):
		for j in range(3):
			It2[i][j] = I2[i][j]-nxn[i][j]
	Ke = matriceordreN(4)
	n2 = [0,0,0]
	for i in range(3):
		n2[i] = n[i]
	for i in range(3):
		n2[i] = 2*n2[i]
	Ketemp = soustrMatordre2(produittensoriel12(n2,n),It2)
	Ke = produittensoriel24(Ketemp,Ketemp)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					Ke[i][j][k][l] = (1./6)*(Ke[i][j][k][l])
	return Ke

def CalculJt(I,n):
	"""
	this function calculates Jt vector, this vector is part of the base of the symmetrical 4rth order transverse isotropic tensor vectorial space
	Jt is a 3x3x3x3 tensor

	it takes as arguments n, which is the Transverse isotropy axis (3x1 vector) and I the Identity (of the 4rth order symmetrical tensors)

	"""
	
	It2 = matriceordreN(2)
	nxn = produittensoriel12(n,n)
	for i in range(3):
		for j in range(3):
			It2[i][j] = I2[i][j]-nxn[i][j]
	Jt = matriceordreN(4)
	Jt = produittensoriel24(It2,It2)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					Jt[i][j][k][l] = (1./2)*(Jt[i][j][k][l])
	return Jt

def CalculIt(n):
	"""
	this function calculates It vector, this vector is part of the base of the symmetrical 4rth order transverse isotropic tensor vectorial space
	It is a 3x3x3x3 tensor

	it takes as arguments n, which is the Transverse isotropy axis (3x1 vector)

	"""
	
	I2 = matriceordreN(2)
	It = matriceordreN(4)
	for i in range(3):
		I2[i][i] = 1
	temp = produittensoriel12(n,n) 
	It = soustrMatordre2(I2,temp)
	return It

def CalculKt(Jt,It):
	"""
	this function calculates Kt vector, this vector is part of the base of the symmetrical 4rth order transverse isotropic tensor vectorial space
	Kt is a 3x3x3x3 tensor

	it takes as arguments J (calculated by CalculJt(I,n) function and It (calculated by CalculIt(n) function

	"""
	Kt = soustrMat(It,Jt)
	return Kt

def CalculF(It,n):
	"""
	this function calculates F vector, this vector is part of the base of the symmetrical 4rth order transverse isotropic tensor vectorial space
	F is a 3x3x3x3 tensor

	it takes as arguments n, which is the Transverse isotropy axis (3x1 vector) and It (calculated by CalculIt(n) function

	"""
	
	Itf = matriceordreN(2)
	for i in range(3):
		for j in range(3):
			Itf[i][j] = (1/sqrt(2.))*It[i][j]
	temp1 = produittensoriel12(n,n)
	F = produittensoriel24(Itf,temp1)
	return F
    

def CalculKl(K,Kt,Ke):
	"""
	this function calculates Kl vector, this vector is part of the base of the symmetrical 4rth order transverse isotropic tensor vectorial space
	Kl is a 3x3x3x3 tensor

	it takes as arguments K (calculated by CalculK(n) function, It (calculated by CalculKt(n) function and It (calculated by CalculKe(n) function

	"""
	
	Kl = matriceordreN(4)
	Kltemp = soustrMat(K,Kt)
	Kl = soustrMat(Kltemp,Ke)
	return Kl

def CalculEl(n):
	"""
	this function calculates F vector, this vector is part of the base of the symmetrical 4rth order transverse isotropic tensor vectorial space
	F is a 3x3x3x3 tensor

	it takes as arguments n, which is the Transverse isotropy axis (3x1 vector)

	"""
	temp1 = produittensoriel12(n,n)
	El = produittensoriel24(temp1,temp1)
	return El
    
def CalculIT(n):
	"""
	this function calculates IT, which is the 4rth order tensor identity restrained in the transverse plan (the normal plan to the n vector)

	it takes as arguments n, which is the Transverse isotropy axis (3x1 vector)
	and returns a 4rth order tensor 3x3x3x3

	"""
	
	IT = matriceordreN(4)
	IT2 = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
	for i in range(3):
		IT2[i][i]= 1 - n[i]
		IT2[i+3][i+3]= n[i]
	IT = Un6x6(IT2)
	return IT

    
def CoordExtract(A,Vecteur):
	"""
	this function takes as entry a tensor A and another named vector because it is supposed to be part of a vectorial space base.
	then it check if A can be written as N*Vector with N being a float.

	arguments: two 4rth order tensors
	returns: N if A = N*Vector
	or 0 if not

	"""
	
	ListeResult = []  
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					if A[i][j][k][l] != 0 and abs(A[i][j][k][l]) >= 1e-10 and abs(Vecteur[i][j][k][l]) >= 1e-10 :
						ListeResult.append(A[i][j][k][l]/Vecteur[i][j][k][l])
	if len(ListeResult)>=1:
		return ListeResult[0]
	else:
		return 0

def ProjectionKJ(C):
	"""
	this function takes as entry a vector n, then it generates El,JT,Kt,Kl and F (see the functions named CalculEl(), CalculJT, etc.) and a stiffnes (of compliance) tensor C
	and finally returns a list of scalar that are the coordinates of C in the base El Jt Kt etc.

	let's consider A = a*EL+b*Jt+g*F+gp*Ft+d*Kt+dp*Kl
	then this function will return:
	[A,d,dp] with A 2x2 matrix 	
	[ a , gp]
	[ g , b ]
	this format will be named projected format

	"""
	
	CoordC = [0,0,0,0,0,0]
	a = np.array([[0,0],[0,0]],dtype = 'f')
	A = []    
	I = definitionI()
	J = CalculJ()
	K = CalculK(J,I)
	mu = produitDoublementContracte44(J,C)/float(10)
	k = produitDoublementContracte44(K,C)/float(3)
	return [k,mu]
    
def ProjectionEl_Jt_Kt_Kl_F(C,n):
	"""
	this function takes as entry a vector n, then it generates El,JT,Kt,Kl and F (see the functions named CalculEl(), CalculJT, etc.) and a stiffnes (of compliance) tensor C
	and finally returns a list of scalar that are the coordinates of C in the base El Jt Kt etc.

	let's consider A = a*EL+b*Jt+g*F+gp*Ft+d*Kt+dp*Kl
	then this function will return:
	[A,d,dp] with A 2x2 matrix 	
	[ a , gp]
	[ g , b ]
	this format will be named projected format

	"""
	
	CoordC = [0,0,0,0,0,0]
	a = np.array([[0,0],[0,0]],dtype = 'f')
	A = []    
	I = definitionI()
	IT = CalculIT(n)
	J = CalculJ()
	K = CalculK(J,I)
	It = CalculIt(n)
	El = CalculEl(n)
	Jt = CalculJt(I,n)
	Ke = CalculKe(n,I)
	Kt = CalculKt(Jt,IT)
	F = CalculF(It,n)
	tF = transpose4(F)
	Kl = CalculKl(K,Kt,Ke)
	CoordC[0] = CoordExtract(produitContracte44(produitContracte44(El,C),El),El) #bon  
	CoordC[1] = CoordExtract(produitContracte44(produitContracte44(Jt,C),Jt),Jt) #bon
	CoordC[2] = CoordExtract(produitContracte44(produitContracte44(F,C),F),F) # bon
	CoordC[3] = CoordExtract(produitContracte44(produitContracte44(tF,C),tF),tF) #bon 
	CoordC[4] = CoordExtract(produitContracte44(produitContracte44(Kt,C),Kt),Kt)  
	CoordC[5] = CoordExtract(produitContracte44(produitContracte44(Kl,C),Kl),Kl)
	a[0][0] = CoordC[0]
	a[1][1] = CoordC[1]
	a[0][1] = CoordC[3]
	a[1][0] = CoordC[2]
	A.append(a)
	A.append(CoordC[4])
	A.append(CoordC[5])
	return A
    
def EngineerNotation(Vecteur,n,PrintOpt):
	""" 
		this function is usefull if you want to see the mechanical properties of a material which stiffness matrix is "Vecteur" and expressed in the El_Jt_Kt_Kl_F base 
		(see homogénéisation en mécanique des mateiraux 1 and 2 by M. Bornert T.Bretheau and P.Gilormini
		
		arguments: 	
		Vecteur => Stiffness tensor expressed in El_Jt_Kt_Kl_F base
		n => axis of transverse isotropy
		PrintOpt => bool type variable, defines if the function will only return El of print the properties
		return: El (this can be changed in the code to be another characteristic and an optional argument will be added in further versions 
		
	"""
	
	CoordC = [0,0,0,0,0,0]   
	a = Vecteur[0]
	CoordC[0] = a[0][0]
	CoordC[1] = a[1][1]
	CoordC[3] = a[0][1]
	CoordC[2] = a[1][0]
	CoordC[4] = Vecteur[1]
	CoordC[5] = Vecteur[2]
	El = CoordC[0]-((CoordC[2]**2)/(CoordC[1]))
	K = CoordC[1]/2
	mut = CoordC[4]/2
	deltap = CoordC[5]/2
	# Test
	Coef1 = (CoordC[0]+CoordC[5])/2.
	Coef2 = (CoordC[0]-CoordC[5])/2.
	Coef3 = CoordC[2]/sqrt(2)
	Coef4 = CoordC[3]/sqrt(2)
	Ctest = np.array([[Coef1,Coef2,Coef3,0,0,0],[Coef2,Coef1,Coef3,0,0,0],[Coef4,Coef4,CoordC[0],0,0,0],[0,0,0,CoordC[5],0,0],[0,0,0,0,CoordC[5],0],[0,0,0,0,0,CoordC[4]]])    
	Stest = np.linalg.inv(Ctest) 
	Et = 1/Stest[1][1]
	El = 1/Stest[0][0]
	#fin Test
	if PrintOpt == 1:
		print 'Et = ', Et 
		print 'K = ', K
		print 'mut = ', mut
		print 'Et = ', El
		print 'G = ', deltap
	return El

def unprojectEl_Jt_Kt_Kl_F(Vecteur,n):
	"""
	this function takes as argument a projected tensor (see ProjectionEl_Jt_Kt_Kl_F(C,n) for more information) and n the transverse isotropy axis
	and it returns the non-projected tensor (dimension 3x3x3x3)

	dependencies: some function of this module are called by unprojectEl_Jt_Kt_Kl_F():
	matriceordreN()
	all the CalculXX functions
	addMat()

	"""

	CoordC = [0,0,0,0,0,0]   
	a = Vecteur[0]
	CoordC[0] = a[0][0]
	CoordC[1] = a[1][1]
	CoordC[3] = a[0][1]
	CoordC[2] = a[1][0]
	CoordC[4] = Vecteur[1]
	CoordC[5] = Vecteur[2]
	C = matriceordreN(4)
	I = definitionI()
	IT = CalculIT(n)
	J = CalculJ()
	K = CalculK(J,I)
	It = CalculIt(n)
	El = CalculEl(n)
	Jt = CalculJt(I,n)
	Ke = CalculKe(n,I)
	Kt = CalculKt(Jt,IT)
	F = CalculF(It,n)
	tF = transpose4(F)
	Kl = CalculKl(K,Kt,Ke)
	Base = [El,Jt,F,tF,Kt,Kl]
	for i in range(len(Base)):
		C = addMat(C,produitMatScal(CoordC[i],Base[i]))
	return C
    
def ProduitIsoTrans(A,B):
	"""
	this function takes as arguments two tensors in projected format, and returns the tensorial product of these two.

	"""
	
	C = []
	C.append(np.dot(A[0],B[0]))
	C.append(A[1]*B[1])
	C.append(A[2]*B[2])
	return C

def SoustrIsoTrans(A,B):
	"""
	this function takes as arguments two tensors in projected format, and returns the substrat of these two (A-B).
	"""
	C = []
	C.append(np.subtract(A[0],B[0]))    
	C.append(A[1]-B[1])
	C.append(A[2]-B[2])
	return C

def SommeIsoTrans(A,B):
	"""
	this function takes as arguments two tensors in projected format, and returns the sum of these two (A+B).

	"""
	
	C = []
	C.append(np.add(A[0],B[0]))    
	C.append(A[1]+B[1])
	C.append(A[2]+B[2])
	return C

def InverseIsoTrans(A):
	"""
	this function takes as argument a tensor A in projected format, and returns the inverse of it (still in projected format)

	"""

	C = []
	C.append(np.linalg.inv(A[0]))
	C.append(1/(A[1]))
	C.append(1/(A[2]))
	return C
    
def ProduitScalVectIsoTrans(Scal,Vect):
	"""
	this function takes as argument a tensor A in projected format and a scalar.
	It returns the the tensor multiplied by the scalar.

	"""
	
	C = []
	C.append(Scal*Vect[0])
	C.append(Vect[1]*Scal)
	C.append(Vect[2]*Scal)
	return C
# Fonctions Mecanique

def EngValuesFromC(C):
	"""
	this function takes as argument a stiffness tensor C (3x3x3x3), and returns a list of it's mechanical properties (with a degree of inhomogeneity up to orthotropic) 
 
	"""
	C2 = SixXSix(C)
	Matrice = np.linalg.inv(C2)
	Ex = 1/Matrice[0][0]
	Ey = 1/Matrice[1][1]
	Ez = 1/Matrice[2][2]
	nuyx = -Matrice[0][1]*Ey
	nuzx = -Matrice[0][2]*Ez
	nuzy = -Matrice[1][2]*Ez
	Gyz = 2*Matrice[3][3]
	Gzx = 2*Matrice[4][4]
	Gxy = 2*Matrice[5][5]
	return [Ex,Ey,Ez,nuyx,nuzx,nuzy,Gyz,Gzx,Gxy] 
    
 

def RemplissageCval(C11,C22,C33,C12,C13,C23,C44,C55,C66):
	"""
	This function takes as arguments the mechanical properties of an orthotropic material (it works for materials with lower degree of inhomogeneity)
	and returns the corresponding stiffness tensor (3x3x3x3)

	!! This function require numpy to be installed on you computer to work properly !!

	"""
	
	C = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],dtype = 'f')
	C[0][0] = C11
	C[0][1] = C12
	C[0][2] = C13
	C[1][0] = C[0][1]
	C[1][1] = C22
	C[1][2] = C23
	C[2][0] = C13
	C[2][1] = C[1][2]
	C[2][2] = C33
	C[3][3] = C44
	C[4][4] = C55
	C[5][5] = C66
	print C
	return C

def RemplissageC(Ex,Ey,Ez,nuyx,nuzx,nuzy,Gyz,Gzx,Gxy):
	"""
	This function takes as arguments the mechanical properties of an orthotropic material (it works for materials with lower degree of inhomogeneity)
	and returns the corresponding stiffness tensor (3x3x3x3)

	!! This function require numpy to be installed on you computer to work properly !!

	"""
	
	Ex = float(Ex)
	Ey = float(Ey)
	Ez = float(Ez)
	C = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],dtype = 'f')
	C[0][0] = 1/Ex
	C[0][1] = -nuyx/Ey
	C[0][2] = -nuzx/Ez
	C[1][0] = C[0][1]
	C[1][1] = 1/Ey
	C[1][2] = -nuzy/Ez
	C[2][0] = C[0][2]
	C[2][1] = C[1][2]
	C[2][2] = 1/Ez
	C[3][3] = 1/(2*Gyz)
	C[4][4] = 1/(2*Gzx)
	C[5][5] = 1/(2*Gxy)
	print C
	Matrice = np.linalg.inv(C)
	return Matrice

# Fonctions d'affichage

def printVoigt(Mat):
	"""
	A simple print function that print tensor in Voigt notation (6x6 matrix)

	arguments: 3x3x3x3 Tensor
	returns: Null

	"""
	
	print Voigt(Mat) 

def print6x6(Mat):
	"""
	A simple print function that print tensor in modified Voigt notation (6x6 matrix)

	arguments: 3x3x3x3 Tensor
	returns: Null

	"""
	print SixXSix(Mat)        
    
  
def checkequal44(A,B):
	"""
	this function was deigned to check the equality between two 4rth order tensors
	it returns a boolean which is set to true if A = B and False else-what.

	"""
	
	Test = 1
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					val1 = A[i][j][k][l]
					val2 = B[i][j][k][l]
				if abs(val1-val2) >= 0.1:
					Test = 0
					return False
	if Test  == 1:
		return True
   
def checkequal0(A):
	"""
	Booléan function that takes a tensor A (3x3x3x3) as entry and returns True if A = 0 and False is A isn't null

	"""

	Test = 1
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					val1 = A[i][j][k][l]
					val2 = 0
				if abs(val1-val2) >= 0.1:
					Test = 0
					return False
	if Test  == 1:
		return True

def checkequal(A,B):
	"""
	Booléan function that takes two 6x6 matrices A and B as entry and returns True if A = B and False if not

	"""
	
	Test = 1
	for i in range(6):
		for j in range(6):
			val1 = A[i][j]
			val2 = B[i][j]
			if val1 != val2:
				Test = 0
				return False
	if Test  == 1:
		return True

# Homogenéisation generale

def NumSup(A):
	"""
	Basic function that returns the upper value of a list (float of int)

	"""
	NumValSup = 0
	for i in range(len(A)):
		if A[i] > A[NumValSup]:
			NumValSup = i
	return NumValSup

def NumInf(A):
	"""
	Basic function that returns the lower value of a list (float of int)

	"""
	
	NumValInf = 0
	for i in range(len(A)):
		if A[i] < A[NumValInf]:
			NumValInf = i
	return NumValInf

def Lielens(MTplus,MTmoins,taux):
	"""
	the Lielens homogenization model is a sort of mean between the two Mori-tanaka bounds, this functions simply takes the two bounds and makes this mean.
	the volume ratio is also needed as second argument.

	"""
	FacteurLiel = (taux[0]+taux[0]**2)/2    
	Cliel = SommeIsoTrans(ProduitScalVectIsoTrans((1-FacteurLiel),MTmoins),ProduitScalVectIsoTrans((FacteurLiel),MTplus))
	return Cliel

def ProjectTableMatIsoTrans(TableMat,n):
	"""
	takes a list of stiffness tensors as entry and the fiber/inclusion axis n 
	and returns a list of these tensor in projected version?

	"""
	
	TableMatprojetee = []
	for i in range(len(TableMat)):
		C = RemplissageC(TableMat[i][0],TableMat[i][1],TableMat[i][2],TableMat[i][3],TableMat[i][4],TableMat[i][5],TableMat[i][6],TableMat[i][7],TableMat[i][8])       
		C4 = Un6x6(C)
		Vecteur = ProjectionEl_Jt_Kt_Kl_F(C4,n)
		TableMatprojetee.append(Vecteur)
	return TableMatprojetee

def ModuleTableGenerator(TableMat):
	"""
	Takes a list of projected stiffness tensors and returns a list of their stiffness module (works for homogeneous materials)

	"""
	Module = []
	for i in range(len(TableMat)):
		Module.append((TableMat[i][0]+TableMat[i][1]+TableMat[i][2])/3)
	return Module
    
def VoigtBound(TableMat,Taux):
	"""
	This function do the homogenization of an heterogeneous material using the voigt homogenization model
	arguments:
	TalbeMat => list of projected stiffness tensors
	Taux: float representing the volume ratio between the two materials
	n: the fiber axis vector

	returns:
	a list of two projected tensors which are the sup and inf bound of Mori Tanaka method.

	"""
	
	Result = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],dtype = 'f')
	for i in range(len(TableMat)):
		if TableMat[i][0] == 0:
			C = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],dtype = 'f')
		else:
		   C = RemplissageC(TableMat[i][0],TableMat[i][1],TableMat[i][2],TableMat[i][3],TableMat[i][4],TableMat[i][5],TableMat[i][6],TableMat[i][7],TableMat[i][8])
		Result += Taux[i]*C
	return Result

def importTensorMathematica(fileName):
	"""
	This function import a tensor exported from mathematica in a text file.
	argument:
	filename, string containing the filename and it's extension

	returns: Table representing the tensor

	"""
	
	file = open(fileName,'r')
	Liste = file.read()
	i = 0
	PythonInstruct = ''
	while i < (len(Liste)):
		if Liste[i] == '{':
			PythonInstruct += '['
		elif Liste[i] == '}':
			PythonInstruct += ']' 
		elif Liste[i] == '^':
			PythonInstruct += '10**' 
		else:
			PythonInstruct += Liste[i]
		i += 1
	Tensor = eval(PythonInstruct)
	return Tensor
    
def MoriTanakaLongFiber(TableMat,Taux,n):
	"""
	This function do the homogenization of an heterogeneous material with cylindrical inclusion
	arguments:
	TalbeMat => list of projected stiffness tensors
	Taux: float representing the volume ratio between the two materials
	n: the fiber axis vector

	returns:
	a list of two projected tensors which are the sup and inf bound of Mori Tanaka method.

	"""
	
	TableMatprojetee = []
	Module = []
	CMTplus = [[[0,0],[0,0]],0,0]
	CMTmoins = [[[0,0],[0,0]],0,0]
	TableMatprojetee = ProjectTableMatIsoTrans(TableMat,n)
	Module = ModuleTableGenerator(TableMat)     
	 # Borne HS+
	PartC = [np.array([[0,0],[0,0]]),0,0]
	NumMatAmbiant = NumSup(Module)
	MatAmbiant = TableMatprojetee[NumMatAmbiant]  
	K0 = MatAmbiant[0][1][1]/2
	Mu0 = MatAmbiant[1]/2
	Cstar = [np.array([[0,0],[0,2*Mu0]]),2*((K0*Mu0)/(K0+2*Mu0)),2*Mu0]
	for i in range(len(TableMat)):
		PartC = SommeIsoTrans(ProduitScalVectIsoTrans(Taux[i],InverseIsoTrans(SommeIsoTrans(Cstar,TableMatprojetee[i]))),PartC)
	CMTplus = SoustrIsoTrans(InverseIsoTrans(PartC),Cstar)



	# Borne HS-
	PartC = [np.array([[0,0],[0,0]]),0,0]
	NumMatAmbiant = NumInf(Module)
	MatAmbiant = TableMatprojetee[NumMatAmbiant]  
	K0 = MatAmbiant[0][1][1]/2
	Mu0 = MatAmbiant[1]/2
	Cstar = [np.array([[0,0],[0,2*Mu0]]),2*((K0*Mu0)/(K0+2*Mu0)),2*Mu0]
	for i in range(len(TableMat)):
		PartC = SommeIsoTrans(ProduitScalVectIsoTrans(Taux[i],InverseIsoTrans(SommeIsoTrans(Cstar,TableMatprojetee[i]))),PartC)
	CMTmoins = SoustrIsoTrans(InverseIsoTrans(PartC),Cstar)
	return [CMTplus,CMTmoins]

def MoriTanakaHillTensor(HillTensor,TableMat,Taux,n):
	"""
	This function do the homogenization of an heterogeneous material with cylindrical inclusion
	arguments:
	Hilltensor => Hill influence tensor (projected form)
	TalbeMat => list of projected stiffness tensors
	Taux: float representing the volume ratio between the two materials
	n: the fiber axis vector

	returns:
	a list of two projected tensors which are the sup and inf bound of Mori Tanaka method.

	"""
	
	TableMatprojetee = []
	Module = []
	CMTplus = [[[0,0],[0,0]],0,0]
	CMTmoins = [[[0,0],[0,0]],0,0]
	TableMatprojetee = ProjectTableMatIsoTrans(TableMat,n)
	Module = ModuleTableGenerator(TableMat)     
	 # Borne HS+
	PartC = [np.array([[0,0],[0,0]]),0,0]
	NumMatAmbiant = NumSup(Module)
	MatAmbiant = TableMatprojetee[NumMatAmbiant]  
	P = HillTensor
	invP = InverseIsoTrans(P)
	Cstar = SoustrIsoTrans(invP,MatAmbiant)
	for i in range(len(TableMat)):
		PartC = SommeIsoTrans(ProduitScalVectIsoTrans(Taux[i],InverseIsoTrans(SommeIsoTrans(Cstar,TableMatprojetee[i]))),PartC)
	CMTplus = SoustrIsoTrans(InverseIsoTrans(PartC),Cstar)



	# Borne HS-
	PartC = [np.array([[0,0],[0,0]]),0,0]
	NumMatAmbiant = NumInf(Module)
	MatAmbiant = TableMatprojetee[NumMatAmbiant]  
	P = HillTensor
	invP = InverseIsoTrans(P)
	Cstar = SoustrIsoTrans(invP,MatAmbiant)
	for i in range(len(TableMat)):
		PartC = SommeIsoTrans(ProduitScalVectIsoTrans(Taux[i],InverseIsoTrans(SommeIsoTrans(Cstar,TableMatprojetee[i]))),PartC)
	CMTmoins = SoustrIsoTrans(InverseIsoTrans(PartC),Cstar)
	return [CMTplus,CMTmoins]

def LielensLongFiberFromProjected(Mat1,Mat2,Taux):
	"""
	This function takes as entry two stiffness tensors and a volume ratio between the two materials
	and returns the Lielens method homogenized tensor

	This function works only for cylindrical inclusions, see MoriTanakaHillTensor() to homogenize inclusion with other shape.
	in further version, Lielens homogenization function with custom hill tensor / Eshelby tensor will be implemented. 

	"""
	
	TableMatprojetee = [Mat1,Mat2]

	# Borne HS+
	PartC = [np.array([[0,0],[0,0]]),0,0]
	NumMatAmbiant = 0
	MatAmbiant = TableMatprojetee[NumMatAmbiant]  
	K0 = MatAmbiant[0][1][1]/2
	Mu0 = MatAmbiant[1]/2
	Cstar = [np.array([[0,0],[0,2*Mu0]]),2*((K0*Mu0)/(K0+2*Mu0)),2*Mu0]
	for i in range(len(TableMatprojetee)):
		PartC = SommeIsoTrans(ProduitScalVectIsoTrans(Taux[i],InverseIsoTrans(SommeIsoTrans(Cstar,TableMatprojetee[i]))),PartC)
	CMTplus = SoustrIsoTrans(InverseIsoTrans(PartC),Cstar)

	# Borne HS-
	PartC = [np.array([[0,0],[0,0]]),0,0]
	NumMatAmbiant = 1
	MatAmbiant = TableMatprojetee[NumMatAmbiant]  
	K0 = MatAmbiant[0][1][1]/2
	Mu0 = MatAmbiant[1]/2
	Cstar = [np.array([[0,0],[0,2*Mu0]]),2*((K0*Mu0)/(K0+2*Mu0)),2*Mu0]
	for i in range(len(TableMatprojetee)):
		PartC = SommeIsoTrans(ProduitScalVectIsoTrans(Taux[i],InverseIsoTrans(SommeIsoTrans(Cstar,TableMatprojetee[i]))),PartC)
	CMTmoins = SoustrIsoTrans(InverseIsoTrans(PartC),Cstar)

	Cliel = Lielens(CMTplus,CMTmoins,Taux)
	return Cliel

def IsoTransCfromsigmaEpsilon(sigma,epsilon):
	"""
	this function calculates the stiffness tensor from the strain and stress.

	Arguments:
	6x1 lists Sigma (stress) and epsilon (Strain)

	returns :
	6x6 matrix (Stiffness matrix in modified voigt notation)

	"""
	
	s11 = sigma[0]
	s22 = sigma[1]
	s33 = sigma[2]
	s12 = sigma[3]
	Matrice = np.array([[s11,0,s22,s33],[s22,0,s11,s33],[0,s33,0,s11+s22],[s12,0,-s12,0]])
	Vecteur = np.array([epsilon[0],epsilon[1],epsilon[2],epsilon[3]])
	MatriceInv = np.linalg.inv(Matrice)
	VectSol = np.dot(MatriceInv,Vecteur)
	VectReturn = []
	for i in VectSol:
		VectReturn.append(i)
	VectReturn.append(epsilon[4]/sigma[4])
	a = VectReturn[0]
	b = VectReturn[1]
	c = VectReturn[2]
	d = VectReturn[3]
	e = VectReturn[4]
	MatriceC = np.array([[a,c,d,0,0,0],[c,a,d,0,0,0],[d,d,b,0,0,0],[0,0,0,a-c,0,0],[0,0,0,0,e,0],[0,0,0,0,0,e]])
	return MatriceC 

    