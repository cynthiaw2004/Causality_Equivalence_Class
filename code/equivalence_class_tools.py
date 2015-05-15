import sys
import os
scriptpath = "/na/homes/cfreeman/Documents/virtual_test_environment/gunfolds/tools"
sys.path.append(os.path.abspath(scriptpath))
import traversal, bfutils, graphkit,unknownrate,comparison
from Levenshtein import hamming
import numpy as np
from numpy.random import randint
import zickle as zkl
import matplotlib.pyplot as plt
from itertools import combinations



#this function is sort of bad and will be made better later after a pattern is found
#input: number of nodes
#output: try to create an H that will have a nonempty equivalence class
def generate_H(num_nodes):
	numextraedge = np.random.randint(low = 1,high = 3) 
	g = bfutils.ringmore(num_nodes,numextraedge) #ground truth 
	gs= bfutils.call_undersamples(g) # all possible undersamples for g
	randomu = np.random.randint(low = 1,high = len(gs))
	#now we can make the H
	H = bfutils.undersample(g,randomu)
	#print H
	return H

#this function is sort of bad and will be made better later after a pattern is found
#input: numH (number of H wanted)
def create_H_and_eqcs(numH):
	H_array = []
	eqc_array = []
	for i in range(0,numH):
		num_nodes = np.random.randint(low = 3, high = 7)
		H = generate_H(num_nodes)
		H_array.append(H)
		eqc_for_H = determine_equivalence_class_iterative_precompute(H)
		eqc_array.append(eqc_for_H)
	zkl.save(H_array,"random_H.zkl")
	zkl.save(eqc_array,"eqcs_for_H")

#input: H
#output: the equivalence class for H 
#NOTE: equivalence class consists of all ground truths that lead to the same H regardless of undersampling rate
#NOTE: to convert these ground truths back to dictionary notation use num2CG in comparison.py
def determine_equivalence_class_iterative_precompute(H):
	return unknownrate.liteqclass(H)

#input:G,u
#output: G^(u+1)
#note: G = G^1 = bfutils.undersample(G,0)
def determine_undersample(G,u):
	return bfutils.undersample(G,u+1)

#input: a graph G (dictionary format)
#output: a dictionary where key = vertex and value = in degree
def determine_in_degrees(G):
	in_degree_dict = {}
	for i in range(len(G)):
		in_degree_dict[str(i+1)] = 0
	for outerkey in G.keys():
		#print "\n"
		#print "from:",outerkey
		for innerkey in G[outerkey].keys():
			#print "to:",innerkey
			#print G[outerkey][innerkey]
			if (0,1) in G[outerkey][innerkey]:
				in_degree_dict[innerkey] = in_degree_dict[innerkey] + 1	
			if (2,0) in G[outerkey][innerkey]: #take into account bidirected edges
				in_degree_dict[outerkey] = in_degree_dict[outerkey] + 1
				in_degree_dict[innerkey] = in_degree_dict[innerkey] + 1
	return in_degree_dict

#input: a graph G (dictionary format)
#output: a dictionary where key = vertex and value =  out degree
def determine_out_degrees(G):
	out_degree_dict = {}
	for i in range(len(G)):
		out_degree_dict[str(i+1)] = 0
	for outerkey in G.keys():
		for innerkey in G[outerkey].keys():
			if (0,1) in G[outerkey][innerkey]:
				out_degree_dict[outerkey] = out_degree_dict[outerkey] + 1	
			if (2,0) in G[outerkey][innerkey]: #take into account bidirected edges
				out_degree_dict[outerkey] = out_degree_dict[outerkey] + 1
				out_degree_dict[innerkey] = out_degree_dict[innerkey] + 1
	return out_degree_dict

#warning: 
#graph2str PRESERVES bidirected edges when converting from dictionary to string
#but there is no function to convert back to dictionary
#graph2str safe to use with undersampled graphs 
#g2num does NOT preserve bidirected edges when converting to string
#but there IS a function to convert back to dictionary (num2CG)
#g2num safe to use with ground truths since ground truths do not have bidirected edges
#TO BE CHANGED SOON WHEN SERGEY PUSHES CHANGES TO GITHUB
def graph2str(G):
    n = len(G)
    d = {((0,1),):'1', ((2,0),):'0',((2,0),(0,1),):'1',((0,1),(2,0),):'1'}
    A = ['0']*(n*n)
    for v in G:
        for w in G[v]:
            A[n*(int(v,10)-1)+int(w,10)-1] = d[tuple(G[v][w])]
    return ''.join(A)

#input: two graphs G1 and G2 (must be dictionary format)
#ouput: the hamming distance (the number of differing characters) between the two graphs
def determine_edit_distance(G1,G2):
	G1str = str(graph2str(G1))
	G2str = str(graph2str(G2))
	#G1str and G2str MUST have same length
	return hamming(G1str,G2str)

#input: H and its equivalence class
#output: an array consisting of all the indegree dictionaries
def obtain_in_degrees_from_eqc(H,eqc):
	eqcgraphs = []
	for graphstr in eqc:
		graph = comparison.num2CG(graphstr,len(H))
		eqcgraphs.append(graph)
	indegrees = []
	for graph in eqcgraphs:
		indegrees.append(determine_in_degrees(graph))
	return indegrees

#input: H and its equivalence class
#output: an array consisting of all the outdegree dictionaries
def obtain_out_degrees_from_eqc(H,eqc):
	eqcgraphs = []
	for graphstr in eqc:
		graph = comparison.num2CG(graphstr,len(H))
		eqcgraphs.append(graph)
	outdegrees = []
	for graph in eqcgraphs:
		outdegrees.append(determine_out_degrees(graph))
	return outdegrees

#input: H and its equivalence class
#output: an array of the hamming distances between members in the eqc
def obtain_edit_distance_from_eqc_only(H,eqc):
	eqcgraphs = []
	for graphstr in eqc:
		graph = comparison.num2CG(graphstr,len(H))
		eqcgraphs.append(graph)
	hamming_distance = []
	for comb in combinations(eqcgraphs,2):
		hamming_distance.append(determine_edit_distance(comb[0],comb[1]))
	return hamming_distance

#input: H and its equivalence class
#output: an array of the hamming distances between members in the eqc and H
def obtain_edit_distance_from_eqc_and_H(H,eqc):
	eqcgraphs = []
	for graphstr in eqc:
		graph = comparison.num2CG(graphstr,len(H))
		eqcgraphs.append(graph)
	hamming_distance = []
	for graph in eqcgraphs:
		hamming_distance.append(determine_edit_distance(graph,H))
	return hamming_distance

#input: a ground truth (no bidirected edges allowed) G in dictionary format
#output: adjacency matrix
def determine_adjacency_matrix(G):
	adj_matrix = np.zeros((len(G),len(G)))
	for outerkey in G.keys():
		#print "from: ",outerkey
		for innerkey in G[outerkey].keys():
			#print "to: ",innerkey
			if (0,1) in G[outerkey][innerkey]:
				adj_matrix[int(outerkey)-1][int(innerkey)-1] = 1
	return adj_matrix			

#input: a matrix
#output: the eigenvalues
def determine_eigenvalues(M):
	return np.linalg.eig(M)


#############testing area###########################


#what doesnt seem like there isnt a pattern: in and out degree, hamming distance between members and between members and H, eigenvalues
#consider monoids now?


Hs = zkl.load("random_H.zkl")
eqcs = zkl.load("eqcs_for_H")

for i in range(len(Hs)):
	H = Hs[i]
	eqc = eqcs[i]
	if len(eqc) == 1 and eqc != set([-1]):
		print "H: ",H
		print eqc
		for graphstr in eqc:
			#print "eqc member: ",graphstr
			graph = comparison.num2CG(graphstr,len(H))
			print graph
		print "\n"





