# Python3 program for Bellman-Ford's single source
# shortest path algorithm.

# Class to represent a graph
from importlib.metadata import entry_points
from operator import index
from platform import node
from pydoc import source_synopsis
from selectors import SelectorKey
from traceback import print_list
from turtle import clear
from numpy import size
import numpy


class Graph:
	def __init__(self, nodes,sources,destinations):
		self.graph = []					# The Graph restore many (u,v,w): the one-way edge from u to v, and the edge's weight is w (normally it is 0)


		self.V = nodes 					# Number of vertices
										# the nodes includes the nodes for pooling, sources and destinations

		self.sources=sources			# Nodes for inputing
		self.destinations=destinations	# Nodes for outputing
	
		self.in_arrows=[[]]*nodes		# in_arrow[i] means the node adjacent to Source i, e.g.in_arrow[6]=[3,5] means Source 6 is adjacent to node 3 and 5
		self.out_arrows=[[]]*nodes		# out_arrow[i] means the node adjacent to Destination i, e.g.out_arrow[6]=[3,5] means Destination 6 is adjacent to node 3 and 5

		self.terminals=[0]*nodes		# terminals[A]=the destination of Source A or the source of Destination A; 
										# A must be a destination of source
										# if A is neither a destination nor a source, then terminals[A]=0

		self.node_cost=[1]*nodes		# The costs of nodes, default value is 1
										# The node cost of a source or a destination is 0
		for i in range(size(sources)):

			# Define the destination of Source i and the source of Destination i
			self.terminals[sources[i]]=destinations[i]
			self.terminals[destinations[i]]=sources[i]

			#Define the node_cost
			self.node_cost[sources[i]]=0
			self.node_cost[destinations[i]]=0

	# function to add an edge to graph
	def addEdgeWithWeight(self, u, v, w):
		self.graph.append([u, v, w])	# the edge from u to v, w stands for the node number between u and v
		self.graph.append([v, u, w])	# the edge from v to u, w stands for the node number between u and v

	# function to add an edge with weight 0 to graph
	def addEdge(self, u, v):
		self.graph.append([u, v, 0])	# the edge from u to v
		self.graph.append([v, u, 0])	# the edge from v to u

	# function to add an edge to isolate two nodes, the weight of edge is a big number
	def addIsolation(self, u, v):
		self.graph.append([u, v, 1000000])	# the isolation from u to v
		self.graph.append([v, u, 1000000])	# the isolation from v to u

	# function to add edges with weight 0 to graph
	# addEdges(self, u, v) will help to add two-way edge (u,u+1),(u+1,u+2),(u+2,u+3)...(v-1,v) to the graph if v-u>1
	def addEdges(self, u, v):
		for i in range(u,v):	
			self.addEdge(i,i+1)
	
	# add the arrows from u to nodes of v
	# e.g. addInArrows(self,3,[5,6]) means add one-way edge (3,5),(3,6) to the graph
	def addInArrows(self,u,v):
		self.in_arrows[u]=v
		for j in v:	
			self.graph.append([u,j,0])

	# add the arrows nodes of u to v
	# e.g. addOutArrows(self,[5,6],3) means add one-way edge (5,3),(6,3) to the graph
	def addOutArrows(self,u,v):
		self.out_arrows[v]=u
		for i in u:	
			self.graph.append([i,v,0])

	def BellmanFord(self, ranking_node):

		# Step 1: Initialize distances from ranking_node to all other vertices
		# as INFINITE

		# Ignore the terminals by ranking it as a big number
		if(ranking_node in self.terminals):
			# print('continue')
			return 100000

		#  Initialize distances from ranking_node to all other vertices
		dist = [1000000] * self.V
		dist[ranking_node] = 0	# dist to itself is 0
		node_set_to_ranking_node = [[]] * self.V	# e.g. node_set_to_ranking_node[3]=[4,5] means Node set from ranking_node to node 3 is [4,5] 

		not_connected_sources=self.sources		# restore the sources not connected and if sources all is connected, it will be empty
		not_connected_destinations=self.destinations

		isTerminalConnect=[False]*self.V		# restore the connectivity of terminal, is all terminal is connected
												# isTerminalConnect[terminal]=True for all sources and destinations

		# Step 2: Calculate the length 
		while (size(not_connected_destinations)>0 or size(not_connected_sources)>0) and  size(node_set_to_ranking_node)<self.V-size(self.terminals):
			#if there are terminals not to be connected yet or all the nodes are not pooled then we continue to pool

			for u, v, w in self.graph:

				if dist[u] != 1000000 and dist[u] + 1 + w < dist[v] and (not v in node_set_to_ranking_node[u]):
					# not v in node_set_to_ranking_node[u] means we don't pool a node more than once

					# let nodeSet[v] be the node collection from ranking_node to Node v (not include Node v)

					node_set_to_ranking_node[v].extend(node_set_to_ranking_node[u])		# add the node set from ranking_node to u to the one to v
					node_set_to_ranking_node[v].append(u)								# add Node u to the node set from ranking_node to v
					node_set_to_ranking_node[v]=list(set(node_set_to_ranking_node[v]))	# Takes the unique elements in case of duplicated nodes

					dist[v] = dist[u] + 1 + w	# If edge has weight w, we add it to dist v


					# Check for the connectivity about sources and destinations
					for source in not_connected_sources:
						if [source,v,0] in self.graph:
							isTerminalConnect[source]=True
							if isTerminalConnect[self.terminals[source]]==True:
								not_connected_sources.remove(source)	# if a source is connected with some nodes, remove it from not_connected_sources
							break

					for destination in not_connected_destinations:
						if [v,destination,0] in self.graph:
							isTerminalConnect[destination]=True
							if isTerminalConnect[self.terminals[destination]]==True:
								not_connected_destinations.remove(destination)
							break
		
		# Calculate the score by our method: add all dists to terminals together
		score=0
		for source in self.sources:
			score+=dist[source]
		for destination in self.destinations:
			score+=dist[destination]
		return score


	# Get the lowest-score node to pool next
	def getBestPoolingNode(self):
		node_scores=[1000000]*self.V
		for node in range(self.V):
			node_scores[node]=self.BellmanFord(node)

		best_pooling_node=node_scores.index(min(node_scores)) # Get the number of the lowest-score node

		self.node_cost[best_pooling_node]=0 # Apply weight-to-0 to let its cost as 0

		return best_pooling_node
	
	# Get the result by our algorithm
	def getBestPoolingResult(self):
		not_connected_sources=self.sources	# record the sources not yet connected with nodes
		not_connected_destinations=self.destinations	# record the destinations not yet connected with nodes
		isTerminalConnect=[False]*self.V	# record the sources or destinations' connectivity

		iteration=0
		best_pooling_node=-1
		while(iteration <=  self.V-size(self.terminals) and ( size(not_connected_destinations)>0 and size(not_connected_sources)>0 )):
			iteration+=1
			best_pooling_node=self.getBestPoolingNode()


		# Check for the connectivity about sources and destinations
		for source in not_connected_sources:
			if [source,best_pooling_node,0] in self.graph:
				isTerminalConnect[source]=True
				if isTerminalConnect[self.terminals[source]]==True:
					not_connected_sources.remove(source)	# if a source is connected with some nodes, remove it from not_connected_sources
				break

		for destination in not_connected_destinations:
			if [best_pooling_node,destination,0] in self.graph:
				isTerminalConnect[destination]=True
				if isTerminalConnect[self.terminals[destination]]==True:
					not_connected_destinations.remove(destination)
				break
		
		# iterate the Bellman-Ford and get the final result
		best_pooling_result=[]
		for node in range(self.V):

			if(self.node_cost[node]==0 and not node in self.terminals):
				best_pooling_result.append(node)
		return best_pooling_result


		
g=Graph(39,[0,2,4],[1,3,5])
g.addEdges(7,15)
g.addEdges(16,26)
g.addEdges(27,38)
g.addEdgeWithWeight(6,16,100)
g.addEdgeWithWeight(15,26,100)
g.addEdgeWithWeight(16,27,100)
g.addEdgeWithWeight(26,38,100)

g.addInArrows(0,[6,16,27])
g.addOutArrows([15,26,38],1)

g.addInArrows(2,[18,29])
g.addOutArrows([26,38],3)

g.addInArrows(4,[31])
g.addOutArrows([38],5)
print('The graph is ',g.graph)
print('The best pooling result is ',g.getBestPoolingResult())