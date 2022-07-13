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
		self.V = nodes # No. of vertices
		self.sources=sources	# nodes for inputing
		self.destinations=destinations	# nodes for outputing
		self.graph = []
		self.in_arrows=[[]]*nodes
		self.out_arrows=[[]]*nodes
		self.terminals=[0]*nodes	# terminals[A]=destination or source of A, A must be a destination of source
		self.node_cost=[1]*nodes	# the costs of nodes

		for i in range(size(sources)):
			self.terminals[sources[i]]=destinations[i]
			self.terminals[destinations[i]]=sources[i]
			self.node_cost[sources[i]]=0
			self.node_cost[destinations[i]]=0

	# function to add an edge to graph
	def addEdgeWithWeight(self, u, v, w):
		self.graph.append([u, v, w])	# the edge from u to v, w stands for the node number between u and v
		self.graph.append([v, u, w])	# the edge from v to u, w stands for the node number between u and v

	# function to add an edge with weight 1 to graph
	def addEdge(self, u, v):
		self.graph.append([u, v, 0])	# the edge from u to v
		self.graph.append([v, u, 0])	# the edge from v to u

	# function to add an edge to isolate two nodes
	def addIsolation(self, u, v):
		self.graph.append([u, v, 1000000])	# the isolation from u to v
		self.graph.append([v, u, 1000000])	# the isolation from v to u

	# function to add edges with weight 1 to graph
	def addEdges(self, u, v):
		for i in range(u,v):	# add the edges from u to v
			self.addEdge(i,i+1)
	
	def addInArrows(self,u,v):
		self.in_arrows[u]=v
		for j in v:	# add the arrows from u to nodes of v
			self.graph.append([u,j,0])

	def addOutArrows(self,u,v):
		self.out_arrows[v]=u
		for i in u:	# add the arrows nodes of u to v
			self.graph.append([i,v,0])

	def BellmanFord(self, src):

		# Step 1: Initialize distances from src to all other vertices
		# as INFINITE

		if(src in self.terminals or src==0):
			# print('continue')
			return 100000

		# print('src is ',src)
		dist = [1000000] * (self.V+1)
		dist[src] = 0
		node_set_to_src = [[]] * (self.V+1)

		not_connected_sources=self.sources
		not_connected_destinations=self.destinations
		isTerminalConnect=[False]*(self.V+1)

		# Step 2: Relax all edges |V| - 1 times. A simple shortest
		# path from src to any other vertex can have at-most |V| - 1
		# edges
		while size(not_connected_destinations)>0 or size(not_connected_sources)>0:
			# Update dist value and parent index of the adjacent vertices of
			# the picked vertex. Consider only those vertices which are still in
			# queue
			if(size(node_set_to_src)==self.V-size(self.terminals)):
				break;

			for u, v, w in self.graph:

				if dist[u] != 1000000 and dist[u] + 1 + w < dist[v] and (not v in node_set_to_src[u]):
					
					# let nodeSet[v] be the node collection from src to Node v (not include Node v)
					print('v is ',v)
					print('u is ',u)
					node_set_to_src[v].extend(node_set_to_src[u])
					node_set_to_src[v].append(u)
					node_set_to_src[v]=list(set(node_set_to_src[v]))
					print('new result node_set_to_src[v] is ',node_set_to_src[v])

					dist[v] = dist[u] + 1 + w

					for source in not_connected_sources:
						if [source,v,0] in self.graph:
							isTerminalConnect[source]=True
							if isTerminalConnect[self.terminals[source]]==True:
								not_connected_sources.remove(source)
							break

					for destination in not_connected_destinations:
						if [v,destination,0] in self.graph:
							isTerminalConnect[destination]=True
							if isTerminalConnect[self.terminals[destination]]==True:
								not_connected_destinations.remove(destination)
							break
						
		score=0
		for source in self.sources:
			score+=dist[source]
		for destination in self.destinations:
			score+=dist[destination]
		return score


	def getBestPoolingNode(self):
		node_scores=[1000000]*(self.V+1)
		for node in range(self.V+1):
			node_scores[node]=self.BellmanFord(node)
		print('node_scores is ',node_scores)
		best_pooling_node=node_scores.index(min(node_scores))
		self.node_cost[best_pooling_node]=0
		print('The new node to be pooled is ',best_pooling_node)

		return best_pooling_node
	
	def getBestPoolingResult(self):
		not_connected_sources=self.sources	# record the sources not yet connected with nodes
		not_connected_destinations=self.destinations	# record the destinations not yet connected with nodes
		isTerminalConnect=[False]*self.V	# record the sources or destinations' connectivity

		iteration=0
		while(iteration <= self.V and ( size(not_connected_destinations)>0 and size(not_connected_sources)>0 )):
			iteration+=1
			self.getBestPoolingNode()
		
		# iterate the Bellman-Ford and get the final result
		best_pooling_result=[]
		for node in range(size(self.node_cost)):
			if(self.node_cost(node)==0):
				best_pooling_result.append(node)
		return best_pooling_result

		
g=Graph(39,[1,3,5],[2,4,6])
g.addEdges(7,16)
g.addEdges(17,27)
g.addEdges(28,39)
g.addEdgeWithWeight(7,17,100)
g.addEdgeWithWeight(16,27,100)
g.addEdgeWithWeight(17,28,100)
g.addEdgeWithWeight(27,39,100)

g.addInArrows(1,[7,17,28])
g.addOutArrows([16,27,39],2)

g.addInArrows(3,[19,30])
g.addOutArrows([27,39],4)

g.addInArrows(5,[32])
g.addOutArrows([39],6)

print(g.getBestPoolingResult())