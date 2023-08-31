load("spermutahedron.py")



def s_oruga_graph(s):
  n = len(s)
  source_edges = []
  if any(s[i]==0 for i in range(n)):
    raise Exception("s cannot be zero")
  for i,x in enumerate(s):
    for j in range(1,x):
      source_edges.append((0,n-i,str(0)+str(n-i)+str(j)))
  bumps_dips = [(n+1-i-1,n+1-i,str(n+1-i-1)+str(n+1-i)+"0") for i in range(1,n+1)]+[(n+1-i-1,n+1-i,str(n+1-i-1)+str(n+1-i)+str(s[i-1])) for i in range(1,n+1)]
  edges = bumps_dips + source_edges
  return DiGraph(edges,multiedges=True)

# Divides a route (path from 1 to n in sections of 2 (i.e. edges) )
def parser_in_2s(route):
  ans = []
  for i in range(len(route)/2):
    ans += [route[2*i]+route[2*i+1]]
  return ans

# Input s: an array of positive integers.
# Input route: a list of strings with possible values "x"+str(i), str(i)+str(j), or "y"+str(i).
# Output: a string following the antilexicografical order with x < any number < y.
def s_next_edge(s,label):
  n = len(s)
  if label[0] == "y":
    raise Exception("There is no next edge for an edge starting with y")
  if label[0] == "2" and int(label[1]) >= s[n-1]+1:
    raise Exception("There is no next edge for the final edge of a source or" +label+"doesn't exist")
  if (label[0] == "x" and int(label[1]) > n) or (label[0] not in ["x","y","2"] and int(label[1]) > s[n-(int(label[0])-1)]-1):
    raise Exception("The edge input "+label+" does not exist")
  if label[0] == "2":
    return "2"+ str(int(label[1])+1)
  if label[0] == "x":
    if s[n-int(label[1])] == 1:
      return "y"+label[1]
    else: # next is a source edge 
      return str(int(label[1])+1)+"1"
  if int(label[1]) != s[n-(int(label[0])-1)]-1: # next is a source edge
    return label[0]+str(int(label[1])+1)
  else: # next is the final edge
    return "y"+str(int(label[0])-1)

# Input s: an array of positive integers.
# Input route: a list of strings with possible values "x"+str(i), str(i)+str(j), or "y"+str(i).
# Input w:  is the element of the permutation by which we are modifying the route
def s_next_route(s,w,route,leftover):
  prefix = route[:-w]
  ix = [s_next_edge(s,route[-w])] # calculates new edge
  if w > 1:
    suffix = route[(-w)+1:]
  else:
    suffix = []
  tentative = prefix+ix+suffix
  num_sources = sum( edge[0] not in ["x","y"] for edge in tentative)
  while( num_sources > 1 or (tentative[0][0] in ["x","y"] and num_sources == 1)):
    leftover += [tentative[0]]
    tentative = tentative[1:]
    if num_sources == 0:
      proof = False
      while(not proof):
        if leftover[-1][0] not in ["x","y"]:
          proof = True
        tentative = [leftover[-1]] + tentative
        leftover = leftover[:-1]
      return(tentative,leftover)
    num_sources = sum( edge[0] not in ["x","y"] for edge in tentative)
  if num_sources == 0:
    proof = False
    while(not proof):
      if leftover[-1][0] not in ["x","y"]:
        proof = True
      tentative = [leftover[-1]] + tentative
      leftover = leftover[:-1]
    return(tentative,leftover)
  return (tentative,leftover)

def s_permutation_to_clique(s,perm):
  n = len(s)
  initial = ["21"] + ["x" + str(i) for i in range(2,n+1)]
  current = [initial,[]]
  clique = [current[0]]
  for w in perm:
    # print (w,current[0],current[1])
    current = s_next_route(s,w,current[0],current[1])
    clique += [current[0]]
  return [''.join(route) for route in clique]

def s_height_delta(s,edges,i,k):
  n=len(s)
  if i+k < len(edges):
    if edges[i][0] == 'x':
      a = 0
    elif edges[i][0] == 'y':
      a = s[n-int(edges[i][1])]
    else:
      a = int(edges[i][1])
    if edges[i+k][0] == 'x':
      b = 0
    elif edges[i+k][0] == 'y':
      b = 1
    else:
      Exception("With our graphs this should not happen. This edge should not be a source edge. i.e. with label starting 'x' or 'y'.")
    return (a+b)**2
  else:
    return 0

def s_route_to_height_k(s,route,k):
  edges = parser_in_2s(route)
  return sum(s_height_delta(s,edges,i,k) for i in range(len(edges)))

def s_route_to_height(s,route,epsilon):
  return sum(s_route_to_height_k(s,route,k)*(epsilon**k) for k in range(1,len(s)+1))

def s_permutation_to_coordinate_i(s,perm,clique,epsilon,i):
  return sum( -s_route_to_height(s,clique[index],epsilon)+s_route_to_height(s,clique[index+1],epsilon) if value == i else 0 for index,value in enumerate(perm))

def s_permutation_to_coordinates(s,perm,epsilon):
  clique = s_permutation_to_clique(s,perm)
  return tuple(s_permutation_to_coordinate_i(s,perm,clique,epsilon,i+1) for i in range(len(s)))

def good_epsilon(s):
  return 2/(sum(2*i+1 for i in s))-0.0000000001

def s_realization_edges(s,epsilon):
  SD = SDecreasingTrees(s)
  L = SD.lattice()
  coordinates = {T:s_permutation_to_coordinates(s,T.to_s_permutation(),epsilon) for T in L}
  if len(s)==4:
    for T in coordinates:
      coordinates[T] = coordinates[T][:-1]
  return [[coordinates[r[0]],coordinates[r[1]]] for r in L.cover_relations()]

def s_realization_points(s,epsilon):
  SD = SDecreasingTrees(s)
  L = SD.lattice()
  return [s_permutation_to_coordinates(s,T.to_s_permutation(),epsilon) for T in L]

def s_realization(s,epsilon):
  edges = s_realization_edges(s,epsilon)
  lines = [line3d(edge,color='blue',thickness=1.5) for edge in edges]
  return sum(lines)