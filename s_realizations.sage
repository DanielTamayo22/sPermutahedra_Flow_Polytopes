load("spermutahedron.py")

def s_oruga_graph(s):
  """
  Generates the s-oruga graph.

  Parameters
  ----------
  s : list of int
      A composition represented as a list of positive integers.

  Returns
  -------
  networkx.DiGraph
      The s-oruga graph represented as a directed graph.

  Raises
  ------
  Exception
      If any element in the composition 's' is zero.

  """
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

def parser_in_2s(route):
  """
  Divide a route of a graph into sections of 2, representing edges in the path from 1 to n.

  Parameters
  ----------
  route : list
      List representing a route, where consecutive elements form edges in the path from 1 to n.

  Returns
  -------
  list
      A list of 2-tuples, each representing an edge in the path.

  Examples
  --------
  >>> parser_in_2s(['x1', 'x2', 'x3'])
  [('x1', 'x2'), ('x2', 'x3')]

  """
  ans = []
  for i in range(len(route)/2):
    ans += [route[2*i]+route[2*i+1]]
  return ans


def s_next_edge(s,label):
  """
  Determine the next edge in "antilexicographical" order based on x<str(i)<y with i being an integer.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  label : str
      A string representing the current edge label. Possible values:
      - "x" + str(i): Represents an edge labeled x followed by an integer i.
      - str(i) + str(j): Represents an edge labeled with two integers i and j.
      - "y" + str(i): Represents an edge labeled y followed by an integer i.

  Returns
  -------
  str
      The label of the next edge in the induced order.

  Raises
  ------
  Exception
      If the provided label is invalid or there is no next edge.

  Examples
  --------
  >>> s_next_edge([2, 3, 1], "x2")
  '21'
  >>> s_next_edge([2, 3, 1], "21")
  '22'
  >>> s_next_edge([2, 3, 1], "22")
  'y2'
  >>> s_next_edge([2, 3, 1], "y2")
  Exception: There is no next edge for an edge starting with y

  """
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

  def s_next_route(s,w,route,leftover):
    """
    Determine the next route based on the element 'w' of a Stirling s-permutation.

    Parameters
    ----------
    s : list of int
        An array of positive integers.
    w : int
        The element of the Stirling s-permutation used to modify the route.
    route : list of str
        A list of strings representing the current route, with possible values
        "x"+str(i), str(i)+str(j), or "y"+str(i).
    leftover : list of str
        A list representing leftover edges from previous steps.

    Returns
    -------
    tuple
        A tuple containing the next route and the updated leftover edges.

    Notes
    -----
    The function modifies the given route based on the element of the Stirling s-permutation 'w'.
    It calculates a new edge using 's_next_edge' and updates the route accordingly.
    The modified route is then checked for validity, and if necessary, leftover edges
    are moved to maintain the validity of the route.

    Examples
    --------
    >>> s_next_route([2, 3, 1], 2, ['x1', 'x2', 'x3'], [])
    (['21', 'x3'], ['x1'])
    >>> s_next_route([2, 3, 1], 1, ['21', 'x3'], ['x1'])
    (['31'], ['x1','21'])
    >>> s_next_route([2, 3, 1], 1, ['31'], ['x1','21'])
    (['21','y3'], ['x1'])
    >>> s_next_route([2, 3, 1], 2, ['21','y3'], ['x1'])
    (['22','y3'], ['x1'])
    >>> s_next_route([2, 3, 1], 2, ['22','y3'], ['x1'])
    (['x1','y2','y3'])
    >>> s_next_route([2, 3, 1], 3, ['x1','y2','y3'])
    (['y1','y2','y3'])

    """
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
  """
  Convert a Stirling s-permutation to a clique in the s-oruga graph.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  perm : list of int
      A Stirling s-permutation used to generate the clique in the s-oruga graph.

  Returns
  -------
  list of str
      A list of strings representing the clique in the s-oruga graph.

  Notes
  -----
  The function takes a Stirling s-permutation 'perm' and generates a clique in the s-oruga graph. It uses 's_next_route' to modify the route for each element in the Stirling s-permutation and collects the resulting routes to
  form the clique.

  Examples
  --------
  >>> s_permutation_to_clique([2, 3, 1], [2,1,1,2,2,3])
  ['21', 'x2', 'y2', 'x1', '31']
  [['x1', 'x2', 'x3'],['21', 'x3'],['31'],['21','y3'],['22','y3'],['x1','y2','y3'],['y1','y2','y3']]

  """
  n = len(s)
  initial = ["21"] + ["x" + str(i) for i in range(2,n+1)]
  current = [initial,[]]
  clique = [current[0]]
  for w in perm:
    current = s_next_route(s,w,current[0],current[1])
    clique += [current[0]]
  return [''.join(route) for route in clique]

def s_height_delta(s,edges,i,k):
  """
  Calculate the delta for the DKK height function for a pair of edges in the s-oruga graph.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  edges : list of str
      A pair of strings representing two edges included in a route of the s-oruga graph.
  i : int
      Index of the first edge in the pair.
  k : int
      The difference in indices between the two edges in the pair.

  Returns
  -------
  int
      The height delta for the given pair of edges following DKK.

  Raises
  ------
  Exception
      If the the second edge is a source edge.

  """
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
  """
  Calculate the delta sum for the DKK height function for a route in the s-oruga graph.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  route : list
      List representing a route, where consecutive elements form edges in the path from 1 to n.
  k : int
      The difference in indices between the pairs of edges for DKK height delta calculations.

  Returns
  -------
  int
      The sum of height deltas for the given route and edge index difference 'k'.

  """
  edges = parser_in_2s(route)
  return sum(s_height_delta(s,edges,i,k) for i in range(len(edges)))

def s_route_to_height(s,route,epsilon):
  """
  Returns the DKK height of a route in the s-oruga graph depending on epsilon.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  route : list
      List representing a route in the s-oruga graph.
  epsilon : float
      A weight factor applied to each height delta term.

  Returns
  -------
  float
      The DKK height for the given route with the specified epsilon.

  """
  return sum(s_route_to_height_k(s,route,k)*(epsilon**k) for k in range(1,len(s)+1))

def s_permutation_to_coordinate_i(s,perm,clique,epsilon,i):
  """
  Calculates the i-th coordinate corresponding to the Stirling s-permutation 'perm' using DKK heigh functions.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  perm : list of int
      A Stirling s-permutation used to generate the clique in the s-oruga graph.
  clique : list of str
      A list of strings representing the clique in the s-oruga graph.
  epsilon : float
      A weight factor applied to the DKK height function.
  i : int
      The index of the coordinate value that is calculated.

  Returns
  -------
  float
      The coordinate value for the specified index 'i'.

  """
  return sum( -s_route_to_height(s,clique[index],epsilon)+s_route_to_height(s,clique[index+1],epsilon) if value == i else 0 for index,value in enumerate(perm))

def s_permutation_to_coordinates(s,perm,epsilon):
  """
  Returns the coordinate corresponding to a Stirling s-permutation following the s-oruga graph.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  perm : list of int
      A Stirling s-permutation used to generate the clique in the s-oruga graph.
  epsilon : float
      A weight factor applied to each term of the DKK height function.

  Returns
  -------
  tuple of float
      The coordinate for the given Stirling s-permutation.

  """
  clique = s_permutation_to_clique(s,perm)
  return tuple(s_permutation_to_coordinate_i(s,perm,clique,epsilon,i+1) for i in range(len(s)))

def good_epsilon(s):
  """
  Calculates a suitable value for epsilon based on the composition 's'.

  Parameters
  ----------
  s : list of int
      An array of positive integers.

  Returns
  -------
  float
      A calculated value for epsilon based on the given array 's'.

  Examples
  --------
  >>> good_epsilon([2, 3, 1])
  0.133333333233333

  """
  return 2/(sum(2*i+1 for i in s))-0.0000000001

def s_realization_edges(s,epsilon):
  """
  Generates the edges for our realization of the s-permutahedron.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  epsilon : float
      A weight factor applied to each height delta term.

  Returns
  -------
  list of list
      A list of pairs of coordinates representing the edges of our realization of the s-permutahedron.

  Notes
  -----
  The function generates edges for the s-permutahedron based on the given array 's' and weight factor 'epsilon'. It uses the SDecreasingTrees class to obtain the lattice and calculates coordinates for each element in the lattice using 's_permutation_to_coordinates'.

  """
  SD = SDecreasingTrees(s)
  L = SD.lattice()
  coordinates = {T:s_permutation_to_coordinates(s,T.to_s_permutation(),epsilon) for T in L}
  if len(s)==4:
    for T in coordinates:
      coordinates[T] = coordinates[T][:-1]
  return [[coordinates[r[0]],coordinates[r[1]]] for r in L.cover_relations()]

def s_realization_points(s,epsilon):
  """
  Generates the vertices for our realization of the s-permutahedron.

  Parameters
  ----------
  s : list of int
      An array of positive integers.
  epsilon : float
      A weight factor applied to each height delta term.

  Returns
  -------
  list of list
      A list coordinates representing the vertices of our realization of the s-permutahedron.

  Notes
  -----
  The function generates vertices for the s-permutahedron based on the given array 's' and weight factor 'epsilon'. It uses the SDecreasingTrees class to obtain the lattice and calculates coordinates for each element in the lattice using 's_permutation_to_coordinates'.

  """
  SD = SDecreasingTrees(s)
  L = SD.lattice()
  return [s_permutation_to_coordinates(s,T.to_s_permutation(),epsilon) for T in L]

def s_realization(s,epsilon):
  """
    Generate a 3D visualization of the s-permutahedron based using DKK height functions.

    Parameters
    ----------
    s : list of int
        An array of positive integers.
    epsilon : float
        A weight factor for the DKK heigh function.

    Returns
    -------
    Sum of 3D lines
        The visualization of the s-permutahedron represented as the sum of 3D lines.

    """
  edges = s_realization_edges(s,epsilon)
  lines = [line3d(edge,color='blue',thickness=1.5) for edge in edges]
  return sum(lines)