## @package prottools
# @TODO: needs full header

# @author Alvin Fagan

from sets import Set
# Not available on dell blades
#from visual import *
import numpy
import math

class Node( vector ):
  def __init__( self, idlabel = '', x = 0.0, y = 0.0, z = 0.0, neighbors = None ):
    self.idlabel = idlabel
    # Might need this later: speed, floating point accuracy
    # Use a tuple to prevent changing values
    self.original = ( x, y, z )
    # Use vpython's vector class
    # Although a vector is an ordered pair of points I'm using it this
    # way due to time constraints
    # vpython not available
    self.position = array( x, y, z )
    # Tries to prevent from pointing to atoms not in the atoms
    # dictionary of the chain but this might be too strict
    self.neighbors = Set() if neighbors is None else neighbors
  
  def rotate( self, theta, vector ):
    m_rot = 
  
  # Any connecting should be encapsulated within Chain class for safety

if __name__ == '__main__':
  # Trying to figure out how to rotate vector into xy plane
  vec = []
  vec.append( vector( 1, 0, 0 ) )
  vec.append( vector( -2, 2, 3 ) )
  vec.append( vec[ 0 ].cross( vec[ 1 ] ) )
  vec.append( vec[ 0 ].cross( vec[ 2 ] ) )
  vec.append( degrees( vec[ 0 ].diff_angle( vec[ 2 ] ) ) ) #5
  # Vectors 0, 2 and 3 are 90 degs to each other
  vec.append( degrees( vec[ 3 ].diff_angle( vec[ 2 ] ) ) )
  vec.append( degrees( vec[ 0 ].diff_angle( vec[ 1 ] ) ) )
  #vec.append( vec[ 1 ].rotate( vec[ 
  
  for v in vec:
    print v

if False:
  print __package__
  import string, random

  def connect( node_a, node_b ):
    node_a.neighbors.append( node_b )
    node_b.neighbors.append( node_a )

  my_nodes = []
  chars = string.ascii_uppercase + string.digits
  char_set = random.sample( chars, len( chars ) )
  print char_set, '\n'
  
  limit = 100
  limit_a = -limit
  limit_b = limit
  for i in range( 0, 5 ):
    x = random.uniform( limit_a, limit_b )
    y = random.uniform( limit_a, limit_b )
    z = random.uniform( limit_a, limit_b )
    label = ''.join( char_set[:4] )
    char_set = char_set[4:]
    my_nodes.append( Node( label, x, y, z ) )
        
    print i, '=', my_nodes[ i ].idlabel
    print my_nodes[ i ].x
    print my_nodes[ i ].y
    print my_nodes[ i ].z, '\n'
    
  print my_nodes[ 1 ].neighbors
    
  connect( my_nodes[ 1 ], my_nodes[ 0 ] )
  connect( my_nodes[ 1 ], my_nodes[ 2 ] )
    
  for node in my_nodes[ 1 ].neighbors:
    print node.idlabel
