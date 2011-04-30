## @package prottools
# @TODO: needs full header

from operator import itemgetter, attrgetter
from visual import *

# @author Alvin Fagan

class Bond:
  def __init__( self, node_a, node_b ):
    # Store list sorted by id so that one pair has one key
    self.nodes = sorted ( [ node_a, node_b ], key=attrgetter('idlabel') )
    self.magnitude = mag( node_a - node_b )
    # Marker so that this bond doesn't have to be looked for again
    # while searching.
    # MAKE SURE TO RESET TO FALSE IF CHANGED
    self.not_found = False
    # Should not connect at this point, connecting requires having
    # a list of all the atoms, hence connecting should be handled by
    # the Chain class
    # self.node_a.connect( node_b )

if __name__ == '__main__':
  class Point:
    def __init__( self, idlabel = '', x = 0.0, y = 0.0 ):
      self.idlabel = idlabel
      self.x = x
      self.y = y
  
  point1 = Point( 'CD34', 1.0, 2.0 )
  point2 = Point( 'AB12', 3.0, 4.0 )
  
  unsorted_points = [ point1, point2 ]
  for point in unsorted_points:
    print point.idlabel
  points = sorted( unsorted_points, key=attrgetter('idlabel') )
  for point in points:
    print point.idlabel, point.x, point.y
  
  numbers = sorted( [ 3, 5, 2, 6, 1 ] )
  print numbers
  for i in range( 0, len( numbers ) ):
    print numbers[ i ]
