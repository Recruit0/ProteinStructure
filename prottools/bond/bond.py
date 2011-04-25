## @package prottools
# @TODO: needs full header

class Bond:
  # This should not point to nodes directly so that it is impossible for
  # it to point to nodes not the chain, only the actual chain knows
  # what nodes are in it
  def __init__( self, node_a, node_b ):
    # Store tuple sorted by id so that one pair has one key
    node_list = sorted ( [ node_a, node_b ] )
    self.node_ids = tuple( node_list[ 0 ].idlabel, node_list[ 1 ].idlabel )
    # Define vector as directed from node with min( idlabel) to node
    # with max( idlabel)
    self.x = node_list[ 1 ].x - node_list[ 0 ].x
    self.y = node_list[ 1 ].y - node_list[ 0 ].y
    self.z = node_list[ 1 ].z - node_list[ 0 ].z
    # BOOKMARK 24.4.2011

    # Connecting requires having a list of all the atoms, hence
    # connecting should be handled by the Chain class

if __name__ == '__main__':
  class Point:
    def __init__( self, idlabel = '', x = 0.0, y = 0.0 ):
      self.idlabel = idlabel
      self.x = x
      self.y = y
  
  point1 = Point( 'AB12', 1.0, 2.0 )
  point2 = Point( 'CD34', 3.0, 4.0 )
  unsorted_points = [ point2, point1 ]
  for point in unsorted_points:
    print point.idlabel
  points = sorted( unsorted_points )
  for point in points:
    print point.idlabel, point.x, point.y
  
  numbers = sorted( [ 3, 5, 2, 6, 1 ] )
  print numbers
  for i in range( 0, len( numbers ) ):
    print numbers[ i ]
