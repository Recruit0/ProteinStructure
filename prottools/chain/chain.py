## @package prottools
# @TODO: needs full header

from sets import Set

class Chain:
  def __init__( self, name = '' ):
    self.name = name
    self.bonds = {}
    self.atoms = {}
    
  ## Makes tuple of nodes become neighbors
  # Each id in node_ids must exist in self.atoms
  # Not sure this should be public function or just stick in add_bond()
  # Using simpler algo. Potentially uses more memory.
  # More optimal algo of O( n^2 / 2 ) is not as readable.
  def connect( self, node_ids ):
    for node_id in node_ids:
      # Make all nodes in node_ids add all nodes in node_ids
      for new_id in node_ids:
        self.atoms[ node_id ].neighbors.add( new_id )
      # Make sure nodes don't have themselves as a neighbor
      atoms[ node_id ].neighbors.remove[ node_id ]
  
  def add_bond( node_a, node_b ):
    new_bond = Bond( node_a.idlabel, node_b.idlabel )
    # BOOKMARK 24.4.2011
    # Inserting same bonds/atoms shouldn't matter at this point
    self.bonds[ new_bond.node_ids ] = new_bond
    # Insert nodes from bond into self.atoms
    self.atoms[ node_ids[ 0 ] ] = node_a
    self.atoms[ node_ids[ 1 ] ] = bond.nodes[ 1 ]
    # Connect nodes in the bond to each other
    connect( node_ids )

if __name__ == '__main__':
  my_nodes = [ 'CD12', 'AB12' ]
  print my_nodes
  
  my_bonds = {}
  key = tuple( sorted( my_nodes ) )
  my_bonds[ key ] = key
  print ( key ) in my_bonds
  print my_bonds[ key ]
