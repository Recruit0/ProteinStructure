## @package prottools
# @TODO: needs full header

# @author Alvin Fagan
# @author W. Cole Davis

from visual import * 

class Chain:
  def __init__( self, name = '' ):
    self.name = name
    self.bonds = {}
    self.atoms = {}
    
  ## Makes tuple of nodes become neighbors
  # Using simpler algo.
  # More optimal algo of O( n^2 / 2 ) is not as readable.
  def connect( self, node_ids ):
    for node_id in node_ids:
      # Make all nodes in node_ids add all nodes in node_ids
      for new_id in node_ids:
        # self.atoms[ node_id ].neighbors.add( new_id )
        self.atoms[ node_id ].neighbors.add( new_id )
    # Make sure nodes don't have themselves as a neighbor
    self.atoms[ node_id ].neighbors.discard( node_id )
    # self.atoms[ node_id ].neighbors.remove[ node_id ]
  
  def add_bond( self, bond ):
    node_ids = ( bond.nodes[ 0 ].idlabel, bond.nodes[ 1 ].idlabel )
    # Fix bond so same atom won't be added more than once
    # Bonds should point to atoms in self.atoms if they already exist
    for i in range( 0, len( node_ids ) ):
      if node_ids[ i ] in self.atoms:
        bond.nodes[ i ] = self.atoms[ node_ids[ i ] ]
    
    # Inserting same bonds/atoms shouldn't matter at this point
    self.bonds[ node_ids ] = bond
    # Insert nodes from bond into self.atoms
    self.atoms[ node_ids[ 0 ] ] = bond.nodes[ 0 ]
    self.atoms[ node_ids[ 1 ] ] = bond.nodes[ 1 ]
    # Connect nodes in the bond to each other
    self.connect( node_ids )
  
  # Rotates whole chain in space about a 3d axis
  def rotate( self, theta, vector ):
    for atom in self.atoms:
      atom.rotate( theta, vector )
  
  # Translates whole chain in space
  def translate( self, vector ):
    for atom in self.atoms:
      atom += vector
  
  ## Returns a list of all the atoms the chains have in common
  # MAYBE this should be a separate function on its own?
  # Order of tuples is self first then whatever order others were in
  # @param chains List of chains that this one will be compared to
  # @param epsilon Max error distance between points to be considered a match
  def intersect( chains, epsilon = 1.0 ):
    matches = []
    for atom_a in self.atoms:
      all_chains_match = True
      # List of atom ids that matched in the other chains
      matching_ids = []
      
      for chain in chains:
        # Assume there is no match until one is found
        match_exists = False
        
        # Look through all atoms
        for atom_b in chain.atoms:
          # If there is a match then note it and stop looking
          # Use magnitude squared to avoid sqrt computations
          if mag2( atom_b - atom_a ) <= epsilon**2:
            match_exists = True
            matching_ids.append( atom_b.idlabel )
            break
            
        # If there was not a match in one of the chains
        # then ignore atom_a
        if not match_exists:
          all_chains_match = False
          break
      
      # If at least one of the chains didn't match then ignore atom_a
      if not all_chains_match:
        break
      # Otherwise we found a match
      else:
        matches.append( matching_ids )
    
    return matches

if __name__ == '__main__':
  my_nodes = [ 'CD12', 'AB12' ]
  print my_nodes
  
  my_bonds = {}
  key = tuple( sorted( my_nodes ) )
  my_bonds[ key ] = key
  print ( key ) in my_bonds
  print my_bonds[ key ]
