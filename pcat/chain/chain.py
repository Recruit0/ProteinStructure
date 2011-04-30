## @package prottools
# @TODO: needs full header

# @author Alvin Fagan
# @author W. Cole Davis

from visual import * 
import math
from operator import itemgetter, attrgetter

## Returns if a list of chains has a magnitude with max d_r error
# Obselete?
def has_mag( chains, mag, d_r = 1.0 ):
  all_match = True
  for chain in chains:
    match_found = False
    # Optimize this later to do binary search or something
    for bond in chain.bonds:
      if abs( bond.magnitude - mag ) <= d_r:
        # Found a match
        match_found = True
        break
    if not match_found:
      # No match in one of the chains means no match in all
      all_match = False
      break
  return all_match

## Returns if list of chains have a specified atom with max epsilon
## error in distance.
# @param epsilon Max difference in position
# @param position This must be a vpython vector
def has_atom( chains, position, epsilon = 1.0 ):
  all_match = True
  for chain in chains:
    match_found = False
    # Optimize this later to do binary search or something
    for atom in chain.atoms:
      # Figure out how to make this a parameter
      if mag2( atom - position ) <= epsilon**2:
        # Found a match
        match_found = True
        break
    if not match_found:
      # No match in one of the chains means no match in all
      all_match = False
      break
  return all_match

## Returns if list of chains have a specified atom with max epsilon
## error in distance.
# @param epsilon Max difference in position
# @param position This must be a vpython vector
def has_atom( chains, position, epsilon = 1.0 ):
  all_match = True
  for chain in chains:
    match_found = False
    # Optimize this later to do binary search or something
    for atom in chain.atoms:
      # Figure out how to make this a parameter
      if mag2( atom - position ) <= epsilon**2:
        # Found a match
        match_found = True
        break
    if not match_found:
      # No match in one of the chains means no match in all
      all_match = False
      break
  return all_match

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
  
  ## Adds bond to chain
  # This needs to be fixed, duplicate ID different chain issue
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
  def intersect( self, chains, epsilon = 1.0 ):
    matches = []
    for atom in self.atoms:
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
          # Magnitude in bond dictionary is only for same chain
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
  
  ## Returns list of bonds with specified magnitude
  # @param d_r Max difference in magnitude for lookup
  def find_length( self, magnitude, d_r = 1.0 ):
    #matches = []
    #for bond in self.bonds:
    #  if abs( bond.magnitude - magnitude ) <= d_r:
    #    matches.append( bond )
    #return matches
    # Use list comprehension instead
    return [bond for bond in self.bonds if abs( bond.magnitude - magnitude ) <= d_r]
  
  ## Returns first atom found near specified position, otherwise None
  # @param epsilon Max error distance from position
  def atom_at( self, position, epsilon = 1.0 ):
    found_atom = None
    for atom in self.atoms:
      # Use magnitude squared to avoid square root ops
      if mag2( atom - position ) <= epsilon**2:
        found_atom = atom
        break
    return found_atom
  
  ## Returns angle between seed and key or None if pair doesn't connect
  ## Both seed and key are bonds
  def sk_angle( self, seed, key ):
    angle = None
    origin = None
    joint = None
    end = None
    for seed_node in seed.nodes:
      if node.idlabel == key[ 0 ].idlabel:
        joint = node
        end = key[ 1 ]
      # Else other way around
      else:
        origin = node
        end = key[ 0 ]
        
    for node in seed.nodes:
      if node is not joint:
        
    if joint is not None:
      if origin is not None:
        v1 = joint - origin
        v2 = 
  
  ## Returns a list of seeds and keys matching in all chains
  #
  def seed_key( chains, d_r = 1.0, d_theta = 1.0 ):
    matches = []
    for seed in self.bonds:
      if has_mag( chains, seed.mag, d_r ):
        # Found a seed
        for node in seed:
          # The system may break here if connecting chains together
          # due to duplicate IDs in different chains
          for neighbor_id in node.neighbors:
            key_ids = sorted ( [ node.idlabel, neighbor_id ], key=attrgetter('idlabel') )
            key = bonds[ key_ids ]
            if has_mag( chains, key.mag, d_r ):
              pass
              #Found key
              # BOOKMARK
      #else:
      #  bond_a.

            # Found a seed with matching magnitude
            # Look through neighbors, these are connected bonds
            # If we can't find same magnitude in other chains then good
            # This will throw out a lot of mismatches
            # An angle requires 3 nodes/2 bonds, use as pruning function
            # Neighbors does not point outside chain so can always
            # lookup magnitude in bond dictionary of chain
            

if __name__ == '__main__':
  my_nodes = [ 'CD12', 'AB12' ]
  print my_nodes
  
  my_bonds = {}
  key = tuple( sorted( my_nodes ) )
  my_bonds[ key ] = key
  print ( key ) in my_bonds
  print my_bonds[ key ]
  # Using list comprehension
  print [x for x in range(0,10) if x > 5]
