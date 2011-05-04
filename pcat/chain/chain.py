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
# Obselete?
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
    for my_atom in self.atoms:
      my_atom.rotate( theta, vector )
  
  # Translates whole chain in space
  def translate( self, vector ):
    for my_atom in self.atoms:
      my_atom += vector
  
  ## Returns first atom found near specified position, otherwise None
  # @param epsilon Max error distance from position
  # Upgrade: Use BSP or something
  def atom_at( self, position, epsilon = 1.0 ):
    found_atom = None
    for my_atom in self.atoms:
      # Use magnitude squared to avoid square root ops
      if mag2( my_atom - position ) <= epsilon**2:
        found_atom = my_atom
        break
    return found_atom
  
  ## Find atoms that are in the same position in all chains and returns
  ## a list of their IDs. Order is self first then whatever order other
  ## chains were in.
  # @param chains List of chains to intersect with this one
  # @param epsilon Max error distance between points to be considered a match
  def intersect( self, chains, epsilon = 1.0 ):
    matches = []
    for my_atom in self.atoms:
      all_chains_match = True
      # List of atom ids that matched in the other chains
      # Make sure to add my_atom so it is in the results too
      matching_ids = [ my_atom.idlabel ]
      # Construct list of ids of matching atoms in other chains
      for chain in chains:
        # Check if chain has atom in same position
        matching_atom = chain.atom_at( my_atom, epsilon )
        if matching_atom is not None:
          # Found a match
          matching_ids.append( matching_atom.idlabel )
        else:
          # One of the chains does not match for this atom
          all_chains_match = False
          break
      # If a match was found in all chains then add it to results
      if all_chains_match:
        matches.append( matching_ids )
    return matches
  
  ## Resets not_found flag for all bonds to false
  #
  def reset_not_found( self ):
    for bond in self.bonds:
      bond.not_found = False
  
  ## Returns a list of bonds with specified length/magnitude.
  # @param d_r Max difference in magnitude for match
  def find_length( self, magnitude, d_r = 1.0 ):
    #matches = []
    #for bond in self.bonds:
    #  if abs( bond.magnitude - magnitude ) <= d_r:
    #    matches.append( bond )
    #return matches
    # Use list comprehension instead
    return [bond for bond in self.bonds if abs( bond.magnitude - magnitude ) <= d_r]
  
  ## Returns angle between seed and key or None if not a seed-key pair
  # Both seed and key are bonds
  # @TODO Remove returning None safety for faster execution
  def sk_angle( self, seed, key ):
    # Vertice called the origin because it will be moved there during
    # seed-key alignment
    origin = None
    joint = None
    end = None
    # Find the joint and set origin accordingly
    # Duplicate node ID issue not resolved
    # Design to restrict node access to self only
    for seed_node in seed.nodes:
      for key_node in key.nodes:
        if seed_node.idlabel == key_node.idlabel:
          joint = seed_node
        else:
          origin = seed_node
    
    # Set end node
    for key_node in key.nodes:
      if key_node.idlabel != joint.idlabel:
        # Self access restriction only not implemented yet
        end = key_node
    
    # All vertices must be defined
    if origin is not None and joint is not None and end is not None:
      # Vector from joint to origin
      v1 = origin - joint
      # Vector from joint to end
      v2 = end - joint
      # Inner angle between v1 and v2
      return v1.diff_angle( v2 )
    # Return None if seed-key pair does not connect
    return None
  
  ## Generates list of neighboring bonds for specified bond
  # @param Bond to get neighbors for
  # @TODO Use self reference only constraint
  def get_neighbors( self, bond ):
    neighbors = []
    for node in bond.nodes:
      # This part uses self reference only constraint
      for neighbor_id in node.neighbors:
        bond_key = sorted ( [ node.idlabel, neighbor_id ] )
        neighbors.append( self.bonds[ bond_key ] )
    return neighbors
  
  ## Finds neighbors of bond with specified length
  # @param d_r Max difference in length
  def find_neighbors( self, bond, length, d_r ):
    neighbors = []
    for node in bond.nodes:
      # This part uses self reference only constraint
      for neighbor_id in node.neighbors:
        bond_key = sorted ( [ node.idlabel, neighbor_id ] )
        neighboring_bond = self.bonds[ bond_key ]
        if abs( neighboring_bond - length ) <= d_r:
          neighbors.append( self.bonds[ bond_key ] )
    return neighbors
  
  ## Returns list of matching SKPs in this chain
  # @param d_r Max difference in magnitude
  # @param d_theta Max difference in angle
  def match_skp( self, seed, key, d_r = 1.0, d_theta = 1.0 ):
    matching_skps = []
    # Function sk_angle does not have self access only safety
    theta_1 = self.sk_angle( seed, key )
    matching_seeds = self.find_length( seed.magnitude, d_r )
    for potential_seed in matching_seeds:
      matching_keys = self.find_neighbors( potential_seed, key.magnitude, d_r )
      for potential_key in matching_keys:
        # Lengths match, now check angle
        theta_2 = self.sk_angle( potential_seed, potential_key )
        if abs( theta_1 - theta_2 ) <= d_theta:
          # Found a matching SKP
          new_skp = [ potential_seed, potential_key ]
          matching_skps.append( new_skp )
    return matching_skps
  
  ## Aligns chain so that SKP is in XY plane with seed on x-axis
  ## No result if seed and key are not connected
  #
  def align( self, seed, key ):
    # Vertice called the origin because it will be moved there during
    # seed-key alignment
    origin = None
    joint = None
    end = None
    # Find the joint and set origin accordingly
    # Duplicate node ID issue not resolved
    # Design to restrict node access to self only
    for seed_node in seed.nodes:
      for key_node in key.nodes:
        if seed_node.idlabel == key_node.idlabel:
          joint = seed_node
        else:
          origin = seed_node
    
    # Set end node
    for key_node in key.nodes:
      if key_node.idlabel != joint.idlabel:
        # Self access restriction only not implemented yet
        end = key_node
    
    # All vertices must be defined
    if origin is not None and joint is not None and end is not None:
      # Move chain so that vertice "origin" is the origin
      self.translate( -1.0 * origin )
      # Vector from origin to joint
      # Subtracting the origin is mostly for illustration purposes
      v1 = joint - origin
      # Vector from joint to end
      v2 = end - joint
      # Orthographic projection into ZY plane, i.e. ignore x
      # Then take angle in this plane
      x_angle = atan2( v1.z, v1.y )
      # Rotate around x-axis so seed lines up with y axis when looking
      # down at x-axis.
      # Seed will be in XY plane
      self.rotate( -x_angle, vector( 1, 0, 0 ) )
      # Rotate around z-axis so seed lines up with x axis
      z_angle = atan2( v1.y, v1.x )
      self.rotate( -z_angle, vector( 0, 0, 1 ) )
      # Ortho project key into ZY plane, i.e. ignore x
      # Look down at X axis, this is ZY plane
      # Rotate so that key lines up with Y axis
      x2_angle = atan2( v2.z, v2.y )
      # Key will be in XY plane
      self.rotate( -x2_angle, vector( 1, 0, 0 ) )
  
  ## Returns a list of seeds and keys matching in all chains
  # @param d_r Max difference in magnitude
  # @param d_theta Max difference in angle, in radians. There is a
  # degrees to radians helper function from vpython.
  def seed_key( self, chains, d_r = 1.0, d_theta = 1.0 ):
    # List of list. Each sub-list has self seed and key pair first then
    # seed-key pairs (SKPs) in same order as chains.
    skp_matches = []
    for my_seed in self.bonds:
      for my_key in self.get_neighbors( my_seed ):
        my_skp = [ my_seed, my_key ]
        potential_skp_match = [ my_skp ]
        skp_missing = False
        for other_chain in chains:
          matching_skps = other_chain.match_skp( my_seed, my_key, d_r, d_theta )
          if len( matching_skps ) > 0:
            # Found match(es)
            # Using only first one
            # Couldn't figure out how to use other ones
            potential_skp_match.append( matching_skps[ 0 ] )
          else:
            # No matching SKPs
            skp_missing = True
            break
        if not skp_missing:
          skp_matches.append( potential_skp_match )
    return skp_matches

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
