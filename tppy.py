import datetime
from ase.io import read, write
from ase.visualize import view
from dlePy.vasp.chgcar import read_chgcar, write_chgcar
from ase.build import add_adsorbate
import numpy as np

def read_chgcar( INDATA, CONTCAR = "CONTCAR" ):
	print ( "Input data: ", INDATA )
	print ( "Read ", INDATA )
	print ( datetime.datetime.now() )
	rho = VaspChargeDensity( INDATA )
	print ( datetime.datetime.now() )
	return rho


# Function to generate the "ICONST" file for collective coordinates
#atoms_list: A list of lists, where each sublist contains two atom indices that define a collective coordinate.
#choice: An optional parameter. If defined, it changes the format of the last line in the file.
def get_ICONST( atoms_list, choice = None ):
	with open("ICONST", "w") as file:
		for i in atoms_list:
			file.write( "R " + str( i[ 0 ] + 1 ) + " " + str( i[ 1 ] + 1 ) + " 0 \n" )
		if choice:
			file.write( "S 1 -1 0 \n" )
		else:
			file.write( "S 1 0 0 0 \n" )


#Shift the entire system so that the atom at index N is centered in the XY-plane of the simulation cell.
#system (Atoms): The atomic system (structure) read using ASE, e.g., via `ase.io.read`.
#N (int): Index of the atom to center. Typically, this is the index of an atom in the dissociating molecule (e.g., Oxygen in H2O or Nitrogen in a cation).
def shift_center( system, N:int ):
	if N < 0 or N >= len(system):
		raise ValueError(f"Invalid atom index {N}. It must be between 0 and {len(system) - 1}.")
	center = system[N].position
	shift = -center + 0.5 * (system.cell[0] + system.cell[1])
	shift[2] = 0
	for atom in system:
		atom.position += shift
	system.wrap()
	return system

#add velocity for MD simulations
def add_velocity( path_to_MD_CONTCAR ):
	with open( path_to_MD_CONTCAR +  "/CONTCAR", 'r' ) as f:
		lines = f.readlines()

	fout =open( 'POSCAR', 'a' )
	w = False
	for line in lines:
		if line.strip() == '' and not w:
			w = True
		if line.strip() == '1':
			w = False
		if w:
			fout.write( line.strip() + '\n')
	fout.close()
	new_struc = read( "POSCAR" )
	return new_struc

#atoms (Atoms): The complete atomic structure.
#group (list or array): Indices of atoms to be rotated.
#axis (tuple or list): Rotation axis as a vector (e.g., (1, 0, 0)).
#angle (float): Rotation angle (degrees or radians).
#center (tuple or list): Point around which the rotation occurs.
def rotate_group( atoms, group, axis, angle, center ):
    subsys = atoms[ group ]
    subsys.rotate( angle, axis, center, rotate_cell = False )
    atoms.positions[ group ] = subsys.positions
    return atoms

def change_atomic_number( system, lst ):
    for i in lst:
        if system[ i ].symbol == "O" or system[ i ].symbol == "N":
            system[ i ].number = 4
        elif system[ i ].symbol == "H":
            system[ i ].number = 2
    view( system )

def delete_from_list(list_, elements_to_delete):
    elements_to_delete.sort(reverse = True)
    for i in elements_to_delete:
        del l[i]
    print(l)
    return l

def total_charge( file ):
    data = read_chgcar( file  )
    system = data.atoms[ 0 ]
    rho = data.chg[ 0 ]
    max_F = np.max( rho )
    min_F = np.min( rho )
    ne = np.sum( rho ) * system.get_volume() / rho.shape[0] / rho.shape[1] / rho.shape[ 2 ]
    return ne, max_F, min_F


def get_Fmax_Fmin( file ):
	values = list()
	with open( file, "r" ) as read_file:
		for i in read_file:
			if len( i.split() ) == 5:
				values.extend( i.split()  )
	values = [ float(i) for i in values ]
	min_value = min( values )
	max_value = max( values )
	print("F max is: ", max_value )
	print("F min is: ", min_value )
	return max_value, min_value


def swap_atoms(system, atom1, atom2):
	aux = list()
	for i in range(0, 3):
		aux.append(system[atom1].position[i])
		system[atom1].position[i] = system[atom2].position[i]
		system[atom2].position[i] = aux[i]
	return system

def change_position(system, atom, x, y, z):
	l = [x, y, z]
	for i in range(0, len(l)):
		system[atom].position[i] = l[i]
	return system

def get_fermi_energy_VASP(path):
	os.chdir(path)
	word = "E-fermi"
	fermi_energy = ""
	with open("OUTCAR", "r") as file:
		for line_number, line in enumerate(file):
			if word in line:
				a = line
	b = a.replace(" E-fermi :  ", "")			
	for i in b:
		if i != " ":
			fermi_energy = fermi_energy + i
		else:
			break
	return float(fermi_energy)

def get_fermi_energy_VASPsol():
        word = "EFERMI "
        with open("OUTCAR", "r") as file:
                for line_number, line in enumerate(file):
                        if word in line:
                                a = line
        b = a.replace("    EFERMI              = ", "")
        c = b.replace("eV", "")
        c = float(c)
        print(c)
        return c


def get_NELECT():
	word = "UPDATED NELECT"
	with open("OUTCAR", "r") as file:
        	for line_number, line in enumerate(file):
                	if word in line:
                        	a = line

	b = re.sub("UPDATED NELECT      =", "", a)
	c = re.sub("electrons", "", b)
	d = float(c)
	#print("Current NELECT = ", d, " electrons")
	return d

#struc = the structure 
#atom1 = integer number, it is the index of the first atom
#atom2 = integer number, it is the index of the second atom
def get_average_positions(struc, atom1, atom2):
	cordinates = ["x", "y", "z"]
	average_positions = {}
	for i in range(0, 3):
		average_positions[str(cordinates[i])] = (struc[atom1].position[i] + struc[atom2].position[i])/2
	keys = list(average_positions.keys())
	average_cordinates = np.array([average_positions[keys[0]],average_positions[keys[1]] ,average_positions[keys[2]] ])
	return average_cordinates


#This function adds atoms to the system
#system = the structure
#symbol = the type of atom we want to add
#x_pos = x_cordinate, it's just a number
#y_pos = y_cordinate, it's just a number 
#z_pos = z_cordinate, it's just a number 
def add_atoms(system, symbol, x_pos, y_pos, z_pos):
	add_adsorbate(system, symbol, -1, position = (x_pos, y_pos), mol_index = 0 )
	system[len(system) - 1 ].position[2] = z_pos 	
	return system


#This function deletes all atoms, their indexes are in the atoms_list
#struc is the POSCAR or CONTCAR
#atoms_list =  list with the indices of atoms we want to remove
def delete_atoms(struc, atoms_list):
	for i in reversed(struc):
		for j in range(0, len(atoms_list)):
			if i.index == atoms_list[j]:
				del struc[i.index]
	return struc


#This function moves a molecule symmetrically on the surface
#system =  the POSCAR or CONTCAR file
#at = the index of an atom of the molecule we want to move.All other atoms will be moved in repsect to this atom
#layer_atom = the index of the atom above of which the at atom will be moved
#atoms_to_move = a list with the atoms we want to move
def move_symmetrically(system, at, layer_atom, atoms_to_move ):
	x_i = system[at].position[0]
	y_i = system[at].position[1]
	for i in range(0, 2):
		system[at].position[i] = system[layer_atom].position[i]
	moving_atoms = [i for i in atoms_to_move if i != at]
	for i in moving_atoms:
  		system[i].position[0] += system[at].position[0] - x_i
  		system[i].position[1] += system[at].position[1] - y_i
	return system

#This function moves atoms symmetrucally on the x-y. It is the same with move_symmetrically
#but instead of a surface atom it gets the new cordinates in a list.
#system =  the POSCAR or CONTCAR file
#at = the index of an atom of the molecule we want to move. All other atoms will be moved in repsect to this atom
#x_y_cords = the new cordinates of which the at atom will be moved. It is a list. First element is x cordinate and second element is the y cordinate
#atoms_to_move = a list with the atoms we want to move
def move_symmetrically_cords(system, at, x_y_cords, atoms_to_move ):
	x_i = system[at].position[0]
	y_i = system[at].position[1]
	for i in range(0, 2):
		system[at].position[i] = x_y_cords[i]
	moving_atoms = [i for i in atoms_to_move if i != at]
	for i in moving_atoms:
  		system[i].position[0] += system[at].position[0] - x_i
  		system[i].position[1] += system[at].position[1] - y_i
	return system

#This function moves two different systems of atoms symmetrically on the x-y plane.
#MOLECULE1 <--> MOLECULE2 symmetrically 
#system =  the POSCAR or CONTCAR file
#at_1 = the index of an atom of the FIRST molecule we want to move. All other atoms will be moved in repsect to this atom
#at_2 = the index of an atom of the SECOND molecule we want to move. All other atoms will be moved in repsect to this atom
#atoms_to_move_1 = a list with the atoms we want to move that belong to the FIRST molecule
#atoms_to_move_2 = a list with the atoms we want to move that belong to the SECOND molecule
def swap_atoms_symmetrically(system, at_1, at_2, atoms_to_move_1, atoms_to_move_2):
	x_y_cords_1 = [system[at_2].position[i] for i in range(0, 2)]
	x_y_cords_2 = [system[at_1].position[i] for i in range(0, 2)]
	for i in range(0, 2):
		system[at_1].position[i] = x_y_cords_1[i]
		system[at_2].position[i] = x_y_cords_2[i]

	moving_atoms_1 = [i for i in atoms_to_move_1 if i != at_1]
	moving_atoms_2 = [i for i in atoms_to_move_2 if i != at_2]
	for i in moving_atoms_1:
		system[i].position[0] += system[at_1].position[0] - x_y_cords_2[0]
		system[i].position[1] += system[at_1].position[1] - x_y_cords_2[1]
	for i in moving_atoms_2:
		system[i].position[0] += system[at_2].position[0] - x_y_cords_1[0]
		system[i].position[1] += system[at_2].position[1] - x_y_cords_1[1]
	return system


#system = the strucwe use
#atom = Index of the atom which we want to add a new atom symmetrically
#bonded_to = Index of the atom where the "atom" above is bonded to
#to_bond = Index of the atom to which the new atom will be bonded.
def add_atom_symmetrically(system, atom, bonded_to, to_bond):
	pos_atom = list()
	pos_bonded_to = list()
	pos_to_bond = list()
	delta_pos = {}
	type_of_atom =  system[atom].symbol
	keys = ["x", "y", "z"]
	for i, j in enumerate(keys):
		pos_atom.append( system[ atom ].position[i] )
		pos_bonded_to.append( system[ bonded_to ].position[i] )
		pos_to_bond.append( system[ to_bond ].position[i] )
		delta_pos[ j ] = pos_to_bond[i] - pos_bonded_to[i]
	print(delta_pos)
	add_atoms(system, type_of_atom, system[ atom ].position[0] + delta_pos["x"], system[ atom ].position[1] + delta_pos["y"], system[ atom ].position[2] + delta_pos["z"] )
	return( system )
