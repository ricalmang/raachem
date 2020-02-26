#############################HOOK on ERRORS###############################
import traceback, sys
def show_exception_and_exit(exc_type, exc_value, tb):
	traceback.print_exception(exc_type, exc_value, tb)
	input("Press key to exit.")
	sys.exit(-1)
sys.excepthook = show_exception_and_exit
############################RAACHEM IMPORTS###############################
from raachem.file_creator.e_analysis import e_analysis, rel_scf, csv_e_analysis, deduplicate
from raachem.file_creator.input import CreateInputs, xyz_insert, validate_input
from raachem.file_creator.svg import vector_graph
from raachem.file_creator.xyz import input_to_xyz, log_to_xyz, log_to_xyz_scan, log_freq_xyz, superimpose_alg, geodes_int, xyz_ent
from raachem.file_creator.deploy_scripts import deploy
from raachem.util.gen_purp import file_weeder, mv_up_folder, preferences
try:import raapbs; pbs = True
except: pbs = False
###############################INTERFACE##################################
#TODO orca interface

def user_interface_1():#INITIAL
	if preferences.comp_software == "orca":
		print("Chose an option:")
		print("0 - Exit")
		print(" ")
		print("2 - Create .inp files from .xyz files")
		if pbs: print("3 - Create .pbs files and a submission script from {} files".format(".inp"))
		else: print("3 - Validate {} files".format(".inp"))
		print(" ")
		print("5 - Insert .xyz geometries into the corresponding .inp files")
		print(" ")
		print("7 - Configure")
		print("8 - Deploy script")
		print("9 - More options")
		option=input()
		if   option == "0": exit(print("ok"))
		elif option == "1": pass
		elif option == "2": CreateInputs()
		elif option == "3":	validate_input(file_weeder(["inp"]))
		elif option == "4":	pass
		elif option == "5":	xyz_insert(file_weeder([".xyz"]))
		elif option == "6": pass
		elif option == "7":	preferences.set_variables()
		elif option == "8": deploy()
		elif option == "9":	user_interface_2()
		else: print("Invalid input. Could not process request!")
		user_interface_1()
	elif preferences.comp_software == "gaussian":
		print("Chose an option:")
		print("0 - Exit")
		print("1 - Create .xyz files from .log files")
		print("2 - Create {} files from .xyz files".format(preferences.gauss_ext))
		if pbs: print("3 - Create .pbs files and a submission script from {} files".format(preferences.gauss_ext))
		else: print("3 - Validate {} files".format(preferences.gauss_ext))
		print("4 - IRC or SCAN analysis")
		print("5 - Insert .xyz geometries into the corresponding {} files".format(preferences.gauss_ext))
		print("6 - Extract free energies from .log files")
		print("7 - Configure")
		print("8 - Deploy script")
		print("9 - More options")
		option=input()
		if   option == "0": exit(print("ok"))
		elif option == "1": log_to_xyz(file_weeder([".log"]))
		elif option == "2": CreateInputs()
		elif option == "3":	validate_input(file_weeder([preferences.gauss_ext]))
		elif option == "4":	log_to_xyz_scan(file_weeder([".log"]))
		elif option == "5":	xyz_insert(file_weeder([".xyz"]))
		elif option == "6": e_analysis(file_weeder([".log"]))
		elif option == "7":	preferences.set_variables()
		elif option == "8": deploy()
		elif option == "9":	user_interface_2()
		else: print("Invalid input. Could not process request!")
		user_interface_1()
def user_interface_2():#MORE OPTIONS
	if preferences.comp_software == "orca":
		print("Chose an option:")
		print("0 - Go back to previous menu")
		print("1 - Create E_profile.svg")
		print("2 - Create enantiomer of xyz files")
		print("3 - Move Files up a folder")
		print(" ")
		print("5 - Superimpose two xyz files")
		print("6 - Generate solvated .xyz's")
		print(" ")
		print("8 - Create .xyz files from .inp files")
		print(" ")
		option=input()
		if   option == '0': user_interface_1()
		elif option == '1': vector_graph()
		elif option == "2":	xyz_ent()
		elif option == "3":	mv_up_folder()
		elif option == "4":	pass
		elif option == "5":	superimpose_alg()
		elif option == "6": geodes_int()
		elif option == "7":	pass
		elif option == "8": input_to_xyz()
		elif option == "9":	pass
		else: print("Invalid input. Could not process request!")
		user_interface_2()
	elif preferences.comp_software == "gaussian":
		print("Chose an option:")
		print("0 - Go back to previous menu")
		print("1 - Create E_profile.svg")
		print("2 - Create enantiomer of xyz files")
		print("3 - Move Files up a folder")
		print("4 - Split frequency in two directions and create corresponding .xyzs")
		print("5 - Superimpose two xyz files")
		print("6 - Generate solvated .xyz's")
		print("7 - Deduplicate .log files")
		print("8 - Create .xyz files from {} files".format(preferences.gauss_ext))
		print("9 - Analyze log files in current and sub directories")
		option=input()
		if   option == '0': user_interface_1()
		elif option == '1': vector_graph()
		elif option == "2":	xyz_ent()
		elif option == "3":	mv_up_folder()
		elif option == "4":	log_freq_xyz(file_weeder([".log"]))
		elif option == "5":	superimpose_alg()
		elif option == "6": geodes_int()
		elif option == "7":	deduplicate()
		elif option == "8": input_to_xyz()
		elif option == "9":	csv_e_analysis()
		else: print("Invalid input. Could not process request!")
		user_interface_2()
user_interface_1()

