#############################HOOK on ERRORS###############################
import traceback, sys
def show_exception_and_exit(exc_type, exc_value, tb):
	traceback.print_exception(exc_type, exc_value, tb)
	input("Press key to exit.")
	sys.exit(-1)
sys.excepthook = show_exception_and_exit
############################RAACHEM IMPORTS###############################
from raachem.file_creator.e_analysis import e_analysis, rel_scf, csv_e_analysis, deduplicate
from raachem.file_creator.gjf import gjf_gen, xyz_insert
from raachem.file_creator.pbs import pbs_creator
from raachem.file_creator.svg import vector_graph
from raachem.file_creator.xyz import gjf_to_xyz, log_to_xyz, log_to_xyz_scan, log_freq_xyz, superimpose_alg, geodes_int
from raachem.file_creator.deploy_scripts import deploy
from raachem.util.gen_purp import file_weeder, mv_up_folder, Var
###############################INTERFACE##################################
def user_interface_1():#INITIAL
	while True:
		print("Chose an option:")
		print("0 - Exit")
		print("1 - Create .xyz files from .log files")
		print("2 - Create .gjf files from .xyz files")
		print("3 - Create .pbs files and a submission script from .gjf files")
		print("4 - IRC or SCAN analysis")
		print("5 - Insert .xyz geometries into the corresponding .gjf files")
		print("6 - Extract free energies from .log files")
		print("7 - Configure")
		print("8 - Deploy script")
		print("9 - More options")
		option=input()
		if   option == "0": exit(print("ok"))
		elif option == "1": log_to_xyz(file_weeder([".log"]))
		elif option == "2": gjf_gen(file_weeder([".xyz"]))
		elif option == "3":	pbs_creator()
		elif option == "4":	log_to_xyz_scan(file_weeder([".log"]))
		elif option == "5":	xyz_insert(file_weeder([".xyz"]))
		elif option == "6": e_analysis(file_weeder([".log"]))
		elif option == "7": Var().set_variables()
		elif option == "8": deploy()
		elif option == "9":	user_interface_2()
		else: print("Invalid input. Could not process request!")
def user_interface_2():#MORE OPTIONS
	while True:
		print("Chose an option:")
		print("0 - Go back to previous menu")
		print("1 - Create E_profile.svg")
		print("2 - SCF energy analysis for optimizations")
		print("3 - Move Files up a folder")
		print("4 - Split frequency in two directions and create corresponding .xyzs")
		print("5 - Superimpose two xyz files")
		print("6 - Generate solvated .xyz's")
		print("7 - Deduplicate .log files")
		print("8 - Create .xyz files from .gjf files")
		print("9 - Analyze log files in current and sub directories")
		option=input()
		if   option == '0': user_interface_1()
		elif option == '1': vector_graph()
		elif option == "2":	rel_scf()
		elif option == "3":	mv_up_folder()
		elif option == "4":	log_freq_xyz(file_weeder([".log"]))
		elif option == "5":	superimpose_alg()
		elif option == "6": geodes_int()
		elif option == "7":	deduplicate()
		elif option == "8": gjf_to_xyz(file_weeder([".gjf"]))
		elif option == "9":	csv_e_analysis()
		else: print("Invalid input. Could not process request!")
user_interface_1()
