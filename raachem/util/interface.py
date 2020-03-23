############################RAACHEM IMPORTS###############################
from raachem.file_creator.e_analysis import e_analysis, csv_e_analysis, deduplicate
from raachem.file_creator.input import CreateInputs, xyz_insert, validate_input
from raachem.file_creator.svg import vector_graph
from raachem.file_creator.xyz import input_to_xyz, log_to_xyz, log_to_xyz_scan, log_freq_xyz, superimpose_alg, geodes_int, xyz_ent
from raachem.file_creator.deploy_scripts import deploy
from raachem.util.gen_purp import file_weeder, mv_up_folder, preferences

###############################INTERFACE##################################

class UserInterface:
	def __init__(self):
		self.menu_a = [self.func_00,self.func_01,self.func_02,
					   self.func_03,self.func_04,self.func_05,
					   self.func_06,self.func_07,self.func_21,
					   self.func_09]
		self.menu_b = [self.func_10,self.func_11,self.func_12,
					   self.func_13,self.func_14,self.func_15,
					   self.func_16,self.func_17,self.func_18,
					   self.func_19,self.func_20]
		try:
			import raapbs
			self.pbs = True
		except ImportError:
			self.pbs = False
	def describe(self):
		for a in dir(self):
			if a.startswith("func"):
				help("raachem.util.interface.UserInterface.{}".format(a))
	def render(self,menu=None):
		if menu is None: menu = self.menu_a
		while True:
			for i,a in enumerate(menu):	a(mode="print",idx=i)
			while True:
				option = input()
				if option in (str(a) for a in range(len(menu))):break
				else: print("Invalid input. Could not process request!")
			menu[int(option)](mode="run")
	@staticmethod
	def func_00(mode="print",idx=None):
		"""Exits the program"""
		if mode == "print": print("{} - Exit".format(idx))
		elif mode == "run": exit(print("ok"))
	@staticmethod
	def func_01(mode="print",idx=None):
		"""Create .xyz files from log files"""
		if preferences.comp_software == "orca":
			if mode == "print":	print(" ")
			elif mode == "run":	pass
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Create .xyz files from .log files".format(idx))
			elif mode == "run":
				log_to_xyz(file_weeder([".log"]))
	@staticmethod
	def func_02(mode="print",idx=None):
		"""Create input files from .xyz files"""
		if preferences.comp_software == "orca":
			if mode == "print":
				print("{} - Create .inp files from .xyz files".format(idx))
			elif mode == "run":
				CreateInputs()
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Create {} files from .xyz files".format(idx,preferences.gauss_ext))
			elif mode == "run":
				CreateInputs()
	def func_03(self,mode="print",idx=None):
		"""Validade inputs and possibly create pbs files"""
		if preferences.comp_software == "orca":
			if mode == "print":
				if self.pbs: print("{} - Create .pbs files and a submission script from {} files".format(idx,".inp"))
				else: print("{} - Validate {} files".format(idx,".inp"))
			elif mode == "run":
				validate_input(file_weeder(["inp"]))
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				if self.pbs: print("{} - Create .pbs files and a submission script from {} files".format(idx,preferences.gauss_ext))
				else: print("{} - Validate {} files".format(idx,preferences.gauss_ext))
			elif mode == "run":
				validate_input(file_weeder([preferences.gauss_ext]))
	@staticmethod
	def func_04(mode="print",idx=None):
		"""Analyze scan, opt and irc files"""
		if preferences.comp_software == "orca":
			if mode == "print":	print(" ")
			elif mode == "run":	pass
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - IRC, OPT or SCAN analysis".format(idx))
			elif mode == "run":
				log_to_xyz_scan(file_weeder([".log"]))
	@staticmethod
	def func_05(mode="print",idx=None):
		"""Insert xyz geometries into input"""
		if preferences.comp_software == "orca":
			if mode == "print":
				print("{} - Insert .xyz geometries into the corresponding .inp files".format(idx))
			elif mode == "run":
				xyz_insert(file_weeder([".xyz"]))
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Insert .xyz geometries into the corresponding {} files".format(idx,preferences.gauss_ext))
			elif mode == "run":
				xyz_insert(file_weeder([".xyz"]))
	@staticmethod
	def func_06(mode="print",idx=None):
		"""Analyze relative free energies"""
		if preferences.comp_software == "orca":
			if mode == "print":	print(" ")
			elif mode == "run": pass
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Extract free energies from .log files".format(idx))
			elif mode == "run":
				e_analysis(file_weeder([".log"]))
	@staticmethod
	def func_07(mode="print",idx=None):
		"""Create inputs from .log files"""
		if preferences.comp_software == "orca":
			if mode == "print":
				print("{} - Create .inp files from .log files".format(idx))
			elif mode == "run":
				CreateInputs(use_logs=True)
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Create {} files from .log files".format(idx,preferences.gauss_ext))
			elif mode == "run":
				CreateInputs(use_logs=True)
	def func_08(self,mode="print",idx=None):
		"""Analyze optimizations"""
		_ , _ , _ = self, mode, idx
		pass
	def func_09(self,mode="print",idx=None):
		"""More options"""
		if mode == "print":	print("{} - More options".format(idx))
		elif mode == "run":	self.render(self.menu_b)
	def func_10(self,mode="print",idx=None):
		"""Go back to previous menu"""
		if mode == "print":	print("{} - Go back to previous menu".format(idx))
		elif mode == "run":	self.render(self.menu_a)
	@staticmethod
	def func_11(mode="print",idx=None):
		"""Create energy profile as .svg"""
		if mode == "print":	print("{} - Create E_profile.svg".format(idx))
		elif mode == "run":	vector_graph()
	@staticmethod
	def func_12(mode="print",idx=None):
		"""Return enantiomer of xyz file"""
		if mode == "print":	print("{} - Create enantiomer of xyz files".format(idx))
		elif mode == "run":	xyz_ent()
	@staticmethod
	def func_13(mode="print",idx=None):
		"""Move files up a folder"""
		if mode == "print":	print("{} - Move Files up a folder".format(idx))
		elif mode == "run":	mv_up_folder()
	@staticmethod
	def func_14(mode="print",idx=None):
		"""Create inputs from .log files"""
		if preferences.comp_software == "orca":
			if mode == "print":	print("")
			elif mode == "run": pass
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Split frequency in two directions and create corresponding .xyzs".format(idx))
			elif mode == "run":
				log_freq_xyz(file_weeder([".log"]))
	@staticmethod
	def func_15(mode="print",idx=None):
		"""Move files up a folder"""
		if mode == "print":	print("{} - Superimpose two xyz files".format(idx))
		elif mode == "run":	superimpose_alg()
	@staticmethod
	def func_16(mode="print",idx=None):
		"""Generate solvated .xyzs"""
		if mode == "print":	print("{} - Generate solvated .xyz's".format(idx))
		elif mode == "run":	geodes_int()
	@staticmethod
	def func_17(mode="print",idx=None):
		"""Deduplicate .log files"""
		if preferences.comp_software == "orca":
			if mode == "print":	print("")
			elif mode == "run": pass
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Deduplicate .log files".format(idx))
			elif mode == "run":
				deduplicate()
	@staticmethod
	def func_18(mode="print",idx=None):
		"""Create .xyz files from input files"""
		if preferences.comp_software == "orca":
			if mode == "print":
				print("{} - Create .xyz files from .inp files".format(idx))
			elif mode == "run":
				input_to_xyz()
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Create .xyz files from {} files".format(idx,preferences.gauss_ext))
			elif mode == "run":
				input_to_xyz()
	@staticmethod
	def func_19(mode="print",idx=None):
		"""Analyze log files in current and sub directories"""
		if preferences.comp_software == "orca":
			if mode == "print":	print("")
			elif mode == "run": pass
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Analyze log files in current and sub directories".format(idx))
			elif mode == "run":
				csv_e_analysis()
	@staticmethod
	def func_20(mode="print",idx=None):
		"""Configure"""
		if mode == "print":	print("{} - Configure".format(idx))
		elif mode == "run":	preferences.set_variables()
	@staticmethod
	def func_21(mode="print",idx=None):
		"""Deploy scripts"""
		if mode == "print":	print("{} - Deploy script".format(idx))
		elif mode == "run":	deploy()

#def user_interface_1():#INITIAL
#	if preferences.comp_software == "orca":
		#print("Chose an option:")
		#print("0 - Exit")
		#print(" ")
		#print("2 - Create .inp files from .xyz files")
		#if pbs: print("3 - Create .pbs files and a submission script from {} files".format(".inp"))
		#else: print("3 - Validate {} files".format(".inp"))
		#print(" ")
		#print("5 - Insert .xyz geometries into the corresponding .inp files")
		#print(" ")
		#print("7 - Create .inp files from .log files")

		#print("9 - More options")
		#option=input()
		#if   option == "0": exit(print("ok"))
		#elif option == "1": pass
		#elif option == "2": CreateInputs()
		#elif option == "3":	validate_input(file_weeder(["inp"]))
		#elif option == "4":	pass
		#elif option == "5":	xyz_insert(file_weeder([".xyz"]))
		#elif option == "6": pass
		#elif option == "7":	CreateInputs(use_logs=True)

		#elif option == "9":	user_interface_2()
		#else: print("Invalid input. Could not process request!")
		#user_interface_1()
#	elif preferences.comp_software == "gaussian":
		#print("Chose an option:")
		#print("0 - Exit")
		#print("1 - Create .xyz files from .log files")
		#print("2 - Create {} files from .xyz files".format(preferences.gauss_ext))
		#if pbs: print("3 - Create .pbs files and a submission script from {} files".format(preferences.gauss_ext))
		#else: print("3 - Validate {} files".format(preferences.gauss_ext))
		#print("4 - IRC or SCAN analysis")
		#print("5 - Insert .xyz geometries into the corresponding {} files".format(preferences.gauss_ext))
		#print("6 - Extract free energies from .log files")
		#print("7 - Create {} files from .log files".format(preferences.gauss_ext))

#		print("9 - More options")
#		option=input()
		#if   option == "0": exit(print("ok"))
		#elif option == "1": log_to_xyz(file_weeder([".log"]))
		#elif option == "2": CreateInputs()
		#elif option == "3":	validate_input(file_weeder([preferences.gauss_ext]))
		#elif option == "4":	log_to_xyz_scan(file_weeder([".log"]))
		#elif option == "5":	xyz_insert(file_weeder([".xyz"]))
		#elif option == "6": e_analysis(file_weeder([".log"]))
		#elif option == "7":	CreateInputs(use_logs=True)

#		elif option == "9":	user_interface_2()
#		else: print("Invalid input. Could not process request!")
#		user_interface_1()
#def user_interface_2():#MORE OPTIONS
#	if preferences.comp_software == "orca":
#		print("Chose an option:")
		#print("0 - Go back to previous menu")
		#print("1 - Create E_profile.svg")
		#print("2 - Create enantiomer of xyz files")
		#print("3 - Move Files up a folder")
		#print(" ")
		#print("5 - Superimpose two xyz files")
		#print("6 - Generate solvated .xyz's")
		#print(" ")
		#print("8 - Create .xyz files from .inp files")
	#	print(" ")
#		print("10 - Configure")
#		print("11 - Deploy script")
		#option=input()
		#if   option == '0': user_interface_1()
		#elif option == '1': vector_graph()
		#elif option == "2":	xyz_ent()
		#elif option == "3":	mv_up_folder()
		#elif option == "4":	pass
		#elif option == "5":	superimpose_alg()
		#elif option == "6": geodes_int()
		#elif option == "7":	pass
#		elif option == "8": input_to_xyz()
#		elif option == "9":	pass
#		elif option == "10": preferences.set_variables()
#		elif option == "11": deploy()
#		else: print("Invalid input. Could not process request!")
#		user_interface_2()
#	elif preferences.comp_software == "gaussian":
#		print("Chose an option:")
#		print("0 - Go back to previous menu")
#		print("1 - Create E_profile.svg")
		#print("2 - Create enantiomer of xyz files")
		#print("3 - Move Files up a folder")
		#print("4 - Split frequency in two directions and create corresponding .xyzs")
		#print("5 - Superimpose two xyz files")
		#print("6 - Generate solvated .xyz's")
		#print("7 - Deduplicate .log files")
		#print("8 - Create .xyz files from {} files".format(preferences.gauss_ext))
		#print("9 - Analyze log files in current and sub directories")
#		print("10 - Configure")
#		print("11 - Deploy script")
#		option=input()
#		if   option == '0': user_interface_1()
#		elif option == '1': vector_graph()
		#elif option == "2":	xyz_ent()
		#elif option == "3":	mv_up_folder()
		#elif option == "4":	log_freq_xyz(file_weeder([".log"]))
		#elif option == "5":	superimpose_alg()
		#elif option == "6": geodes_int()
		#elif option == "7":	deduplicate()
		#elif option == "8": input_to_xyz()
		#elif option == "9":	csv_e_analysis()
#		elif option == "10": preferences.set_variables()
#		elif option == "11": deploy()
#		else: print("Invalid input. Could not process request!")
#		user_interface_2()


