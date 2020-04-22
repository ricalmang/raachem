############################RAACHEM IMPORTS###############################
from raachem.file_creator.e_analysis import e_analysis, csv_e_analysis, deduplicate
from raachem.file_creator.input import CreateInputs, xyz_insert, validate_input
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
		self.menu_b = [self.func_10,self.func_12,self.func_13,
					   self.func_14,self.func_15,self.func_16,
					   self.func_17,self.func_18,self.func_20,
					   self.func_19,]
		try:
			import raapbs
			self.pbs = True
		except ImportError:
			self.pbs = False
	def describe(self):
		"""Call help for all functions"""
		for a in dir(self):
			if a.startswith("func"):
				help("raachem.util.interface.UserInterface.{}".format(a))
	def render(self,menu=None):
		"""Render menu"""
		if menu is None: menu = self.menu_a
		while True:
			if preferences.comp_software == "gaussian":
				print("Choose an option:")
			else:
				print("Choose an option (option* = functions for gaussian software files):")
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
			if mode == "print":	print("{} - Create .xyz files from .log files*".format(idx))
			elif mode == "run":	log_to_xyz(file_weeder([".log"]))
		elif preferences.comp_software == "gaussian":
			if mode == "print":	print("{} - Create .xyz files from .log files".format(idx))
			elif mode == "run": log_to_xyz(file_weeder([".log"]))
	@staticmethod
	def func_02(mode="print",idx=None):
		"""Create input files from .xyz files"""
		if preferences.comp_software == "orca":
			if mode == "print":	print("{} - Create .inp files from .xyz files".format(idx))
			elif mode == "run":	CreateInputs()
		elif preferences.comp_software == "gaussian":
			if mode == "print":	print("{} - Create {} files from .xyz files".format(idx,preferences.gauss_ext))
			elif mode == "run":	CreateInputs()
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
			if mode == "print":	print("{} - IRC, OPT or SCAN analysis*".format(idx))
			elif mode == "run":	log_to_xyz_scan(file_weeder([".log"]))
		elif preferences.comp_software == "gaussian":
			if mode == "print": print("{} - IRC, OPT or SCAN analysis".format(idx))
			elif mode == "run":	log_to_xyz_scan(file_weeder([".log"]))
	@staticmethod
	def func_05(mode="print",idx=None):
		"""Insert xyz geometries into input"""
		if preferences.comp_software == "orca":
			if mode == "print":	print("{} - Insert .xyz geometries into the corresponding .inp files".format(idx))
			elif mode == "run":	xyz_insert(file_weeder([".xyz"]))
		elif preferences.comp_software == "gaussian":
			if mode == "print":	print("{} - Insert .xyz geometries into the corresponding {} files".format(idx,preferences.gauss_ext))
			elif mode == "run":	xyz_insert(file_weeder([".xyz"]))
	@staticmethod
	def func_06(mode="print",idx=None):
		"""Analyze relative free energies"""
		if preferences.comp_software == "orca":
			if mode == "print":	print("{} - Extract free energies from .log files*".format(idx))
			elif mode == "run": e_analysis(file_weeder([".log"]))
		elif preferences.comp_software == "gaussian":
			if mode == "print":	print("{} - Extract free energies from .log files".format(idx))
			elif mode == "run": e_analysis(file_weeder([".log"]))
	@staticmethod
	def func_07(mode="print",idx=None):
		"""Create inputs from .log files"""
		if preferences.comp_software == "orca":
			if mode == "print":	print("{} - Create .inp files from .log files".format(idx))
			elif mode == "run":	CreateInputs(use_logs=True)
		elif preferences.comp_software == "gaussian":
			if mode == "print":	print("{} - Create {} files from .log files".format(idx,preferences.gauss_ext))
			elif mode == "run":	CreateInputs(use_logs=True)
	def func_08(self,mode="print",idx=None):
		"""Placeholder"""
		_ , _ , _ = self, mode, idx
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
		if mode == "print":	pass
		elif mode == "run":	pass
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
		"""Create to xyzs from a given vibration"""
		if preferences.comp_software == "orca":
			if mode == "print":
				print("{} - Split frequency in two directions and create corresponding .xyzs*".format(idx))
			elif mode == "run":
				log_freq_xyz(file_weeder([".log"]))
		elif preferences.comp_software == "gaussian":
			if mode == "print":
				print("{} - Split frequency in two directions and create corresponding .xyzs".format(idx))
			elif mode == "run":
				log_freq_xyz(file_weeder([".log"]))
	@staticmethod
	def func_15(mode="print",idx=None):
		"""Superimpose xyz structures"""
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
			if mode == "print":
				print("{} - Deduplicate .log files*".format(idx))
			elif mode == "run":
				deduplicate()
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
			if mode == "print":
				print("{} - Analyze log files in current and sub directories*".format(idx))
			elif mode == "run":
				csv_e_analysis()
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

