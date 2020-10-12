import os, math, shutil, re
from difflib import SequenceMatcher
from raachem.file_class.xyz import XyzFile
from raachem.file_class.gjf import GjfFile
from raachem.file_class.log import LogFile
from raachem.file_class.inp import InpFile
from raachem.util.constants import element_radii, keywords, elements
from raachem.util.gen_purp import file_weeder, read_item, preferences, sel_files, is_str_float

class CreateInputs:
	def __init__(self,use_logs=False):
		self.use_logs = use_logs
		self.original_ext = ".log" if use_logs else ".xyz"
		self.template_ext = ".GAUSSIAN" if preferences.comp_software == "gaussian" else ".ORCA"
		self.extension = preferences.gauss_ext if preferences.comp_software == "gaussian" else ".inp"
		self.weeded_list = file_weeder([self.original_ext])
		self.p_files_dir = os.path.join(os.path.dirname(__file__),"builtin_templates")
		self.user_files_dir = os.path.join(os.path.dirname(__file__), "par_templates")
		self.submit_parameters()
	def submit_parameters(self):
		l_files = file_weeder([self.template_ext])
		if len(l_files) == 1:
			print("Reading parameters from {} in current folder".format(l_files[0]))
			parameters = read_item(l_files[0])
			self.find_idx(parameters)
		elif len(l_files) > 1:
			l_files.insert(0, None)
			option = True
			print("Chose a parameter file:")
			for idx, i in enumerate(l_files):
				if i is None: print("  0 - Cancel")
				else: print("{:>3} - {}".format(idx, i))
			while option == True: option = {str(idx): i for idx, i in enumerate(l_files)}.get(input(), True)
			if option is None: return
			else: parameters = read_item(option); self.find_idx(parameters)
		else:
			print("Could not find {} file on current file".format(self.template_ext))
			option = True
			while option == True:
				print("Chose an option:")
				print("0 - To cancel")
				print("1 - Use builtins")
				print("2 - User defined")
				option = {str(i): str(i) for i in range(3)}.get(input(), True)
			if option == "0": return
			else:
				folder = self.p_files_dir if option == "1" else self.user_files_dir
				p_files = file_weeder([self.template_ext], cf=folder)
				p_files.insert(0, None)
				print(f"Reading {self.template_ext} files from:\n{folder}\nChoose a parameter file:")
				option_2 = True
				for idx, i in enumerate(p_files):
					if i == None: print("  0 - Cancel")
					else: print("{:>3} - {}".format(idx, i))
				while option_2 == True: option_2 = {str(idx): i for idx, i in enumerate(p_files)}.get(input(), True)
				if option_2 == None: return
				else:
					shutil.copy(os.path.join(folder, option_2), os.path.join(os.getcwd(), "PARAMETERS{}".format(self.template_ext)))
					self.try_again()
	def try_again(self):
		print("A PARAMETERS{} file was added to the current folder!".format(self.template_ext))
		print("You can edit it now")
		print("0 - To cancel")
		print("1 - File is ready!")
		option = True
		while option == True: option = {str(i): str(i) for i in range(2)}.get(input(), True)
		if option == "0": return
		elif option == "1":	self.__init__(use_logs=self.use_logs)
	def find_idx(self,parameters):
		lookup_str = "=====WILL BE GEOMETRY BLOCK====="
		for idx, line in enumerate(parameters[1:]):
			if line.strip() == lookup_str and self.template_ext == ".GAUSSIAN": self.save_gjf(parameters[1:], idx); return
			if line.strip() == lookup_str and self.template_ext == ".ORCA": self.save_inp(parameters[1:], idx); return
		print("The following line is required on the {} file:".format(self.template_ext))
		print(lookup_str)
		return
	def save_inp(self,parameters,index):
		if not preferences.folder_op: self.weeded_list = sel_files(self.weeded_list)
		if not self.weeded_list: return
		for i in self.weeded_list:
			if os.path.isfile(os.path.join(os.getcwd(), (i.replace(self.original_ext, ".inp")))):
				print(i.replace(self.original_ext,".inp") + " already exist on current directory!")
				continue
			xyz = LogFile(read_item(i)).last_xyz_obj() if self.use_logs else XyzFile(read_item(i))
			inp_out = []
			for idx_a, line in enumerate(parameters[0:index]):
				if self.use_logs and idx_a + 1 == len(parameters[0:index]):
					line = "* " + " ".join(LogFile(read_item(i)).charge_mult)
				inp_out.append(line.replace("FILENAME", i.replace(self.original_ext, "")) + "\n")
			for line in xyz.form_cord_block():
				inp_out.append(line + "\n")
			for line in parameters[index + 1:]:
				inp_out.append(line + "\n")
			with open(os.path.join(os.getcwd(), (i.replace(self.original_ext,".inp"))), "w",newline="\n") as file:
				for line in inp_out: file.write(line)
			print(i.replace(self.original_ext, ".inp"), " created!")
		return
	def save_gjf(self,parameters,index):
		heavy_e, gjf_overwrite, folder_op = [preferences.heavy_atom,preferences.gjf_overwrite,preferences.folder_op]
		ecps, basis = None, None
		if not folder_op: self.weeded_list = sel_files(self.weeded_list)
		if not self.weeded_list: return
		if any([b in parameters for b in ["INSERT_GBS_BASIS","INSERT_GBS_ECP"]]):
			print("This parameter file requires an aditional gbs file.")
			print("Where should it be taken from?")
			print("0 - Current folder")
			print("1 - Builtins")
			print("2 - User defined")
			while True:
				option = input()
				if option in ["0","1","2"]: break
			if option == "0": files_dir = os.getcwd()
			elif option == "1": files_dir = self.p_files_dir
			elif option == "2": files_dir = self.user_files_dir
			gbs_files = file_weeder([".gbs"], cf=files_dir)
			print("Chosse one basis/ecp set\n0 - Cancel")
			for i, file in enumerate(gbs_files): print("{} - {}".format(i+1, file))
			while True:
				option = input()
				if option in [str(a) for a in range(len(gbs_files)+1)]: break
				print("Invalid option")
			if option == "0": return
			with open(os.path.join(files_dir,gbs_files[int(option)-1])) as file:
				gbs_file = file.read().splitlines()
			basis, ecps = self.read_gbs(gbs_file)
			assert type(basis) == dict, "Basis was not read"
			assert type(ecps) == dict, "Ecps were not read"
		for i in self.weeded_list:
			try:
				gjf_out, rm_lines, = [], []
				if os.path.isfile(os.path.join(os.getcwd(),(i.replace(self.original_ext,preferences.gauss_ext)))) and not gjf_overwrite:
					print(i.replace(self.original_ext,preferences.gauss_ext) + " already exist on current directory!")
					continue
				xyz = LogFile(read_item(i)).last_xyz_obj() if self.use_logs else XyzFile(read_item(i))
				for idx_a,line in enumerate(parameters[0:index]):
					if self.use_logs and idx_a+1 == len(parameters[0:index]): line = " ".join(LogFile(read_item(i)).charge_mult)
					gjf_out.append(line.replace("FILENAME",i.replace(self.original_ext,""))+"\n")
				for line in xyz.form_cord_block():
					gjf_out.append(line+"\n")
				possible = [str(b + 1) for b in range(xyz.n_atoms())]
				for line in parameters[index+1:]:
					# SCAN like: "B S APROX" or "B S DIST"
					is_mod_red = any("modredundant" in i.lower() for i in parameters[0:index])
					s_line = line.split()
					while len(s_line) == 3:
						if s_line[0] != "B" or s_line[1] != "S": break
						if s_line[2] not in ["APROX","DIST"]: break
						if not is_mod_red:
							print("Missing 'modredundant' keyword?")
							break
						if xyz.n_atoms() <2:
							raise Exception(f"At least two atoms are neded in struture {xyz.name()} to perform a scan")
						atom=["a","b"]
						if s_line[-1] == "APROX":
							print(f"Detected bond scan (APROX) for xyz file {i}.")
							while not all(atom[a] in possible if atom[0] != atom[1] else False for a in [0,1]):
								atom = [input("Enter atom A: "),input("Enter atom B: ")]
							line = f"B {atom[0]} {atom[1]} S APROX"
						elif s_line[-1] == "DIST":
							print(f"Detected bond scan (DIST) for xyz file {i}.")
							while not all(atom[a] in possible if atom[0] != atom[1] else False for a in [0,1]):
								atom = [input("Enter atom A: "),input("Enter atom B: ")]
							line = f"B {atom[0]} {atom[1]} S DIST"
						s_line = line.split()
						break
					# SCAN like: "B 1 2 S APROX" or "B 1 2 S DIST"
					while len(s_line) == 5:

						if s_line[0] != "B" or s_line[3] != "S": break
						if not s_line[1].isdigit(): break
						if not s_line[2].isdigit(): break
						if s_line[4] not in ["APROX", "DIST"]: break
						if not is_mod_red:
							print("Missing 'modredundant' keyword?")
							break
						if xyz.n_atoms() < 2:
							raise Exception(f"At least two atoms are neded in struture {xyz.name()} to perform a scan")
						atoms_idx = [int(s_line[1])-1, int(s_line[2])-1]
						if any(True for a in atoms_idx if not a in range(xyz.n_atoms())):
							raise Exception(f"Scan atom numbers are larger than the number of atoms for: {xyz.name()}")
						atoms_cord = [i for idx,i in enumerate(xyz.cord_block()) if idx in atoms_idx]
						dist = math.sqrt(sum((float(i)-float(atoms_cord[1][idx]))**2 for idx,i in enumerate(atoms_cord[0]) if idx > 0))
						ideal_dist=sum(b[1] for idx,b in enumerate(element_radii) if b[0] in [i[0] for i in atoms_cord])/100
						if s_line[-1] == "APROX":
							gjf_out.append(line.replace("APROX",f"{1+int((dist-ideal_dist)/0.075)} -0.075\n"))
						elif s_line[-1] == "DIST":
							gjf_out.append(line.replace("DIST",f"{1+int(ideal_dist/0.075)} 0.075\n"))
						s_line = line.split()
						break
					# SCAN like: "D S"
					while len(s_line) == 2:
						if not is_mod_red: break
						if s_line[0] != "D" or s_line[1] != "S": break
						print(f"Detected dihedral scan for xyz file {i}.")
						if xyz.n_atoms() < 4:
							raise Exception(f"At least two four atoms are neded in struture {xyz.name()} to perform a dihedral scan")
						while True:
							atom = [input(f"Enter atom {a}: ") for a in ["A","B","C","D"]]
							atom = [a.strip() for a in atom]
							if not all([len(a.split()) == 1 and a.isdigit() for a in atom]):
								print("No non-numeric characters or spaces are allowed for atom numbers")
								continue
							if not len(set(atom)) == 4:
								print("No atom should be repeated!")
								continue
							if not all([a in possible for a in atom]):
								print("Atoms with incorrect number were found!")
								continue
							break
						while True:
							n_steps = input("Please give the number of steps to be taken: ").strip()
							if n_steps.isdigit(): break
							else: print("The number of steps must be an integer")
						while True:
							print("Please give the step size to be taken in degrees")
							size = input("The number must include a dot character as a decimal place (eg. '10.0'): ").strip()
							if is_str_float(size) and size.count(".") == 1: break
							else: print("Please make sure the number complies with formating rules!")
						line = f"D {' '.join(atom)} S {n_steps} {size}"
						s_line = line.split()
						break
					# READOPTIMIZE + NOTATOMS
					if line.replace(" ", "").lower() == "notatoms=":
						print("Please enter the atoms you want to freeze on structure {} separated by comma".format(xyz.name()))
						line = line + input(line)
					# READOPTIMIZE + ATOMS
					elif line.replace(" ", "").lower() == "noatomsatoms=":
						print("Please enter the atoms you want to optimize on structure {} separated by comma".format(xyz.name()))
						line = line + input(line)
					# BASIS like: "LIGHT_ELEMENT_BASIS 0" or "HEAVY_ELEMENT_BASIS 0" and ECP like "HEAVY_ELEMENT_ECP 0"
					elif len(s_line) == 2:
						if any(True for a in ["/gen","gen ","genecp"] if a in " ".join(parameters[0:index-3]).lower()):
							if s_line == ["LIGHT_ELEMENT_BASIS","0"]:
								elm = [a for a in xyz.elements() if elements.index(a) < heavy_e+1]
								if elm: line = line.replace("LIGHT_ELEMENT_BASIS"," ".join(elm))
								else:
									for a in range(3): rm_lines.append(len(gjf_out)+a)
							elif s_line == ["HEAVY_ELEMENT_BASIS","0"]:
								elm = [a for a in xyz.elements() if elements.index(a) > heavy_e]
								if elm: line = line.replace("HEAVY_ELEMENT_BASIS"," ".join(elm))
								else:
									for a in range(3): rm_lines.append(len(gjf_out)+a)
						if "pseudo=read" in " ".join(parameters[0:index - 3]).lower().replace(" ",""):
							if s_line == ["HEAVY_ELEMENT_ECP","0"]:
								elm = [a for a in xyz.elements() if elements.index(a) > heavy_e]
								if elm: line = line.replace("HEAVY_ELEMENT_ECP"," ".join(elm))
								else:
									pattern = re.compile(r'pseudo.{0,3}=.{0,3}read',re.IGNORECASE)
									eval_lines = gjf_out[0:index-3]
									for idx_a,a in enumerate(eval_lines):
										gjf_out[idx_a] = pattern.sub("",a)
									rm_lines.append(len(gjf_out))
									rm_lines.append(len(gjf_out)+1)
					# BASIS like: "INSERT_GBS_BASIS"
					elif "INSERT_GBS_BASIS" in line:
						for element in sorted(xyz.elements(),key=lambda x: elements.index(x)):
							for line_a in basis[element]:
								gjf_out.append(line_a+"\n")
						continue
					# ECP like: "INSERT_GBS_ECP"
					elif "INSERT_GBS_ECP" in line:
						if "pseudo=read" in " ".join(parameters[0:index - 3]).lower().replace(" ", ""):
							need_ecp = [a for a in xyz.elements() if elements.index(a) > preferences.heavy_atom]
							for element in sorted(need_ecp,key=lambda x: elements.index(x)):
								for line_a in ecps[element]:
									gjf_out.append(line_a+"\n")
							if not need_ecp:
								pattern = re.compile(r'pseudo.{0,3}=.{0,3}read', re.IGNORECASE)
								eval_lines = gjf_out[0:index - 3]
								for idx_a, a in enumerate(eval_lines):
									gjf_out[idx_a] = pattern.sub("", a)
							continue
					gjf_out.append(line.replace("FILENAME",i.replace(self.original_ext,""))+"\n")
				gjf_out.append("\n")
				with open(os.path.join(os.getcwd(),(i.replace(self.original_ext,preferences.gauss_ext))),"w") as gjf_file:
					for line in [i for idx,i in enumerate(gjf_out) if idx not in rm_lines]: gjf_file.write(line)
				print(i.replace(self.original_ext,preferences.gauss_ext)," created!")
			except TypeError as e:
				print("Error on file {}".format(i))
				print(e)
		return
	@staticmethod
	def read_gbs(gbs_file):
		gbs = [a for a in gbs_file if not a.startswith("!")]
		gbs_starts = []
		for i,a in enumerate(gbs):
			if len(a.split()) != 2: continue
			if not a.split()[0].capitalize() in elements: continue
			if a.split()[1] != "0":continue
			else: gbs_starts.append(i)
		gbs_basis_ends = [i+1 for i,a in enumerate(gbs) if a.startswith("****")]
		basis_dict = {gbs[a].split()[0].capitalize():gbs[a:b] for a,b in zip(gbs_starts,gbs_basis_ends) if a < b}
		ecp_starts = [a for a in gbs_starts if a > max(gbs_basis_ends)]
		ecp_ends = []
		for start in ecp_starts:
			for i,line in enumerate(gbs[start:]):
				if i < 1:continue
				elif len(line.split()) == 0: ecp_ends.append(i+start); break
				elif i + 1 == len(gbs[start:]):  ecp_ends.append(i+start); break
				elif len(line.split()) != 2: continue
				elif line.split()[0].capitalize() in elements and line.split()[1] == "0": ecp_ends.append(i+start); break
		ecp_dict = {gbs[a].split()[0].capitalize():gbs[a:b] for a,b in zip(ecp_starts,ecp_ends) if a < b}
		return basis_dict,ecp_dict
def xyz_insert(weeded_list):
	"""Inserts geometries into both orca and gaussian input files"""
	extension = ".inp" if preferences.comp_software == "orca" else preferences.gauss_ext
	sub_folder = os.path.join(os.getcwd(),"inserted_input_files")
	if os.path.exists(sub_folder):
		print("'inserted_input_files' directory already exists in current directory!")
		print("Please remove it and try again!")
		return
	os.mkdir(sub_folder)
	for i in weeded_list:
		try:
			xyz = XyzFile(read_item(i))
			if preferences.comp_software == "orca": comp_input = InpFile(read_item(i.replace(".xyz",extension)))
			else: comp_input = GjfFile(read_item(i.replace(".xyz",extension)))
			comp_input = comp_input.replace_cord(xyz)
			with open(os.path.join(sub_folder,i.replace(".xyz",extension)),"w") as file:
				file.write(comp_input.return_print)
		except FileNotFoundError: print("file " + i.replace(".xyz",extension) + " could not be found!")
	print("\nJob done!\nPlease lookup the inserted_input_files directory\n")
	return

def validate_input(weeded_list):
	print("---------------------------------------------------------------------------")
	print("{:^30}{:^15}{:^10}{:^10}{:^10}".format("File","e- number","charge","multip","Validated"))
	print("---------------------------------------------------------------------------")
	novel_keys = []
	for item in weeded_list:
		if preferences.comp_software == "orca":	comp_input = InpFile(read_item(item))
		else: comp_input = GjfFile(read_item(item))
		a = (comp_input.name(), comp_input.n_electrons(), comp_input.charge(), comp_input.multiplicity(), comp_input.c_m_validate_txt())
		print(" {:<30}{:>10}{:>10}{:>10}{:^15}".format(*a))
		if preferences.comp_software == "gaussian":
			split_list = [i for i in comp_input.list[1:comp_input.title_idx()] if not i.lower().startswith("%chk")]
			for x in [None,"/","(",")",",","=","%",":"]:
				split_list = [a for b in [i.split(x) for i in split_list] for a in b if len(a) > 3]
			no_match = [i for i in split_list if i.lower() not in [j.lower() for j in keywords] and not i[0].isdigit()]
			typo_error = []
			for b in no_match:
				p_matches = [[b,i,SequenceMatcher(None,b,i).ratio()] for i in keywords if SequenceMatcher(None,b,i).ratio() > 0.6]
				if p_matches:
					p_matches.sort(key=lambda x: x[2])
					typo_error.append([p_matches[0][0],p_matches[0][1]])
			novel_keys.append([i for i in no_match if i not in [a[0] for a in typo_error]])
			if any([comp_input.basis_errors(),typo_error,comp_input.ecp_errors(preferences.heavy_atom),len(comp_input.name().split()) != 1]):
				print("{:>7}+----------------------------ALERT---------------------------+".format(" "))
				for error in comp_input.basis_errors(): print("{:>8}{:>60}".format("|",error+" |"))
				for error in comp_input.ecp_errors(preferences.heavy_atom): print("{:>8}{:>60}".format("|",error+" |"))
				for i in typo_error: print("{:>8}{:>60}".format("|", "Is '{}' a typo of '{}'?".format(i[0],i[1]) + " |"))
				if len(comp_input.name().split()) != 1: print("{:>8}{:>60}".format("|", "Filename must not contain spaces!" + " |"))
				print("{:>7}+-----------------------------------------------------------+".format(" "))
	print("---------------------------------------------------------------------------\n")
	novel_keys = list(dict.fromkeys([a for b in novel_keys for a in b]))
	if novel_keys and preferences.comp_software == "gaussian":
		print("The following keys were not recognized:")
		print(novel_keys)
		print("---------------------------------------------------------------------------\n")
	try:
		import raapbs
		raapbs.option()
	except ImportError:
		pass