import os, math, shutil, re
from difflib import SequenceMatcher
from raachem.file_class.xyz import XyzFile
from raachem.file_class.gjf import GjfFile
from raachem.file_class.inp import InpFile
from raachem.util.constants import element_radii, keywords,cf, elements
from raachem.util.gen_purp import file_weeder, read_item, preferences, sel_files

class CreateInputs:
	p_files_dir = os.path.join(os.path.dirname(__file__), "par_templates")
	def __init__(self):
		self.template_ext = ".GAUSSIAN" if preferences.comp_software == "gaussian" else ".ORCA"
		self.extension = preferences.gauss_ext if preferences.comp_software == "gaussian" else ".inp"
		self.weeded_list = file_weeder(["xyz"])
		self.submit_parameters()
	def submit_parameters(self,p_files_dir=p_files_dir):
		p_files = file_weeder([self.template_ext], cf=p_files_dir)
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
				option = {str(i): str(i) for i in range(2)}.get(input(), True)
			if option == "0": return
			elif option == "1" and p_files:
				p_files.insert(0, None)
				option_2 = True
				print("Reading {} files from:\n{}\nChose a parameter file:".format(self.template_ext,p_files_dir))
				for idx, i in enumerate(p_files):
					if i == None: print("  0 - Cancel")
					else: print("{:>3} - {}".format(idx, i))
				while option_2 == True: option_2 = {str(idx): i for idx, i in enumerate(p_files)}.get(input(), True)
				if option_2 == None: return
				else:
					shutil.copy(os.path.join(p_files_dir, option_2), os.path.join(cf, "PARAMETERS{}".format(self.template_ext)))
					self.try_again()
			else:
				print("No builtins")
	def try_again(self):
		print("A PARAMETERS {} file was added to the current folder!".format(self.template_ext))
		print("You can edit it now")
		print("0 - To cancel")
		print("1 - File is ready!")
		option = True
		while option == True: option = {str(i): str(i) for i in range(2)}.get(input(), True)
		if option == "0": return
		elif option == "1":	self.__init__()
	def find_idx(self,parameters):
		lookup_str = "=====WILL BE GEOMETRY BLOCK====="
		for idx, line in enumerate(parameters[1:]):
			if line.strip() == lookup_str and self.template_ext == ".GAUSSIAN": self.save_gjf(parameters[1:], idx); return
			if line.strip() == lookup_str and self.template_ext == ".ORCA": pass; return
		print("The following line is required on the {} file:".format(self.template_ext))
		print(lookup_str)
		return
	def save_inp(self,parameters,index):
		if not preferences.folder_op: self.weeded_list = sel_files(self.weeded_list)
		if self.weeded_list == False: return
		for i in self.weeded_list:
			if os.path.isfile(os.path.join(cf, (i.replace(".xyz", ".inp")))):
				print(i.replace(".xyz",".inp") + " already exist on current directory!")
				continue
			xyz = XyzFile(read_item(i))
			inp_out = []
			for line in parameters[0:index]:
				inp_out.append(line.replace("FILENAME", i.replace(".xyz", "")) + "\n")
			for line in xyz.form_cord_block():
				inp_out.append(line + "\n")
			for line in parameters[index + 1:]:
				inp_out.append(line + "\n")
			with open(os.path.join(cf, (i.replace(".xyz",".inp"))), "w") as file:
				for line in inp_out: file.write(line)
			print(i.replace(".xyz", ".inp"), " created!")
		return
	def save_gjf(self,parameters,index):
		heavy_e, gjf_overwrite, folder_op = [preferences.heavy_atom,preferences.gjf_overwrite,preferences.folder_op]
		if not folder_op: self.weeded_list = sel_files(self.weeded_list)
		if self.weeded_list == False: return
		for i in self.weeded_list:
			if os.path.isfile(os.path.join(cf,(i.replace(".xyz",preferences.gauss_ext)))) and not gjf_overwrite:
				print(i.replace(".xyz",preferences.gauss_ext) + " already exist on current directory!")
				continue
			xyz = XyzFile(read_item(i))
			gjf_out=[]
			rm_lines = []
			for line in parameters[0:index]:
				gjf_out.append(line.replace("FILENAME",i.replace(".xyz",""))+"\n")
			for line in xyz.form_cord_block():
				gjf_out.append(line+"\n")
			for line in parameters[index+1:]:
				# SCAN like: "B S APROX" or "B S DIST"
				if len(line.split()) == 3 and any("modredundant" in i.lower() for i in parameters[0:index]):
					if xyz.n_atoms() <2:
						raise Exception("At least two atoms are neded in struture {} to perform a scan".format(xyz.name()))
					atom=["a","b"]
					possible = [str(b+1) for b in range(xyz.n_atoms())]
					if line.split() == ["B","S","APROX"]:
						print("Detected bond scan (APROX) for xyz file {}.".format(i))
						while not all(atom[a] in possible if atom[0] != atom[1] else False for a in [0,1]):
							atom = [input("Enter atom A: "),input("Enter atom B: ")]
						line=" ".join(["B", atom[0], atom[1], "S","APROX"])
					elif line.split() == ["B", "S", "DIST"]:
						print("Detected bond scan (DIST) for xyz file {}.".format(i))
						while not all(atom[a] in possible if atom[0] != atom[1] else False for a in [0,1]):
							atom = [input("Enter atom A: "),input("Enter atom B: ")]
						line =" ".join(["B", atom[0], atom[1], "S", "DIST"])
				# SCAN like: "B 1 2 S APROX" or "B 1 2 S DIST"
				if len(line.split()) == 5 and any("modredundant" in i.lower() for i in parameters[0:index]):
					if all([line.split()[0] == "B",line.split()[1].isnumeric(),line.split()[2].isnumeric(), line.split()[3] == "S"]):
						atoms_idx = [int(i)-1 for idx,i in enumerate(line.split()) if idx in [1,2]]
						if xyz.n_atoms() < 2:
							raise Exception("At least two atoms are neded in struture {} to perform a scan".format(xyz.name()))
						if any(True for a in atoms_idx if not a in range(xyz.n_atoms())):
							raise Exception("Scan atom numbers are larger than the number of atoms for: {}".format(xyz.name()))
						atoms_cord = [i for idx,i in enumerate(xyz.cord_block()) if idx in atoms_idx]
						dist = math.sqrt(sum((float(i)-float(atoms_cord[1][idx]))**2 for idx,i in enumerate(atoms_cord[0]) if idx > 0))
						ideal_dist=sum(b[1] for idx,b in enumerate(element_radii) if b[0] in [i[0] for i in atoms_cord])/100
						if all([line.split()[0] == "B", line.split()[3] == "S",line.split()[-1] == "APROX"]):
							gjf_out.append(line.replace("APROX",str(1+int((dist-ideal_dist)/0.075))+" -0.075\n"))
							continue
						elif all([line.split()[0] == "B", line.split()[3] == "S",line.split()[-1] == "DIST"]):
							gjf_out.append(line.replace("DIST",str(1+int(ideal_dist/0.075))+" 0.075\n"))
							continue
				elif len(line.split()) == 2:
					if any(True for a in ["/gen","gen ","genecp"] if a in " ".join(parameters[0:index-3]).lower()):
						if line.split() == ["LIGHT_ELEMENT_BASIS","0"]:
							elm = [a for a in xyz.elements() if elements.index(a) < heavy_e+1]
							if elm: line = line.replace("LIGHT_ELEMENT_BASIS"," ".join(elm))
							else:
								for a in range(3): rm_lines.append(len(gjf_out)+a)
						elif line.split() == ["HEAVY_ELEMENT_BASIS","0"]:
							elm = [a for a in xyz.elements() if elements.index(a) > heavy_e]
							if elm: line = line.replace("HEAVY_ELEMENT_BASIS"," ".join(elm))
							else:
								for a in range(3): rm_lines.append(len(gjf_out)+a)
					if "pseudo=read" in " ".join(parameters[0:index - 3]).lower().replace(" ",""):
						if line.split() == ["HEAVY_ELEMENT_ECP","0"]:
							elm = [a for a in xyz.elements() if elements.index(a) > heavy_e]
							if elm: line = line.replace("HEAVY_ELEMENT_ECP"," ".join(elm))
							else:
								pattern = re.compile(r'pseudo.{0,3}=.{0,3}read',re.IGNORECASE)
								eval_lines = gjf_out[0:index-3]
								for idx_a,a in enumerate(eval_lines):
									gjf_out[idx_a] = pattern.sub("",a)
								rm_lines.append(len(gjf_out))
								rm_lines.append(len(gjf_out)+1)
				gjf_out.append(line.replace("FILENAME",i.replace(".xyz",""))+"\n")
			gjf_out.append("\n")
			with open(os.path.join(cf,(i.replace(".xyz",preferences.gauss_ext))),"w") as gjf_file:
				for line in [i for idx,i in enumerate(gjf_out) if idx not in rm_lines]: gjf_file.write(line)
			print(i.replace(".xyz",preferences.gauss_ext)," created!")
		return

def xyz_insert(weeded_list):
	"""Inserts geometries into both orca and gaussian input files"""
	extension = ".inp" if preferences.comp_software == "orca" else preferences.gauss_ext
	sub_folder = os.path.join(cf,"inserted_input_files")
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

def validate_gjf(weeded_list):
	print("---------------------------------------------------------------------------")
	print("{:^30}{:^15}{:^10}{:^10}{:^10}".format("File","e- number","charge","multip","Validated"))
	print("---------------------------------------------------------------------------")
	novel_keys = []
	for item in weeded_list:
		gjf = GjfFile(read_item(item))
		split_list = [i for i in gjf.list[1:gjf.title_idx()] if not i.lower().startswith("%chk")]
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
		print(" {:<30}{:>10}{:>10}{:>10}{:^15}".format(gjf.name(),gjf.n_electrons(),gjf.charge(),gjf.multiplicity(),gjf.c_m_validate_txt()))
		if any([gjf.basis_errors(),typo_error,gjf.ecp_errors(preferences.heavy_atom),len(gjf.name().split()) != 1]):
			print("{:>7}+----------------------------ALERT---------------------------+".format(" "))
			for error in gjf.basis_errors(): print("{:>8}{:>60}".format("|",error+" |"))
			for error in gjf.ecp_errors(preferences.heavy_atom): print("{:>8}{:>60}".format("|",error+" |"))
			for i in typo_error: print("{:>8}{:>60}".format("|", "Is '{}' a typo of '{}'?".format(i[0],i[1]) + " |"))
			if len(gjf.name().split()) != 1: print("{:>8}{:>60}".format("|", "Filename must not contain spaces!" + " |"))
			print("{:>7}+-----------------------------------------------------------+".format(" "))
	print("---------------------------------------------------------------------------\n")
	novel_keys = list(dict.fromkeys([a for b in novel_keys for a in b]))
	if novel_keys:
		print("The following keys were not recognized:")
		print(novel_keys)
		print("---------------------------------------------------------------------------\n")
	try: import raapbs; raapbs.option()
	except: return

def validate_inp(weeded_list):
	print("---------------------------------------------------------------------------")
	print("{:^30}{:^15}{:^10}{:^10}{:^10}".format("File","e- number","charge","multip","Validated"))
	print("---------------------------------------------------------------------------")
	for item in weeded_list:
		inp = InpFile(read_item(item))
		print(" {:<30}{:>10}{:>10}{:>10}{:^15}".format(inp.name,inp.n_electrons,inp.charge,inp.mult,inp.c_m_validate_txt))
	print("---------------------------------------------------------------------------\n")
	try: import raapbs; raapbs.option()
	except: return