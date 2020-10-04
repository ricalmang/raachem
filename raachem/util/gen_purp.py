import os, shutil, time, itertools


def read_item(file_name=None, promp=False, extension=None, cf=os.getcwd()):
	"""Reads an .xyz, .gjf, .com or .log item and returns a list of its contents ready for class instantiation"""
	if promp != False:
		if extension is None: extension = [".xyz"]
		print(promp)
		files = file_weeder(extension)
		if len(files) == 0: print("Sorry no such files in current directory!"); return []
		for idx,item in enumerate(files): print("{:<5}{:<25}".format(idx+1,item))
		while True:
			try: file_name = files[int(input())-1]; break
			except: print("Invalid input!")
	with open(os.path.join(cf,file_name),"r") as in_file: in_content = in_file.read().splitlines()
	in_content.insert(0,file_name)
	return in_content

def file_weeder(ext_to_weed,cf=os.getcwd(), promp=True):
	"""Looks up files with the extensions provided in current directory"""
	if type(ext_to_weed) == str: ext_to_weed = [ext_to_weed]
	weeded_list = []
	for a,b in itertools.product(ext_to_weed,os.listdir(cf)):
		if a in b: weeded_list.append(b)
	if promp and len(weeded_list) == 0:
		print("No '{}' files found in current directory!".format("' or '".join(ext_to_weed)))
		print(cf)
		return []
	return sorted(weeded_list)

def is_str_float(i):
	"""Check if a string can be converted into a float"""
	try: float(i); return True
	except ValueError: return False

def mv_up_folder():
	"""Move files with a chossen extension up a folder"""
	while True:
		extensions = [[0, "return"], ["1", ".log"], ["2", ".gjf"],["2", ".com"],["3", ".xyz"]]
		print("Select a file extension to move up a folder")
		for idx, ext in enumerate(extensions):
			if idx == 0: print("0   - Cancel request")
			else: print("{:<3} - {:>10}".format(ext[0], ext[1]))
		extension = {str(i[0]): i[1] for i in extensions}.get(input(), None)
		print(extension)
		if extension == "return": return
		if extension != None: break
		print("Invalid Option")
	folders = [x[0] for x in os.walk(os.getcwd())]
	for folder in folders[1:]:
		for file in file_weeder(extension, folder, promp=False):
			if os.path.isfile(os.path.join(folders[0], file)):
				print("Filename: {}\nFrom: {}\nAlready exists in uper directory:\n{}".format(file, folder, folders[0]))
				break
			shutil.copy2(os.path.join(folder, file), os.path.join(folders[0], file))
def sel_files(weeded_list):
	print("Choose the files you want to operate on:")
	print("(Multiple files can be separated by a space; 'a' for all files)")
	print("{:>3}{:>30}".format("0", "None"))
	for i,a in enumerate(weeded_list): print("{:>3}{:>30}".format(str(i+1),a))
	while True:
		option = input().lower().split()
		if all(b in ["a","0",*[str(a) for a in range(len(weeded_list)+1)]] for b in option): break
	if "0" in option: return False
	else: return weeded_list if "a" in option else [weeded_list[int(i)-1] for i in option]

def timeit(method):
	def timed(*args, **kw):
		ts = time.time()
		result = method(*args, **kw)
		te = time.time()
		print('{}:{:.2f} ms'.format(method.__name__, (te - ts) * 1000))
		return result
	return timed

class Var:
	conf_dir = os.path.dirname(__file__)
	conf_file = os.path.join(conf_dir, "user.txt")
	def __init__(self):
		self.heimdall_user = "heimdall_user"
		self.heimdall_mail = "heimdall_mail"
		self.heimdall_notification = False
		self.aguia_user = "aguia_user"
		self.athene_user = "athene_user"
		self.sub_s_name = "sub_s_name"
		self.heavy_atom = 36
		self.gjf_overwrite = False
		self.folder_op = True
		self.gauss_ext = ".com"
		self.comp_software = "gaussian"
		self.read_variables()

	def read_variables(self,conf_file=conf_file):
		if not os.path.isfile(conf_file): return
		with open(conf_file,mode="r") as file: options = file.readlines()
		options = [a.replace("=", " ").split() for a in options]
		options = [a for a in options if len(a) == 2]
		for a in options:
			if a[0] not in vars(self): continue
			elif a[0] == "heimdall_notification": setattr(self,a[0], a[1].lower() == "true")
			elif a[0] == "heavy_atom": setattr(self,a[0], int(a[1]) if a[1].isdigit() else 36)
			elif a[0] == "gjf_overwrite": setattr(self,a[0], a[1].lower() == "true")
			elif a[0] == "folder_op": setattr(self, a[0], a[1].lower() == "true")
			elif a[0] == "gauss_ext": setattr(self, a[0], a[1].lower() if a[1].lower() in [".gjf",".com"] else ".com")
			elif a[0] == "comp_software": setattr(self, a[0], a[1].lower() if a[1].lower() == "orca" else "gaussian")
			#elif a[0] == "menu_a":setattr(self,a[0]," ".join(a[1:]))
			else: setattr(self,a[0],a[1])
		return
	def set_variables(self):
		print("Please type 'chave' and press enter")
		if input().strip().lower() != "chave": return
		variables = [
			["To return", None],
			["HEIMDALL USER:", "heimdall_user"],
			["HEIMDALL MAIL:", "heimdall_mail"],
			["HEIMDALL SEND EMAIL (true/false):", "heimdall_notification"],
			["AGUIA USER:", "aguia_user"],
			["ATHENE USER:", "athene_user"],
			["SUBMISSION SCRIPT NAME:", "sub_s_name"],
			["ECP IS ADVISED FOR ELEMENTS LARGER THAN:", "heavy_atom"],
			["OVERWRITE .GJF FILES WITH NO PROMP (true/false):", "gjf_overwrite"],
			["AUTO OPERATE ON ALL FILES IN THE CWD (true/false):", "folder_op"],
			["GAUSSIAN INPUT FILE EXTENSION ('.gjf' or '.com'):", "gauss_ext"],
			["COMPUTATIONAL CHEMISTRY SOFTWARE (orca/gaussian):", "comp_software"]]
		try:
			import raapbs
			var = variables
		except ImportError:
			var = [variables[0]]
			var += [a for i, a in enumerate(variables) if i > 6]
		while True:
			print("Which variables do you want to set?")
			for i, v in enumerate(var):
				if v[1] is None: print("{} - {}".format(i, v[0]))
				else: print("{} - {} {}".format(i, v[0], getattr(self, v[1])))
			while True:
				option = input().strip()
				if option.isdigit():
					option = int(option)
					if option in range(len(var)): break
					else: print("Invalid Input!")
				else: print("Invalid Input!")
			if option == 0: break
			print("Enter variable '{}' value:".format(var[option][1]))
			while True:
				value = input()
				if len(value.split()) != 1:
					print("Variable name can neither be empty nor contain spaces!")
					continue
				elif var[option][1] in ["gjf_overwrite","folder_op","heimdall_notification"]:
					value = True if str(value).lower() in ["yes", "true"] else False
				elif var[option][1] == "heavy_atom":
					value = int(value) if value.isdigit() else self.heavy_atom
				elif var[option][1] == "comp_software":
					value = "orca" if value.lower() == "orca" else "gaussian"
				setattr(self,var[option][1],value)
				break
			self.write_save()
	def write_save(self,conf_file=conf_file):
		output = "\n".join(["{} = {}".format(a,getattr(self,a,)) for a in vars(self)])
		with open(conf_file, mode="w", newline="\n") as file: file.write(output)
		global preferences
		preferences = Var()
preferences = Var()

