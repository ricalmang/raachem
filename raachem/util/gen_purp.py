import os, shutil, time
cf=os.getcwd()
conf_dir = os.path.dirname(__file__)
conf_file = os.path.join(conf_dir, "user.txt")

def read_item(file_name=None, promp=False, extension=None, cf=cf):
	"""Reads an xyz, gjf or log item and returns a list of its contents ready for class instantiation"""
	if extension is None:
		extension = [".xyz"]
	if promp != False:
		print(promp)
		files = file_weeder(extension)
		if len(files) == 0:
			print("Sorry no such files in current directory!")
			return []
		for idx,item in enumerate(files):
			print("{:<5}{:<25}".format(idx+1,item))
		while True:
			try:
				file_name = files[int(input())-1]
				break
			except: print("Invalid input!")
	with open(os.path.join(cf,file_name),"r") as in_file:
		in_content = in_file.read().splitlines()
	in_content.insert(0,file_name)
	return in_content
def file_weeder(ext_to_weed,cf=cf, promp=True):
	"""Looks up files with the extensions provided in current directory"""
	fulllist=os.listdir(cf)
	weeded_list=[]
	if type(ext_to_weed) == str:
		ext_to_weed = [ext_to_weed]
	for extension in ext_to_weed:
		matching = [s for s in fulllist if extension in s]
		for match in matching:
			weeded_list.append(match)
	if promp:
		if len(weeded_list) == 0:
			print("No {} files found in current directory!".format("or".join(ext_to_weed)))
			print(cf)
			return []
	return sorted(weeded_list)
def w_any(lines,write_mod="a",filename="alg_log.txt",folder=cf):
	"""Generally wites a file"""
	with open((os.path.join(folder,filename)),write_mod,newline= "\n") as activity_log:
		for line in lines:
			activity_log.write(str(line))
			activity_log.write("\n")
def is_str_float(i):
	"""Check if a string can be converted into a float"""
	try:
		float(i)
		return True
	except ValueError:
		return False

def mv_up_folder():
	"""Move files with a chossen extension up a folder"""
	while True:
		extensions = [[0, "return"], ["1", ".log"], ["2", ".gjf"], ["3", ".xyz"]]
		print("Select a file extension to move up a folder")
		for idx, ext in enumerate(extensions):
			if idx == 0:
				print("0   - Cancel request")
			else:
				print("{:<3} - {:>10}".format(ext[0], ext[1]))
		extension = {str(i[0]): i[1] for i in extensions}.get(input(), None)
		print(extension)
		if extension == "return": return
		if extension != None: break
		print("Invalid Option")
	folders = [x[0] for x in os.walk(cf)]
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
	for i,a in enumerate(weeded_list):
		print("{:>3}{:>30}".format(str(i+1),a))
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
		print('{}:{} ms'.format(method.__name__, (te - ts) * 1000))
		return result
	return timed

class Var:
	conf_dir = os.path.dirname(__file__)
	conf_file = os.path.join(conf_dir, "user.txt")
	def __init__(self):
		self.heimdall_user = "heimdall_user"
		self.heimdall_mail = "heimdall_mail"
		self.heimdall_notification = "heimdall_notification"
		self.aguia_user = "aguia_user"
		self.athene_user = "athene_user"
		self.sub_s_name = "sub_s_name"
		self.heavy_atom = 36
		self.gjf_overwrite = False
		self.folder_op = True
		self.read_variables()
	def read_variables(self):
		if not os.path.isfile(conf_file): return
		with open(conf_file) as file:
			options = file.readlines()
		options = [a.replace("=", " ").split() for a in options]
		options = [a for a in options if len(a) == 2]
		for a in options:
			if a[0] not in vars(self): continue
			elif a[0] == "heimdall_notification": setattr(self,a[0], a[1].lower() == "true")
			elif a[0] == "heavy_e": setattr(self,a[0], int(a[1]) if a[1].isdigit() else 36)
			elif a[0] == "gjf_overwrite": setattr(self,a[0], a[1].lower() == "true")
			elif a[0] == "folder_op": setattr(self, a[0], a[1].lower() == "true")
			else: setattr(self,a[0],a[1])
		return
	def set_variables(self):
		print("Please type 'chave' and press enter")
		if input().strip().lower() != "chave": return
		while True:
			print("Which variables do you want to set?")
			print("0 - To return")
			print("1 - HEIMDALL USER: {}".format(self.heimdall_user))
			print("2 - HEIMDALL MAIL: {}".format(self.heimdall_mail))
			print("3 - HEIMDALL SEND EMAIL?: {}".format("Yes" if self.heimdall_notification else "No"))
			print("4 - AGUIA USER: {}".format(self.aguia_user))
			print("5 - ATHENE USER: {}".format(self.athene_user))
			print("6 - SUBMISSION SCRIPT NAME: {}".format(self.sub_s_name))
			print("7 - ECP IS ADVISED FOR ELEMENTS LARGER THAN: {}".format(self.heavy_atom))
			print("8 - OVERWRITE .GJF FILES WITH NO PROMP: {}".format("Yes" if self.gjf_overwrite else "No"))
			print("9 - AUTO OPERATE ON ALL FILES IN THE CWD: {}".format("Yes" if self.folder_op else "No"))
			variables = {"0":None,"1":"heimdall_user","2":"heimdall_mail","3":"heimdall_notification","4":"aguia_user",
						"5":"athene_user","6":"sub_s_name","7":"heavy_atom","8":"gjf_overwrite","9":"folder_op"}
			while True:
				option = input().strip()
				if option in variables:	break
				else: print("Invalid Input!")
			if option == "0": break
			print("Enter variable '{}' value:".format(variables[option]))
			while True:
				value = input()
				if len(value.split()) != 1:
					print("Variable name can neither be empty nor contain spaces!")
					continue
				elif variables[option] in ["gjf_overwrite","folder_op"]:
					value = "true" if str(value).lower() in ["yes", "true"] else "false"
				elif variables[option] == "heavy_atom":
					value = value if value.isdigit() else str(self.heavy_atom)
				setattr(self,variables[option],value)
				break
			self.write_save()
	def write_save(self):
		otput = []
		otput.append("#HEIMDALL VARIABLES")
		otput.append("heimdall_user = {}".format(self.heimdall_user))
		otput.append("heimdall_mail = {}".format(self.heimdall_mail))
		otput.append("heimdall_notification = {}".format(self.heimdall_notification))
		otput.append("#AGUIA VARIABLES")
		otput.append("aguia_user = {}".format(self.aguia_user))
		otput.append("#ATHENE VARIABLES")
		otput.append("athene_user = {}".format(self.athene_user))
		otput.append("#SUBMISSION SCRIPT NAME")
		otput.append("sub_s_name = {}".format(self.sub_s_name))
		otput.append("#ATOMS WITH RECOMMENDED ECP ARE HEAVIER THAN:")
		otput.append("heavy_e = {}".format(self.heavy_atom))
		otput.append("#SHOULD .GJF FILES BE OVERWRITTEN WITH NO PROMP?")
		otput.append("gjf_overwrite = {}".format(self.gjf_overwrite))
		otput.append("#OPERATE ON ALL FILES IN FOLDER?")
		otput.append("folder_op = {}".format(self.folder_op))
		w_any(otput, write_mod="w", filename="user.txt", folder=conf_dir)
