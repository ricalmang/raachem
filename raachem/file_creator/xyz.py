import os, random, math
from raachem.util.gen_purp import read_item, file_weeder, preferences
from raachem.file_class.log import LogFile
from raachem.file_class.xyz import XyzFile
from raachem.file_class.inp import InpFile
from raachem.file_class.gjf import GjfFile

def xyz_ent():
	for i in file_weeder([".xyz"]):
		xyz = XyzFile(read_item(i)).enantiomer()
		xyz.save_file()
def input_to_xyz():
	extension = preferences.gauss_ext if preferences.comp_software == "gaussian" else ".inp"
	for i in file_weeder([extension]):
		if preferences.comp_software == "orca":
			xyz = InpFile(read_item(i)).xyz_obj()
			xyz.save_file()
		elif preferences.comp_software == "gaussian":
			xyz = GjfFile(read_item(i)).xyz_obj()
			xyz.save_file()
def log_to_xyz(weeded_list):
	while True:
		print("Which geometry do you want?")
		print("0 - Cancel")
		print("1 - Last")
		print("2 - Lowest energy")
		print("3 - First")
		geom = input()
		if geom == "0": return
		if geom in ["1","2","3"]: break
		else: print("Invalid option!")
	for i in weeded_list:
		try:
			if geom == "1":	LogFile(read_item(i)).last_xyz_obj().save_file()
			elif geom == "2": LogFile(read_item(i)).low_e_xyz_obj().save_file()
			elif geom == "3": LogFile(read_item(i)).first_xyz_obj().save_file()
		except Exception as e:
			print("Error on file: {}".format(i))
			print(e)

def log_to_xyz_scan(weeded_list):
	print("\nWhich file(s) do you want to analyse? (Separate multiple entries by space | Enter 0 to give up)\n")
	for idx,file in enumerate(weeded_list):
		print("{:<5}{:<40}".format(idx+1,file))
	while True:
		option=input().split()
		if all(True if a in [str(b) for b in range(len(weeded_list)+1)] else False for a in option):
			if option: break
	if "0" in option: print("No files will be analyzed!");return
	else:
		for item in [int(i)-1 for i in option]:
			log = LogFile(read_item(weeded_list[item]))
			if log.calc_type.lower() == "red": xyzs = log.scan_geoms(); ext_name = "_scan_traject.xyz"
			elif log.calc_type.lower() == "irc": xyzs = log.irc(); ext_name = "_irc_traject.xyz"
			elif log.calc_type.lower() in ["opt","ts"]: xyzs = log.opt(); ext_name = "_opt_traject.xyz"
			else: print("The file is not a scan, opt or irc log file"); continue
			if len(xyzs) < 2: print("The file contains less than 2 steps"); return
			max_e = max(float(i.title()) for i in xyzs)
			min_e = min(float(i.title()) for i in xyzs)
			c=[0,"percentage","kcal"]
			print("{:5}{:^90}{:>8}".format("Entry"," >---> Relative Energy >---> "," kcal"))
			for entry in [float(i.title()) for i in xyzs]:
				c[0]=c[0]+1
				c[1]= "|"*int(90*((float(entry)-min_e)/(max_e-min_e)) if max_e != min_e else 0)
				c[2]= (float(entry)-min_e)*627.5
				print("{:<5}{:<90}".format(c[0],c[1]),"{:>6.5f}".format(c[2]) if ext_name == "_opt_traject.xyz" else "{:>6.2f}".format(c[2]))
			file_name = weeded_list[item].replace(".log",ext_name)
			with open(file_name,mode="w",newline="\n") as file: file.write("\n".join([a for i in xyzs for a in i.return_print()]))

def log_freq_xyz(weeded_list):
	for i in weeded_list:
		log = LogFile(read_item(i))
		if not log.frequencies():
			print("No frequencies found in file:{}!".format(i))
			continue
		option = None
		while True:
			print("Analizing file: {}".format(i))
			print("Chose a frequency to split (0-To cancel/m-For more)")
			for idx,item in enumerate(log.frequencies()):
				print("{:>5}:{:>15}".format(idx+1,item))
				if option == "m": continue
				if idx > 3: break
			option=input()
			if option == "m": continue
			elif option in [str(a+1) for a in range(len(log.frequencies()))]: break
			elif option == "0": break
			else: print("Invalid option, try again!")
		if option == "0": continue
		print("Give a multiplier (recomended multiplier = 1) ")
		while True:
			try:
				mult=float(input().replace(",","."))
				break
			except ValueError:
				print("Not a valid number")
		left = log.last_xyz_obj().displace(-mult,log.displ_for_freq_idx(int(option)-1))
		left.list[0] = i.replace(".log","_l.xyz")
		left.save_file()
		right = log.last_xyz_obj().displace(mult,log.displ_for_freq_idx(int(option)-1))
		right.list[0] = i.replace(".log","_r.xyz")
		right.save_file()
		print("\nJob Done!\n")
	print("Finished !")

def superimpose_alg():
	while True:
		print("Chosse an option:")
		print("0 - To cancel")
		print("1 - For batch mode")
		print("2 - For single structure(In screen printing)")
		mode = input()
		if mode in [str(n) for n in range(3)]: break
	if mode == "0":
		return
	elif mode == "1":
		if os.path.exists(os.path.join(os.getcwd(), "rotated")):
			print("Rotated directory already exists in current directory!")
			print("Please remove it and try again.")
			return
		os.mkdir(os.path.join(os.getcwd(), "rotated"))

		xyz_1 = read_item(None, "Which item do you wish to compare to?", [".xyz"])
		if xyz_1: xyz_1 = XyzFile(xyz_1)
		else:return

		while True:
			num_atoms = input("compare first N atoms (0 for all):\n")
			if num_atoms in [str(a) for a in range(xyz_1.n_atoms())]:
				num_atoms = int(num_atoms)
				break
			else:
				print("Invalid input for {}".format(xyz_1.name()))
		for file in file_weeder(["xyz"]):
			if not file == xyz_1.name():
				xyz_2 = XyzFile(read_item(file))
				xyz_3 = xyz_2.superimpose(xyz_1, num_atoms, True)
				xyz_3.save_file(os.path.join(os.getcwd(), "rotated"))
		xyz_1.std_cord(num_atoms).save_file(os.path.join(os.getcwd(), "rotated"))
	elif mode == "2":
		xyz_1 = read_item(None, "Which item do you wish to compare to?", [".xyz"])
		if xyz_1: xyz_1 = XyzFile(xyz_1)
		else: return

		xyz_2 = read_item(None, "Which item do you wish to rotate?", [".xyz"])
		if xyz_2: xyz_2 = XyzFile(xyz_2)
		else:return

		while True:
			num_atoms = input("compare first N atoms (0 for all):\n")
			if num_atoms in [str(a) for a in range(min([xyz_1.n_atoms(),xyz_2.n_atoms()]))]:
				num_atoms = int(num_atoms)
				break
			else:
				print("Invalid input for either or both {} and {}".format(xyz_1.name(),xyz_2.name()))
		xyz_2.superimpose(xyz_1, num_atoms, True).print_file()
		xyz_1.std_cord(num_atoms).print_file()
def geodes_int():
	"""Requests input for solvation algorithm"""
	def rot_alg(xyz_1, xyz_2, angles, p_angles, geodes, save):
		"""Rotate xyz_2 by angles and p_angles then aproximate it to xyz_1 on geodes trajectories"""
		trajectory = []
		trajectory_mult = []
		xyz_old = xyz_2
		iteration = 0
		for axis in p_angles:
			xyz_2 = xyz_old.rotate(axis[1], axis[0])
			for angle in angles:
				xyz_2 = xyz_2.rotate(angle, "y")
				trajectory_other = []
				for dot in geodes:
					xyz_3 = aproximate(xyz_1, xyz_2, dot)
					xyz_4 = xyz_1 + xyz_3
					trajectory.append(xyz_4.return_print())
					trajectory_other.append(xyz_3.return_print())
					iteration += 1
					if save == "y":
						xyz_4.list[0] = "{}_{}.xyz".format(os.path.splitext(xyz_4.name())[0], str(iteration))
						xyz_4.save_file()
				trajectory_other.append(xyz_1.return_print())
				cord = [j for i in trajectory_other for j in i[2:]]
				trajectory_mult.append(str(len(cord)) + "\n\n")
				for line in cord:
					trajectory_mult.append(str(line) + "\n")
		# save suporposed trajectory
		file_path = os.path.join(os.getcwd(), "superposed.xyz")
		with open(file_path, "w") as file:
			for geom in trajectory_mult:
				file.write(str(geom))
		file_path = os.path.join(os.getcwd(), "trajectory.xyz")
		# save parwise trajectory
		with open(file_path, "w") as file:
			for xyz in trajectory:
				for line in xyz:
					file.write(str(line) + "\n")

	def aproximate(xyz_1, xyz_2, vertex):
		"""Aproximate xyz_2 to xyz_1 from trajectory vertex. Returns [xyz_1+xyz_2, xyz_2]"""

		def min_dist(a, b):
			return min(math.sqrt(sum(pow(c - d, 2) for c, d in zip(bb, aa))) for bb in b for aa in a)

		xyz_a = [[float(a) for a in i] for i in xyz_1.cord_strip()]
		xyz_b = [[float(a) + b * 50 for a, b in zip(i, vertex)] for i in xyz_2.cord_strip()]
		for num in range(5):
			step_size = 1 / 2 ** num
			while True:
				if min_dist(xyz_a, xyz_b) < 2:
					xyz_b = [[a + b * step_size for a, b in zip(i, vertex)] for i in xyz_b]
					break
				else:
					xyz_b = [[a - b * step_size for a, b in zip(i, vertex)] for i in xyz_b]
		cordinates = [" ".join([a, *[str(c) for c in b]]) for a, b in zip(xyz_2.all_elements(), xyz_b)]
		return XyzFile([xyz_2.name(), xyz_2.n_atoms(), " ", *cordinates])
	def start():
		angles, p_angles, geodes = None, None, None
		xyz_1 = read_item(None, "Choose the molecule you want to solvate?", ["xyz"])
		if xyz_1: xyz_1 = XyzFile(xyz_1).std_cord()
		else: return

		while True:
			print("Choose a trajectory englobing method")
			print("0 - To cancel")
			print("1 - Cubic gobling")
			print("2 - Icosahedron gobling")
			print("3 - Give a fulerene")
			usr_opt = input()
			if usr_opt == "0": return
			if usr_opt in ["1","2","3"]: break
		if usr_opt == "1":
			geodes = [[a,b,c] for a in [-1,1] for b in [-1,1] for c in [-1,1]]
		elif usr_opt == "2":
			k = (1 + math.sqrt(5)) / 2
			geodes = [[0, 1, k], [0, -1, k], [0, -1, -k],
					[0, 1, -k], [k, 0, 1], [-k, 0, 1],
					[-k, 0, -1],[k, 0, -1], [1, k, 0],
					[-1, k, 0], [-1, -k, 0], [1, -k, 0]]
		elif usr_opt == "3":
			geodes = read_item(None, "Chose a fulerene", ["xyz"])
			if geodes:
				geodes = XyzFile(geodes).std_cord().cord_strip()
				geodes = [[float(n) for n in i] for i in geodes]
			else:
				print("You need a to choose a trajectory engobling method")
				return

		xyz_2 = read_item(None, "Chose the solvent or counterion", ["xyz"])
		if xyz_2: xyz_2 = XyzFile(xyz_2).std_cord()
		else: return

		while True:
			print("Orientation algorithm:")
			print("0 - To cancel")
			print("1 - Espherical ion")
			print("2 - 12 Random cubic orientations")
			print("3 - 24 Sistematic cubic rotations")
			usr_opt = input()
			if usr_opt in [str(a) for a in range(4)]: break
		if usr_opt == "0":
			return
		elif usr_opt == "1":
			angles = [0]
			p_angles = [["x", 0]]
		elif usr_opt == "2":
			angles = random.choices([math.pi * i for i in [0, 0.5, 1, 1.5]], k=3)
			rot = [["x", 0], ["x", .5], ["x", 1], ["x", 1.5], ["z", 0.5], ["z", 1.5]]
			p_angles = random.choices([[i[0], math.pi * i[1]] for i in rot], k=4)
		elif usr_opt == "3":
			angles = [math.pi * i for i in [0, 0.5, 1, 1.5]]
			rot = [["x", 0], ["x", .5], ["x", 1], ["x", 1.5], ["z", 0.5], ["z", 1.5]]
			p_angles = [[i[0], math.pi * i[1]] for i in rot]

		while True:
			save = input("Save individual files? (y/n)")
			if save in ["y", "n"]: break
		rot_alg(xyz_1, xyz_2, angles, p_angles, geodes, save)
	start()

