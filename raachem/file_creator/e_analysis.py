import os, shutil
from raachem.file_class.log import LogFile
from raachem.util.gen_purp import file_weeder, read_item, preferences, timeit


def e_analysis(weeded_list):
	"""Analyses Gaussian log files in weeded_list and prints out a txt file with results in current directory"""
	file_name = "Extracted_data.txt"
	out = []
	out.append("{} File(s) were analyzed in total.".format(len(weeded_list)))
	out.append("{:^29}|{:^20}|{:^20}|{:^20}".format("Name", "Sum of Elet an TC", "Thermal Correct", "last SCF"))
	rel_e = []
	for i in weeded_list:
		log = LogFile(read_item(i))
		t = log.thermal
		try: last_SCF = log.scf_done[-1][-1]
		except: last_SCF = "None"
		if t == False: out.append("{:<30} no_thermo_data_found_Last_SCF: {:>20}".format(i, last_SCF)); continue
		out.append("{:<30}{:>20}{:>20}{:>20}".format(str(i), str(t[7]), str(t[3]), str(last_SCF)))
		rel_e.append([i, float(t[7]) * 627.5095, float(last_SCF) * 627.5095, float(t[3]) * 627.5095])
	if len(rel_e) > 1:
		print("Do you want the free energies to be reported relative to which item?")
		print("Enter 0 if you don't want to analyze them)")
		for idx, entry in enumerate(rel_e):
			print(" {:<4}{:<20}{:>25}".format(idx + 1, entry[0], round(entry[1], 2)))
		while True:
			option = input()
			if option in [str(a) for a in range(len(rel_e)+1)]: option = int(option); break
			else: print("Could notunderstand request")
		if option == 0: print("Leaving analysis\n")
		else:
			out.append("\nFree energies relative to {}, (Name,G,H,TC):".format(rel_e[option-1][0]))
			for i in rel_e:
				a = [i[0], *[round(i[n] - rel_e[option-1][n], 2) for n in [1,2,3]]]
				out.append("{:<26}{:>26}{:>10}{:>10}".format(*a))
			out.append("\n")
	with open(file_name,mode="w",newline="\n") as file: file.write("\n".join(out))
	print("\nDone! \nPlease lookup:\n\n" + os.path.join(os.getcwd(), file_name), "\n")

def rel_scf(list=False):
	energies = []
	for i in file_weeder([".log"]):
		log = LogFile(read_item(i))
		energies.append([log.name, log.scf_done[-1][1], log.normal_termin])
	if len(energies) == 0: print("No .log files in current folder"); return
	energies = [[i[0], float(i[1]), i[2]] for i in energies]
	energies.sort(key=lambda x: x[1])
	min_e = energies[0][1]
	energies = [[i[0], (i[1] - min_e) * 627.509, i[2]] for i in energies]
	if list == False: print("\n".join(["{:>30}{:>15f}{:>5}".format(*l) for l in energies]))
	elif list == True: return energies
@timeit
def csv_e_analysis():
	def evaluate_list(folder, logs=[]):
		for file in os.listdir(folder):
			if os.path.isdir(os.path.join(folder, file)): evaluate_list(os.path.join(folder, file), logs)
			elif file.endswith(".log"): logs.append(os.path.join(folder, file))
		return logs
	csv_list = []
	files = evaluate_list(os.getcwd())
	last = len(files)
	for i,a in enumerate(files):
		try:
			log = LogFile(read_item(os.path.relpath(a, os.getcwd())))
			# line = ["Termination","Free energy","Enthalphy","last_SCF","negativeFreq","TYP","File","Folder"]
			line = ["Yes" if log.normal_termin else "No",
					str(log.thermal[7]) if log.frequencies() else "No data",
					str(log.thermal[6]) if log.frequencies() else "No data",
					str(log.scf_done[-1][-1]) if log.normal_termin else "No data",
					str(len([a for a in log.frequencies() if float(a) < 0])) if log.frequencies() else "No data",
					log.calc_type,
					a,
					os.path.dirname(a),
					log.error_msg]

			if i+1 < last: print("\rEvaluating... {}/{}".format(i+1,last),end="")
			else: print("\rEvaluation done ({}/{}), saving '.csv' file...".format(i+1,last))
		except Exception as e:
			print("\nError on file:\n{}\n".format(a));print(e,"\n")
			line = ["-",
					"-",
					"-",
					"-",
					"-",
					"-",
					a,
					os.path.dirname(a),
					"-"]
		finally:
			csv_list.append(line)
	if not csv_list: return print("No .log files in {} directory".format(os.getcwd()))
	csv_list.sort(key=lambda x: x[0], reverse=True)
	csv_code = ["Free energy, +A, +B, +C, +D, -E, -F, Complex, Rel_E,-Freq ," +
				"TYP , Filename , INP , LOG , FOLD , Folder name, Done?, last_SCF, Hentalphy,Error msg"]
	for idx,line in enumerate(csv_list):
		row = [line[1],
			*["-" for _ in range(7)],
			"=(SUM($A{0}:$E{0})-SUM($F{0}:$G{0})-SUM($A${0}:$E${0})+SUM($F${0}:$G${0}))*627.509474".format(idx+2),
			line[4],
			line[5],
			os.path.basename(line[6]),
			'=HYPERLINK("{}";"Link")'.format(line[6]).replace(".log", preferences.gauss_ext) if os.path.isfile(line[6].replace(".log", preferences.gauss_ext)) else "-",
			'=HYPERLINK("{}";"Link")'.format(line[6]),
			'=HYPERLINK("{}";"Link")'.format(line[7]),
			os.path.relpath(line[7], os.getcwd()),
			line[0],
			line[3],
			line[2],
			line[8]]
		csv_code.append(",".join(row))
	try:
		with open("linked_analysis.csv", mode="w",newline="\n") as file: file.write("\n".join(csv_code))
		print("Done, please lookup:\n{}".format(os.path.join(os.getcwd(), "linked_analysis.csv")))
	except PermissionError:
		print("Error while saving file!\nIs the file '{}' already open?".format("linked_analysis.csv"))

def deduplicate():
	print("Analyzing energies...")
	energies = [[b.name,float(b.scf_done[-1][1]),b.normal_termin,b.last_xyz_obj()] for i in file_weeder([".log"]) for b in [LogFile(read_item(i))]]
	unique = energies
	if not unique: print("No log files to be analyzed"); return
	black_list, folder_mov = [], []
	print("Starting analysis...")
	for file in sorted(unique,key=lambda x: (x[2],-x[1]),reverse=True):
		if file[0] in black_list: continue
		black_list.append(file[0])
		sim_en = [i for i in unique if i[0] not in black_list and  i[1] + 1 > file[1] > i[1] - 1]
		if sim_en:
			duplicates = []
			for obj in sim_en:
				if obj[3].superimpose(file[3], ret="max_d",conv=6):
					print("{} is a duplicate of {}".format(obj[3].name(), file[3].name()))
					duplicates.append(obj[3].name())
					black_list.append(obj[3].name())
			if duplicates: folder_mov.append([file[3].name(),duplicates])
	for folder in folder_mov:
		subfolder = os.path.join(os.getcwd(),"duplicates_of_{}".format(folder[0].replace(".log","")))
		try: os.mkdir(subfolder)
		except FileExistsError:	pass
		print("Moving duplicates of {} to the following directory:\n{}".format(folder[0], subfolder))
		for file in folder[1]:
			for alt_file in file_weeder([file.replace(".log",".")]):
				try:
					shutil.move(os.path.join(os.getcwd(),alt_file),os.path.join(subfolder,alt_file))
					print("Moved: {}".format(alt_file))
				except PermissionError:
					print("Error while moving log files:")
					print("Maybe {} already exists in the following directory:\n{}".format(alt_file,subfolder))
				except FileNotFoundError:
					print("File {} not found in the following directory:\n{}".format(alt_file,subfolder))
	print("Done!")
