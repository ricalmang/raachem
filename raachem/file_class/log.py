from raachem.file_class.xyz import XyzFile
from raachem.util.constants import elements
from raachem.util.gen_purp import is_str_float
import functools

class LogFile:
	calc_types = ["TS","Red","IRC","Opt","SP"]
	def __init__(self,file_content):
		self.list = file_content
		self.lenght = len(self.list)
		self.name = self.list[0].strip()
		self.empty_line_idxs = []
		self.charge_mult = None
		self.input_geom_idx = None
		self.start_xyz_idxs = []
		self.end_resume_idxs = []
		self.start_resume_idxs = []
		self.linked_job_idxs = []
		self.multi_dash_idxs =[]
		self.scf_done = []
		####.thermal = ["ZPC","TCE","TCH","TCG","SZPE","STE","STH","STG"]
		self.thermal = [None , None, None, None, None , None, None, None]
		self.start_displ_block = None
		self.oc_orb_energies = []
		self.uno_orb_energies = []
		self.hash_line_idxs = []
		self.norm_term_idxs = []
		self.errors = []
		self.irc_points = []
		self.scan_points = []
		self.opt_points = []
		for i,a in enumerate(a.strip() for a in self.list):
			# i = index
			# a = line.strip()
			# b = line.split()
			# c = len(b)
			if a == "":                                                         self.empty_line_idxs.append(i); continue
			if a[-1] == "@":                                                    self.end_resume_idxs.append(i); continue
			elif a[0] == "1":
				if a.startswith(r"1\1"):                                      self.start_resume_idxs.append(i); continue
			if a[0].isdigit() or a[0].islower():                                                                continue
			elif a[0] == "-":
				if a.startswith("------"):                                      self.multi_dash_idxs.append(i); continue
			elif a[0] == "!":
				b = a.split(); c = len(b)
				if c == 4:
					condition_a = all(x in y for x,y in zip(b,("!",["Optimized","Non-Optimized"],"Parameters","!")))
					if condition_a:                                          self.scan_points.append([i,b[1]]); continue
			elif a[0] == "A":
				text_a = "Alpha  occ. eigenvalues --"
				text_b = "Alpha virt. eigenvalues --"
				if a.startswith(text_a):
					text = a.replace(text_a,"").split();                     self.oc_orb_energies.extend(text); continue
				elif a.startswith(text_b):
					text = a.replace(text_b,"").split();                    self.uno_orb_energies.extend(text); continue
			elif a[0] == "C":
				b = a.split(); c = len(b)
				if all((a.startswith("Charge"),self.charge_mult is None, c == 6)):
					pattern = ("Charge", "=", "Multiplicity", "=")
					if all(x == b[n] for x,n in zip(pattern,(0,1,3,4))):
						self.input_geom_idx = i;                                    self.charge_mult = b[2::3]; continue
			elif a[0] == "E":
				if a.startswith("Error"):                                                self.errors.append(i); continue
			elif a[0] == "F":
				condition_a = a.startswith("Frc consts  -- ") and self.start_displ_block is None
				if condition_a:                                                 self.start_displ_block = i - 4; continue
			elif a[0] == "I":
				if a == "Input orientation:":                                self.start_xyz_idxs.append(i + 5); continue
			elif a[0] == "L":
				if a.startswith("Link1:"):                                      self.linked_job_idxs.append(i); continue
			elif a[0] == "N":
				if a.startswith("Normal termination of Gaussian"):               self.norm_term_idxs.append(i); continue
			elif a[0] == "P":
				b = a.split(); c = len(b)
				if c != 6 or any(x != b[n] for x,n in zip(["Point","Number:","Path","Number:"],[0,1,3,4])):     continue
				if any(not b[n].isnumeric() for n in [2, 5]):                                                   continue
				else:                                                  self.irc_points.append([i, b[5], b[2]]); continue
			elif a[0] == "S":
				b = a.split(); c = len(b)
				if a == "Standard orientation:":                             self.start_xyz_idxs.append(i + 5); continue
				elif a.startswith("SCF Done:") and c > 5:                       self.scf_done.append([i,b[4]]); continue
				elif a.startswith("Sum of electronic and zero-point Energies="):       self.thermal[4] = b[-1]; continue
				elif a.startswith("Sum of electronic and thermal Energies="):          self.thermal[5] = b[-1]; continue
				elif a.startswith("Sum of electronic and thermal Enthalpies="):        self.thermal[6] = b[-1]; continue
				elif a.startswith("Sum of electronic and thermal Free Energies="):     self.thermal[7] = b[-1]; continue
				elif a.startswith("Step") and c == 9:
					x = ["Step", "number", "out", "of", "a", "maximum", "of"]
					y = [0, 1, 3, 4, 5, 6, 7]
					z = all(b[n].isnumeric() for n in [2, 8])
					if all(d == b[n] for d,n in zip(x,y)) and z:                     self.opt_points.append(i); continue
			elif a[0] == "T":
				b = a.split()
				if a.startswith("Thermal correction to Energy="):                      self.thermal[1] = b[-1]; continue
				elif a.startswith("Thermal correction to Enthalpy="):                  self.thermal[2] = b[-1]; continue
				elif a.startswith("Thermal correction to Gibbs Free Energy="):         self.thermal[3] = b[-1]; continue
			elif a[0] == "Z":
				b = a.split()
				if a.startswith("Zero-point correction="):                             self.thermal[0] = b[-2]; continue
			elif a[0] == "#":                                                    self.hash_line_idxs.append(i); continue
		#--------------------------------------------POST PROCESSING----------------------------------------------------
		x = None if self.start_xyz_idxs is None else [min(a for a in self.multi_dash_idxs if a > b) for b in self.start_xyz_idxs]
		self.end_xyz_idxs = x
		self.scan_end = [min(a for a in self.multi_dash_idxs if a > b[0]) for b in self.scan_points]
		x = None if self.start_displ_block is None else min(a for a in self.empty_line_idxs if a > self.start_displ_block)
		self.end_displ_block = x
		self.displ_block = None if self.end_displ_block is None else self.list[self.start_displ_block:self.end_displ_block]
		try:
			x = None if self.hash_line_idxs is None else min(a for a in self.multi_dash_idxs if a > self.hash_line_idxs[0])
			x = None if self.hash_line_idxs is None else "".join([a.lstrip() for a in self.list[self.hash_line_idxs[0]:x]])
			self.raw_route = x
		except IndexError as e:
			print("Error while finding route section of log file")
			print(e)
			print(self.name)
		try:
			self.oc_orb_energies = [float(b) for b in self.oc_orb_energies]
			self.uno_orb_energies = [float(a) for a in self.uno_orb_energies]
			self.homo = max(self.oc_orb_energies) if self.oc_orb_energies else None
			self.lumo = min(self.uno_orb_energies) if self.uno_orb_energies else None
			self.homo_lumo_gap = self.lumo - self.homo if not any([self.lumo is None, self.homo is None]) else None
		except ValueError as e:
			print("Error finding homo and lumo energies")
			print(e)
			print(self.name)
			self.homo, self.lumo, self.homo_lumo_gap = None, None, None
		if all([self.start_resume_idxs,self.end_resume_idxs]):
			x = ["".join([x.strip() for x in self.list[a:b]]).split("\\") for a,b in zip(self.start_resume_idxs,self.end_resume_idxs)]
			self.resumes = x
		else: self.resumes = None
		# --------------------------------------------------ASSURANCE---------------------------------------------------
		self.init_errors = []
		if self.charge_mult is None:
			self.init_errors.append("Charge and multiplicity could not be identified!")
		if len(self.start_resume_idxs) != len(self.end_resume_idxs):
			self.init_errors.append("Inconsistent resumes")
		if len(self.name.split()) != 1:
			self.init_errors.append("Name must not contain empty spaces or be empty")
		if not self.list[1].strip().startswith("Entering Gaussian System"):
			self.init_errors.append("Is this a Gaussian log file?")
		if not self.start_xyz_idxs is None:
			if len(self.start_xyz_idxs) != len(self.end_xyz_idxs):
				self.init_errors.append("Found an inconsistent number of geometries")
		if not any([self.homo is None, self.lumo is None]):
			if self.homo > self.lumo:
				self.init_errors.append("Lumo is lower than homo?")
		if self.init_errors:
			for a in self.init_errors: print(a)
			print("Errors above were found on file\n{}".format(self.name))
		#self.loghelp()
		#print(self.raw_route)
		#print(self.first_xyz_obj())
	@functools.lru_cache(maxsize=1)
	def loghelp(self):
		for a in vars(self):
			if a != "list":
				print(a.upper(),"--->",getattr(self,a))
	@functools.lru_cache(maxsize=1)
	def xyz_cord_block(self,start_idx,end_idx):
		data = [a.split() for a in self.list[start_idx:end_idx]]
		return [[elements[int(l[1])],*[l[i] for i in [3,4,5]]] for l in data]
	@functools.lru_cache(maxsize=1)
	def last_cord_block(self):
		if not all([self.xyz_cord_block, self.end_xyz_idxs]):
			if self.resumes:
				print("WARNING: Coordinates will be drawn from the last job abstract:")
				print("lines {} - {} of file:".format(self.start_resume_idxs[-1],self.end_resume_idxs[-1]))
				print("{}".format(self.name))
				return LogAbstract(self.resumes[0]).xyz_object().cord_block()
			else: return None
		else:
			return self.xyz_cord_block(self.start_xyz_idxs[-1],self.end_xyz_idxs[-1])
	@functools.lru_cache(maxsize=1)
	def first_cord_block(self):
		if not all([self.start_xyz_idxs,self.end_xyz_idxs]):
			if self.input_geom_idx:
				coordinates = []
				for i,a in enumerate(self.list[self.input_geom_idx:]):
					if i > 5 and not coordinates: break
					a = a.split()
					if len(a) == 4:
						if a[0] in elements and all(is_str_float(a[n]) for n in [1, 2, 3]):
							coordinates.append(a)
						elif coordinates: break
					elif coordinates: break
				return coordinates
			else: return None
		else:
			return self.xyz_cord_block(self.start_xyz_idxs[0],self.end_xyz_idxs[0])
	@functools.lru_cache(maxsize=1)
	def _n_atoms(self):
		if self.last_cord_block():
			return len(self.last_cord_block())
		elif self.first_cord_block():
			return len(self.first_cord_block())
	def any_xyz_obj(self,a_idx,b_idx,title=" ",name=False):
		if name == False: name = self.name
		return XyzFile([name, self.n_atoms, title, *(" ".join(l) for l in self.xyz_cord_block(a_idx,b_idx))])
	@functools.lru_cache(maxsize=1)
	def last_xyz_obj(self):
		return XyzFile([self.name,self.n_atoms," ",*(" ".join(l) for l in self.last_cord_block())])
	@functools.lru_cache(maxsize=1)
	def first_xyz_obj(self):
		return XyzFile([self.name,self.n_atoms," ",*(" ".join(l) for l in self.first_cord_block())])
	@functools.lru_cache(maxsize=1)
	def low_e_xyz_obj(self):
		if self.calc_type == "SP": return None
		else:
			xyzs = {"TS":self.opt,"Red":self.scan_geoms,"IRC":self.irc,"Opt":self.opt}[self.calc_type]()
			if len(xyzs) == 0: return None
			else: return sorted(xyzs,key= lambda x: float(x.title()) if is_str_float(x.title()) else 1)[0]
	@functools.lru_cache(maxsize=1)
	def frequencies(self):
		if self.displ_block is None: return False
		else: return [j for a in [i.split()[2:] for i in self.displ_block[2::self.n_atoms+7]] for j in a]
	@functools.lru_cache(maxsize=1)
	def displ_for_freq_idx(self,freq_idx):
		if self.displ_block is None: return []
		displ = []
		for num in range(self.n_atoms):
			displ.append([j for a in [i.split()[2:] for i in self.displ_block[7+num::self.n_atoms+7]] for j in a])
		displ_for_freq_str = [a[freq_idx*3:freq_idx*3+3] for a in displ]
		displ_for_freq_float = [[float(i) for i in b] for b in displ_for_freq_str]
		return displ_for_freq_float
	@functools.lru_cache(maxsize=1)
	def _calc_type(self):
		if self.raw_route:
			r_sect = [self.raw_route]
			for x in [None, "/", "(", ")", ",", "=", "%", ":"]:
				r_sect = [a for b in [i.split(x) for i in r_sect] for a in b if len(a) > 1]
			r_sect = [a.lower() for a in r_sect]
			if "ts" in r_sect: return "TS"
			elif any(True for a in r_sect if a in ("modredundant", "readoptimize", "readfreeze")): return "Red"
			elif "irc" in r_sect: return "IRC"
			elif "opt" in r_sect: return "Opt"
			else: return "SP"
		else: return "No data"
	@functools.lru_cache(maxsize=1)
	def _normal_termin(self):
		return any(True if "Normal termination of Gaussian" in l else False for l in self.list[-5:])
	def _error_msg(self):
		error_idxs = [a for a in self.errors if a + 5 > self.lenght]
		if error_idxs: return " | ".join([self.list[n] for n in error_idxs])
		else: return "No data"
	@functools.lru_cache(maxsize=1)
	def irc(self):
		if not all([self.start_xyz_idxs,self.end_xyz_idxs,self.irc_points,self.scf_done]): return []
		points = self.irc_points
		scf = [max(self.scf_done,key=lambda x: x[0] if x[0] < a[0] else 0)[1] for a in points]
		a_idx = [max(self.start_xyz_idxs,key=lambda x: x if x < a[0] else 0) for a in points]
		b_idx = [max(self.end_xyz_idxs,key=lambda x: x if x < a[0] else 0) for a in points]
		points = [[*d[1:],c,self.any_xyz_obj(a,b,title=c)] for a,b,c,d in zip(a_idx,b_idx,scf,points)]
		path_a = sorted([a for a in points if a[0] == "1"], key = lambda x: int(x[1]), reverse=True)
		path_b = [a for a in points if a[0] == "2"]
		return [a[3] for a in [*path_a,*path_b]]
	@functools.lru_cache(maxsize=1)
	def opt(self):
		if not all([self.start_xyz_idxs,self.end_xyz_idxs,self.opt_points,self.scf_done]): return []
		points = self.opt_points
		scf = [max(self.scf_done,key=lambda x: x[0] if x[0] < a else 0)[1] for a in points]
		a_idx = [max(self.start_xyz_idxs,key=lambda x: x if x < a else 0) for a in points]
		b_idx = [max(self.end_xyz_idxs,key=lambda x: x if x < a else 0) for a in points]
		return [self.any_xyz_obj(a,b,title=c) for a,b,c in zip(a_idx,b_idx,scf)]
	@functools.lru_cache(maxsize=1)
	def scan_geoms(self):
		if not all([self.start_xyz_idxs, self.end_xyz_idxs, self.scan_points, self.scf_done]): return []
		geoms = []
		points = self.scan_points
		start_idx = [min(i for i in self.start_xyz_idxs if i > b[0]) for b in points]
		end_idx = [min(i for i in self.end_xyz_idxs if i > b[0]) for b in points]
		scf_idx = [max(i for i in self.scf_done if i[0] < b[0]) for b in points]
		for i,(a,b,c,d) in enumerate(zip(start_idx,end_idx,scf_idx,self.scan_points)):
			name = self.name.replace(".log","_" + str(i+1)+".xyz")
			if d[1] == "Optimized": print("Optimized geometry found at line {}!".format(d[0]))
			elif d[1] == "Non-Optimized": print("Non-Optimized1 geometry found at line {}!".format(d[0]))
			geoms.append(self.any_xyz_obj(a,b,title=str(c[1]), name=name))
		if len(geoms) == 0:
			print("No Optimized geometries found for {} file".format(self.name()))
		return geoms

	n_atoms = property(_n_atoms)
	normal_termin = property(_normal_termin)
	calc_type = property(_calc_type)
	error_msg = property(_error_msg)

class LogAbstract:
	def __init__(self,content):
		self.list = content
		self.version = None
		self.dipole = None
		self.img_freq = None
		self.hash_line = None
		for i,a in enumerate(self.list):
			a = a.lstrip()
			if a.lstrip == "": print("Empty!");continue
			elif a.startswith("Version="): self.version = a.replace("Version=","")
			elif a.startswith("#"): self.hash_line = i
			elif a.startswith("NImag="): self.img_freq = a.replace("NImag=0","")
			elif a.startswith("DipoleDeriv="): self.img_freq = a.replace("DipoleDeriv=","")
			else: continue
	def __str__(self):
		return "\n".join(self.list)
	def read_strucure(self):
		charge_mult = None
		title = None
		coordinates = []
		for i,a in enumerate(self.list[self.hash_line:]):
			if i > 5 and not coordinates: break
			a = a.split(",")
			if len(a) == 2 and not coordinates:	charge_mult = a; continue
			if len(a) == 4:
				if a[0] in elements and all(is_str_float(a[n]) for n in [1,2,3]):
					coordinates.append("   ".join(a))
				elif coordinates: break
			elif coordinates: break
		return charge_mult, XyzFile([self.list[0],str(len(coordinates)),title,*coordinates])
	def charge_mult(self):
		return self.read_strucure()[0]
	def xyz_object(self):
		return self.read_strucure()[1]

