import os, math, functools
import numpy as np
from raachem.util.constants import elements, element_radii
from raachem.util.gen_purp import is_str_float

class XyzFile:
	def __init__(self,file_content):
		self.list = file_content

		if len(self.list) < 2: raise Exception(".xyz Object is empty?")
		elif not (str(self.list[1]).strip().isdigit() and len(str(self.list[1]).split()) == 1):
			print("{} is not a proper .xyz file\nAttempting to read it anyway!".format(self.list[0]))
			try_xyz = []
			for line in self.list:
				line = line.split()
				if len(line) != 4: continue
				if not all(is_str_float(line[i]) for i in range(1, 4)): continue
				if line[0] in elements[0:]:
					try_xyz.append(" ".join(line))
					continue
				try:
					line[0] = elements[int(line[0])]
					try_xyz.append(" ".join(line))
				except:
					raise Exception("Could not understand file {}".format(self.list[0]))
			try_xyz.insert(0,len(try_xyz))
			try_xyz.insert(1," ")
			try_xyz.insert(0,self.list[0])
			self.list = try_xyz
		self.list_l = [str(a).split() for a in self.list]
	def __add__(self,other):
		assert type(self) == type(other), "Operation '+' allowed only for two XYZ objects"
		new = [os.path.splitext(self.name())[0]+"_"+other.name(), str(self.n_atoms()+other.n_atoms()),
			   self.title()+" "+other.title(),*(self.form_cord_block() + other.form_cord_block())]
		return XyzFile(new)
	def __sub__(self, other):
		el_a = self.all_elements()
		el_b = other.all_elements()
		assert len(el_a) > len (el_b), "Can't subtract a larger structure from a smaller one"
		assert type(self) == type(other), "Operation '-' allowed only for two XYZ objects"
		idxs_to_rem = []
		for n in range(len(el_a) - len(el_b)):
			if all([True if el_a[n+i] == a else False for i,a in enumerate(el_b)]):
				idxs_to_rem = [c+n for c in range(len(el_b))]
				break
		if len(idxs_to_rem) ==  0: print("Could not subtract value!")
		xyz_cord = [a for idx,a in enumerate(self.form_cord_block()) if idx not in idxs_to_rem]
		new = [os.path.splitext(self.name())[0]+"-"+other.name(), str(self.n_atoms()-other.n_atoms()),
			   self.title()+"-"+other.title(),*xyz_cord]
		return XyzFile(new)
	def __str__(self):
		return "\n".join(self.return_print())
	@functools.lru_cache(maxsize=1)
	def name(self):
		if len(self.list[0]) == 0: raise Exception(".xyz Object has no name")
		return self.list[0]
	@functools.lru_cache(maxsize=1)
	def n_atoms(self):
		if any([len(str(self.list[1]).split()) != 1, not str(self.list[1]).isnumeric()]):
			raise Exception("First line of {} (.xyz type) file should contain only the number of atoms in the geometry!".format(self.name()))
		return int(self.list[1])
	@functools.lru_cache(maxsize=1)
	def title(self):
		return self.list[2]
	@functools.lru_cache(maxsize=1)
	def cord_block(self):
		cordinates = []
		for idx,line in enumerate(self.list_l):
			if idx <= 2: continue
			if idx >= self.n_atoms() + 3: continue
			if line[0] in elements:	cordinates.append(line)
			else: cordinates.append([elements[int(line[0])],*line[0:]])
		return cordinates
	@functools.lru_cache(maxsize=1)
	def form_cord_block(self):
		return ["{:<5}{:>20f}{:>20f}{:>20f}".format(x[0], *[round(float(x[a]), 6) for a in [1, 2, 3]]) for x in self.cord_block()]
	@functools.lru_cache(maxsize=1)
	def cord_strip(self):
		return [line[1:] for line in self.cord_block()]
	@functools.lru_cache(maxsize=1)
	def all_elements(self):
		return [line[0] for line in self.cord_block()]
	@functools.lru_cache(maxsize=1)
	def elements(self):
		return list(dict.fromkeys(self.all_elements()))
	@functools.lru_cache(maxsize=1)
	def n_electrons(self):
		return sum(elements.index(e) for e in self.all_elements())
	@functools.lru_cache(maxsize=1)
	def return_print(self):
		return [str(self.n_atoms()),self.title(),*[l for l in self.form_cord_block()]]
	def print_file(self):
		print("======={}=======".format(self.name()))
		print("=======START=======")
		print("\n".join([l for l in self.return_print()]))
		print("========END========")
	def save_file(self,directory=None):
		if directory is None:
			file_path = os.path.splitext(os.path.join(os.getcwd(),self.name().replace(" ","")))[0]+".xyz"
		else:
			file_path = os.path.splitext(os.path.join(directory,self.name().replace(" ","")))[0]+".xyz"
		if os.path.exists(file_path):
			print("File {} already exists!".format(os.path.splitext(os.path.basename(file_path))[0] + ".xyz"))
			return
		with open(file_path,"w") as file:
			for line in self.return_print():file.write(str(line)+"\n")
		print("File {} saved!".format(os.path.splitext(os.path.basename(file_path))[0] + ".xyz"))
	def print_all(self):
		print("\n".join([l for l in self.list]))
	def displace(self,mult,displacement):
		cord_block = [[a,*[float(b[n])-c[n]*mult for n in range(3)]] for a,b,c in zip(self.all_elements(),self.cord_strip(),displacement)]
		cord_block = [" ".join([str(i) for i in l]) for l in cord_block]
		return XyzFile([self.name(),self.n_atoms(),self.title(),*cord_block])
	def rotate(self, angle, axis):
		"takes xyz object and returns xyz object rotated by angle over axis"
		if axis == "x":
			m_mat = [[1., 0., 0.], [0., math.cos(angle), -math.sin(angle)], [0., math.sin(angle), math.cos(angle)]]
		if axis == "y":
			m_mat = [[math.cos(angle), 0., math.sin(angle)], [0., 1., 0.], [-math.sin(angle), 0., math.cos(angle)]]
		if axis == "z":
			m_mat = [[math.cos(angle), -math.sin(angle), 0.], [math.sin(angle), math.cos(angle), 0.], [0., 0., 1.]]
		m_mat = np.array(m_mat, np.float64)
		rotated = np.array([i[1:4] for i in self.cord_block()], np.float64).transpose()
		rotated = np.matmul(m_mat,rotated).transpose()
		rotated = np.ndarray.tolist(rotated)
		rotated = [[b,*[str(n) for n in a]] for b,a in zip(self.all_elements(),rotated)]
		xyz_mat = [self.name(), self.n_atoms()," ",*[" ".join(a) for a in rotated]]
		return XyzFile(xyz_mat)
	def superimpose(self, other, num_atoms=0, print_step=False, ret = "geom",conv=12):
		"""Takes xyz object and returns xyz object rotated by angle over axis.
		Returns either the max_distance 'max_d' or final geometry 'geom' after rotations and superpositions"""
		def rotate(xyz,angle,axis):
			if axis == "x":
				m_mat = [[1., 0., 0.], [0., math.cos(angle), -math.sin(angle)], [0., math.sin(angle), math.cos(angle)]]
			if axis == "y":
				m_mat = [[math.cos(angle), 0., math.sin(angle)], [0., 1., 0.], [-math.sin(angle), 0., math.cos(angle)]]
			if axis == "z":
				m_mat = [[math.cos(angle), -math.sin(angle), 0.], [math.sin(angle), math.cos(angle), 0.], [0., 0., 1.]]
			m_mat = np.array(m_mat, np.float64)
			rotated = np.array(xyz, np.float64).transpose()
			rotated = np.matmul(m_mat,rotated).transpose()
			return np.ndarray.tolist(rotated)
		def calc_err(xyz_1, xyz_2, n_atms):
			n_atms = len(xyz_1) if n_atms == 0 else n_atms
			sq_dist = sum(sum(math.pow(c-d,2) for c,d in zip(a,b)) for a,b in zip(xyz_1[:n_atms],xyz_2))
			return round(sq_dist / n_atms, 5)
		def max_dist(xyz_a, xyz_b):
			return max(math.sqrt(sum(pow(c-d,2) for c,d in zip(a,b))) for a,b in zip(xyz_a,xyz_b))
		#----------------------
		last_error = None
		xyz_1 = [[float(a) for a in b] for b in other.std_cord(num_atoms).cord_strip()]
		xyz_2 = [[float(a) for a in b] for b in self.std_cord(num_atoms).cord_strip()]
		if print_step: print("======ACTIONS======")
		for num in range(conv):
			step_size = 1 / 2 ** num
			while True:
				rot = [[1, "x"], [1, "y"], [1, "z"], [-1, "x"], [-1, "y"], [-1, "z"]]
				movements = [rotate(xyz_2, step_size * i[0], i[1]) for i in rot]
				if ret == "max_d":
					last_error = max_dist(xyz_2, xyz_1)
					errors = [max_dist(i, xyz_1) for i in movements]
				else:
					last_error = calc_err(xyz_2, xyz_1, num_atoms)
					errors = [calc_err(i, xyz_1, num_atoms) for i in movements]
				best_m = errors.index(min(errors))
				if min(errors) < last_error:
					xyz_2 = movements[best_m]
					last_error = min(errors)
					if print_step:
						msg = [round(step_size * rot[best_m][0], 5), rot[best_m][1], last_error]
						print("Rotating {} radian in {}. Avg. Error (A) = {}".format(*msg))
					continue
				else:
					if ret == "max_d" and max_dist(xyz_1, xyz_2) < 0.1:
						return True
					break
		if print_step: print("Final Avg. Error (A) = {}".format(round(last_error, 5)))
		if print_step: print("========END========")
		if ret == "geom":
			cord_block = [" ".join([a,*[str(n) for n in b]]) for a,b in zip(self.all_elements(),xyz_2)]
			return XyzFile([self.name(),self.n_atoms(),self.title(),*cord_block])
		elif ret == "max_d":
			return False
	def std_cord(self, n_atoms=0):
		pure_cord = self.cord_strip() if n_atoms == 0 else self.cord_strip()[0:n_atoms]
		xyz_avg = [[float(n) for n in i] for i in pure_cord]
		xyz_avg = [sum([i[n] for i in xyz_avg]) / len(xyz_avg) for n in range(3)]
		xyz_avg = [[float(i[n]) - xyz_avg[n] for n in range(3)] for i in self.cord_strip()]
		xyz_avg = [[str(n) for n in a] for a in xyz_avg]
		xyz_avg = [" ".join([b,*a]) for b,a in zip(self.all_elements(),xyz_avg)]
		xyz_mat = [self.name(), self.n_atoms(), " ", *xyz_avg]
		return XyzFile(xyz_mat)
	def enantiomer(self):
		xyz = [" ".join([*a[0:-1],str(-float(a[-1]))]) for a in self.cord_block()]
		xyz_mat = [os.path.splitext(self.name())[0]+"_ent.xyz", self.n_atoms(), " ", *xyz]
		return XyzFile(xyz_mat)
	@functools.lru_cache(maxsize=1)
	def connectivity(self):
		atoms = self.all_elements()
		rad = [dict(element_radii)[a]*0.01 for a in atoms]
		xyz = [[float(b) for b in a] for a in self.cord_strip()]
		def dist(a, b): return sum((d - c) ** 2 for c, d in zip(a, b)) ** (0.5)
		con = [[i for i,b in enumerate(xyz) if dist(a,b)<1.2*(rad[y]+rad[i]) and i!=y] for y,a in enumerate(xyz)]
		return [a for a in zip(atoms,con)]

	#TODO
	@functools.lru_cache(maxsize=1)
	def angles(self):
		aa = self.cord_strip()
		a = [[[i,aa[i]],[[n,aa[n]] for n in a[1]]] for i,a in enumerate(self.connectivity()) if len(a[1]) > 1]
		for n in a:
			b_a = np.array([float(i) for i in n[0][-1]])
			#print(n[0][-1])
			for n1 in n[1]:
				a_a = np.array([float(i) for i in n1[1]])
				#print(n1[1])
				for n2 in n[1]:
					if n2 == n1: continue
					#print(n2[1])
					c_a = np.array([float(i) for i in n2[1]])
					ba = a_a - b_a
					bc = c_a - b_a
					cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
					angle = np.arccos(cosine_angle)
					print("n1\n",n1)
					print("n\n",n)
					print("n2\n",n2)
					print(self.elements()[n1[0][0]],self.elements()[n[0][0]],self.elements()[n2[0][0]], np.degrees(angle))

	def alkene(self):
		con = self.connectivity()
		csp2 = [a for a in con if all([a[1]=="C",len(a[2])==3])]
		def only_reciprocal(a):
			csp2_idxs = [a[0] for a in csp2]
			try_one = True if len([b for b in a[2] if b in csp2_idxs]) == 1 else False
			if try_one == True:
				f = [b for b in a[2] if b in csp2_idxs]
				b = con[f[0]]
				try_two = True if len([c for c in b[2] if c in csp2_idxs]) == 1 else False
				return all([try_one,try_two])
			else: return False
		alkene = [a for a in csp2 if only_reciprocal(a)]
		if len(alkene) == 2:
			print("Double bond:",alkene)
		else: print("No double bond")
	def dist_matrix(self):
		xyz = [[float(b) for b in a] for a in self.cord_strip()]
		def dist(a, b): return sum((d - c) ** 2 for c, d in zip(a, b)) ** (0.5)
		return [[dist(a,b) for a in xyz] for b in xyz]




'''		
		for a in elms:print(a)
		black_list = []
		def build_branches(branch):
			black_list.append(branch.data)
			if branch.parent != None:
				for b in [a for a in elms[branch.data][2] if a != branch.parent.data]:
					if b in black_list: print("Atom {} is in a cycle".format(str(branch.data)))
					else: branch.add_child(elms[b][1],b)
			else:
				for b in [a for a in elms[branch.data][2] if a not in black_list]:
					branch.add_child(elms[b][1],b)
			for a in branch.child:
				build_branches(a)
			return branch
		tree = build_branches(Tree(elms[0][1],12))
		tree.print_tree()'''

