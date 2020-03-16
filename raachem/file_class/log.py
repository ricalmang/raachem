from raachem.file_class.xyz import XyzFile
from raachem.util.constants import elements
import functools

class LogFile:
	def __init__(self,file_content):
		self.list = file_content
		self.s_list = [a.split() for a in self.list]
		assert any("Gaussian" in l for l in self.list[:10]),\
			"Are you sure {} is a Gaussian log file?".format(self.name())
	@functools.lru_cache(maxsize=1)
	def name(self):
		if len(self.list[0]) == 0:raise Exception(".log Object has no name")
		else: return self.list[0]
	@functools.lru_cache(maxsize=1)
	def charge_mult(self):
		pattern = ("Charge", "=","value", "Multiplicity","=", "value")
		nums = (0,1,3,4)
		for i,a in enumerate(self.list):
			if i > 150: print("CM not found\nIt will be assumed to be 0 1");return ["0", "1"]
			elif len(a.split()) != 6: continue
			elif not all([b == c for d,(b,c) in enumerate(zip(a.split(),pattern)) if d in nums]): continue
			else: return a.split()[2::3]
	@functools.lru_cache(maxsize=1)
	def start_xyz_idxs(self):
		indexes = [i+5 for i,l in enumerate(self.s_list) if len(l) == 2 and "orientation:" == l[-1]]
		if len(indexes) == 0: raise Exception("No cartesian coordinates found for {}!".format(self.name()))
		return indexes
	@functools.lru_cache(maxsize=1)
	def end_xyz_idxs(self):
		one_el = [i for i,l in enumerate(self.s_list) if len(l) == 1]
		end_idxs = [min(a for a in one_el if a > num) for num in self.start_xyz_idxs()]
		if len(end_idxs) != len(self.start_xyz_idxs()): raise Exception("Inconsistent number of geometries in object {}!".format(self.name))
		return end_idxs
	@functools.lru_cache(maxsize=1)
	def xyz_cord_block(self,start_idx,end_idx):
		return [[elements[int(l[1])],*[l[i] for i in [3,4,5]]] for l in self.s_list[start_idx:end_idx]]
	@functools.lru_cache(maxsize=1)
	def last_cord_block(self):
		return self.xyz_cord_block(self.start_xyz_idxs()[-1],self.end_xyz_idxs()[-1])
	@functools.lru_cache(maxsize=1)
	def first_cord_block(self):
		return self.xyz_cord_block(self.start_xyz_idxs()[0],self.end_xyz_idxs()[0])
	@functools.lru_cache(maxsize=1)
	def n_atoms(self):
		return len(self.first_cord_block())
	def any_xyz_obj(self,a_idx,b_idx,title=" ",name=False):
		if name == False: name = self.name()
		return XyzFile([name, self.n_atoms(), title, *(" ".join(l) for l in self.xyz_cord_block(a_idx,b_idx))])
	@functools.lru_cache(maxsize=1)
	def last_xyz_obj(self):
		return XyzFile([self.name(),self.n_atoms()," ",*(" ".join(l) for l in self.last_cord_block())])
	@functools.lru_cache(maxsize=1)
	def first_xyz_obj(self):
		return XyzFile([self.name(),self.n_atoms()," ",*(" ".join(l) for l in self.first_cord_block())])
	@functools.lru_cache(maxsize=1)
	def scf_done(self):
		scf = [[idx,l.split()[4]] for idx,l in enumerate(self.list) if "SCF Done:" in l]
		if len(scf) == 0: raise Exception("Could not find string 'SCF Done:' in {} file".format(self.name()))
		return scf
	@functools.lru_cache(maxsize=1)
	def first_thermal(self):
		thermal = ["ZPC","TCE","TCH","TCG","SZPE","STE","STH","STG"]
		freq = False
		for line in self.list:
			if "Zero-point correction=" in line: freq = True
			if not freq: continue
			if "Zero-point correction=" in line: thermal[0] = line.split()[-2]
			elif "Thermal correction to Energy=" in line:thermal[1] = line.split()[-1]
			elif "Thermal correction to Enthalpy=" in line:thermal[2] = line.split()[-1]
			elif "Thermal correction to Gibbs Free Energy=" in line:thermal[3] = line.split()[-1]
			elif "Sum of electronic and zero-point Energies=" in line:thermal[4] = line.split()[-1]
			elif "Sum of electronic and thermal Energies=" in line:thermal[5] = line.split()[-1]
			elif "Sum of electronic and thermal Enthalpies=" in line:thermal[6] = line.split()[-1]
			elif "Sum of electronic and thermal Free Energies=" in line:thermal[7] = line.split()[-1]
			if all(a != i for a,i in zip(("ZPC","TCE","TCH","TCG","SZPE","STE","STH","STG"),thermal)):
				return [float(i) for i in thermal]
		return False
	@functools.lru_cache(maxsize=1)
	def displ_block(self):
		block = ["start","end"]
		for idx_s,line in enumerate(self.list):
			if not "Frc consts  -- " in line: continue
			block[0] = idx_s-4
			for idx_e,line in enumerate(self.list[idx_s:]):
				if len(line.split()) != 0: continue
				block[1] = idx_e+idx_s
				break
			break
		if block[0] == "start" or block[1] == "end": return False
		else: return self.list[block[0]:block[1]]
	@functools.lru_cache(maxsize=1)
	def frequencies(self):
		if self.displ_block() == False:	return False
		return [j for a in [i.split()[2:] for i in self.displ_block()[2::self.n_atoms()+7]] for j in a]
	@functools.lru_cache(maxsize=1)
	def displ_for_freq_idx(self,freq_idx):
		if self.displ_block() == False:	return False
		displ = []
		for num in range(self.n_atoms()):
			displ.append([j for a in [i.split()[2:] for i in self.displ_block()[7+num::self.n_atoms()+7]] for j in a])
		displ_for_freq_str = [a[freq_idx*3:freq_idx*3+3] for a in displ]
		displ_for_freq_float = [[float(i) for i in b] for b in displ_for_freq_str]
		return displ_for_freq_float
	@functools.lru_cache(maxsize=1)
	def oc_orb_energies(self):
		oc_orb_energies = []
		for i in [x for l in [line.split() for line in self.list if "occ. eigenvalues" in line] for x in l]:
			try: oc_orb_energies.append(float(i))
			except: pass
		return oc_orb_energies
	@functools.lru_cache(maxsize=1)
	def uno_orb_energies(self):
		uno_orb_energies = []
		for i in [x for l in [line.split() for line in self.list if "virt. eigenvalues" in line] for x in l]:
			try: uno_orb_energies.append(float(i))
			except: pass
		return uno_orb_energies
	@functools.lru_cache(maxsize=1)
	def homo_lumo_gap(self):
		return min(self.uno_orb_energies())-max(self.oc_orb_energies())
	@functools.lru_cache(maxsize=1)
	def calc_type(self):
		up_to = 100 if len(self.list) > 100 else len(self.list)
		dashes = [i for i,a in enumerate(self.list[:up_to]) if a.lstrip().startswith("--------")]
		route = [i for i,a in enumerate(self.list[:up_to]) if a.lstrip().startswith("#")][0]
		r_sect = self.list[max(a for a in dashes if a < route)+1:min(a for a in dashes if a > route)]
		for x in [None, "/", "(", ")", ",", "=", "%", ":"]:
			r_sect = [a for b in [i.split(x) for i in r_sect] for a in b if len(a) > 1]
		r_sect = [a.lower() for a in r_sect]
		if "ts" in r_sect: return "TS"
		elif any(True for a in r_sect if a in ("modredundant", "readoptimize", "readfreeze")): return "Red"
		elif "irc" in r_sect: return "IRC"
		elif "opt" in r_sect: return "Opt"
		else: return "SP"
	@functools.lru_cache(maxsize=1)
	def normal_termin(self):
		return any(True if "Normal termination of Gaussian" in l else False for l in self.list[-5:])
	@functools.lru_cache(maxsize=1)
	def irc(self):
		points = []
		for i,a in enumerate(self.s_list):
			if len(a) != 6: continue
			if any(c != a[b] for c,b in zip(["Point","Number:","Path","Number:"],[0,1,3,4])): continue
			if any(not a[b].isnumeric() for b in [2,5]): continue
			else: points.append([i,a[5],a[2]])
		scf = [max(self.scf_done(),key=lambda x: x[0] if x[0] < a[0] else 0)[1] for a in points]
		a_idx = [max(self.start_xyz_idxs(),key=lambda x: x if x < a[0] else 0) for a in points]
		b_idx = [max(self.end_xyz_idxs(),key=lambda x: x if x < a[0] else 0) for a in points]
		points = [[*d[1:],c,self.any_xyz_obj(a,b,title=c)] for a,b,c,d in zip(a_idx,b_idx,scf,points)]
		path_a = sorted([a for a in points if a[0] == "1"], key = lambda x: int(x[1]), reverse=True)
		path_b = [a for a in points if a[0] == "2"]
		return [a[3] for a in [*path_a,*path_b]]
	@functools.lru_cache(maxsize=1)
	def scan_geoms(self):
		geoms = []
		for idx,line in enumerate(self.s_list):
			if len(line) != 4: continue
			if any(a != line[b] for a,b in zip(("!","Parameters","!"),(0,2,3))):continue
			if all(line[1] != i for i in ("Optimized","Non-Optimized")):continue
			try:
				start_idx = min(i for i in self.start_xyz_idxs() if i > idx)
				end_idx = min(i for i in self.end_xyz_idxs() if i > idx)
				scf_idx = max(i for i in self.scf_done() if i[0] < idx)
				name = self.name().replace(".log","_" + str(len(geoms)+1)+".xyz")
				if line[1] == "Optimized": print("Optimized geometry found at line {}!".format(idx+1))
				elif line[1] == "Non-Optimized": print("Non-Optimized1 geometry found at line {}!".format(idx+1))
				geoms.append(self.any_xyz_obj(start_idx,end_idx,title=str(scf_idx[1]), name=name))
			except ValueError:
				print("WARNING: Linked Jobs may or may not be appended")
				pass
		if len(geoms) == 0:
			print("No Optimized geometries found for {} file".format(self.name()))
		return geoms