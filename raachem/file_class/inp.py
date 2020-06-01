from raachem.util.constants import elements
from raachem.file_class.xyz import XyzFile
import functools
class InpFile:
	possible_blocks = ("autoci", "basis", "casscf", "cipsi", "cis", "cim", "coords", "cpcm", "elprop", "eprnmr",
					   "freq", "geom", "loc", "md", "mdci", "method", "mp2", "mrci", "mrcc", "numgrad", "nbo",
					   "output", "pal", "paras", "rel", "plots", "rocis", "rr", "scf", "vpt2")

	def __init__(self,file_content):
		self.list = file_content
		self.list_l = [a.split() for a in file_content]
		self.str_l = [a.replace(" ","") for a in self.list]
		self.return_print = "\n".join(self.list[1:])
		#########################
		#########################
		self.comt_ls = [i for i,a in enumerate(self.str_l) if a.startswith("#")]
		self.keys_ls = [i for i,a in enumerate(self.str_l) if a.startswith("!")]
		self.ljob_ls = [i for i,a in enumerate(self.str_l) if a.startswith("%base")]
		self.njob_ls = [i for i,a in enumerate(self.str_l) if a.startswith("$new_job")]
		self.file_ls = [i for i,a in enumerate(self.str_l) if a.startswith("*xyzfile")]
		self.astk_ls = [i for i,a in enumerate(self.str_l) if a == "*"]
		self.cord_ls = [i for i,a in enumerate(self.str_l) if a.startswith("*xyz") and not i in self.file_ls]
		self.perc_ls = [i for i,a in enumerate(self.str_l) if a.startswith("%") and not i in self.ljob_ls]
		self.ends_ls = [i for i,a in enumerate(self.str_l) if a == "end"]
	@functools.lru_cache(maxsize=1)
	def name(self):
		if len(self.list[0]) == 0: raise Exception(".inp object has no name")
		return self.list[0]
	@functools.lru_cache(maxsize=1)
	def charge(self):
		return self._charge_multiplicity()[0]
	@functools.lru_cache(maxsize=1)
	def multiplicity(self):
		return self._charge_multiplicity()[1]
	@functools.lru_cache(maxsize=1)
	def n_electrons(self):
		return sum(elements.index(e) for e in self.all_elements()) - self.charge()
	@functools.lru_cache(maxsize=1)
	def n_atoms(self):
		return len(self.all_elements())
	@functools.lru_cache(maxsize=1)
	def all_elements(self):
		return [a[0] for a in self.list_l[self.c_m_idx()+1:self.end_cord_idx()]] if self.cord_ls else None
	@functools.lru_cache(maxsize=1)
	def elements(self):
		return list(dict.fromkeys(self.all_elements()))
	@functools.lru_cache(maxsize=1)
	def c_m_validate(self):
		return not self.n_electrons()%2 == self.multiplicity()%2
	@functools.lru_cache(maxsize=1)
	def c_m_validate_txt(self):
		return "Yes" if self.c_m_validate() else "--NO!--"
	@functools.lru_cache(maxsize=1)
	def n_proc(self):
		n_proc = None
		for a in self.route_text().split():
			if not a.lower().startswith("pal"):continue
			if all(b.isdigit() for b in a.lower()[3:]): n_proc = int(a.lower()[3:])
		if n_proc == 1: print("ALERT: For single core jobs the number of parallel processors 'PAL' should be omited!")
		while True:
			if "pal" not in self.block_keys(): break
			if "nprocs" not in self.block_keys()["pal"]: break
			print("WARNING: Proccessor count seems to have been provided twice: {}".format(self.name))
			if self.block_keys()["pal"].index("nprocs") + 1  == len(self.block_keys()["pal"]): break
			value = self.block_keys()[self.block_keys()["pal"].index("nprocs") + 1]
			if value.isdigit(): n_proc = int(value)
			break
		if n_proc is None: n_proc = 1
		return n_proc
	@functools.lru_cache(maxsize=1)
	def cord_block(self):
		return [*self.list_l[self.c_m_idx()+1:self.end_cord_idx()]]
	@functools.lru_cache(maxsize=1)
	def route_text(self):
		return " ".join([self.list[a].replace("!", "", 1) for a in self.keys_ls])
	@functools.lru_cache(maxsize=1)
	def c_m_idx(self):
		return min(self.cord_ls)
	@functools.lru_cache(maxsize=1)
	def end_cord_idx(self):
		return min(a for a in self.astk_ls if a > self.c_m_idx())
	#########################
	#########################
	@functools.lru_cache(maxsize=1)
	def _charge_multiplicity(self):
		try:
			charge = int(self.list[self.c_m_idx()].split()[-2])
			mult = int(self.list[self.c_m_idx()].split()[-1])
		except IndexError:
			charge = None
			mult = None
			print("Did you provide charge and multiplicity data?")
			print("Line {}:\n'{}'".format(self.c_m_idx()+1,self.list[self.c_m_idx()]))
		except ValueError:
			charge = None
			mult = None
			print("Did you properly provide charge and multiplicity data?")
			print("Line {}:\n'{}'".format(self.c_m_idx()+1,self.list[self.c_m_idx()]))
		finally:
			return [charge,mult]
	@functools.lru_cache(maxsize=1)
	def block_keys(self,possible_blocks=possible_blocks):
		d = {}
		for a in self.perc_ls:
			first_key = self.list[a].replace("%","",1).split()
			while True:
				if len(first_key) == 0: a += 1; first_key = self.list[a].split()
				else: break
			if first_key[0].lower() in possible_blocks:
				end_line = self._hunt_end(a)
				if end_line != None: d[first_key[0].lower()] = " ".join(self.list[a:end_line+1]).replace("%","",1).split()
		return d

	#########################
	#########################
	def replace_cord(self,xyz_obj):
		new = []
		for line in self.list[0:self.c_m_idx()+1]: new.append(line)
		for line in xyz_obj.form_cord_block(): new.append(line)
		for line in self.list[self.end_cord_idx():]: new.append(line)
		return InpFile(new)
	def _hunt_end(self,start):
		for idx,line in enumerate(self.list[start:]):
			if "end" in line.lower().split(): return idx + start
		return None
	def xyz_obj(self):
		return XyzFile([self.name(),self.n_atoms()," ",*[" ".join(a) for a in self.cord_block()]])



