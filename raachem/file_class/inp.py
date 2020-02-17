from raachem.util.constants import elements
from raachem.file_class.xyz import XyzFile
class InpFile:
	possible_blocks = ("autoci", "basis", "casscf", "cipsi", "cis", "cim", "coords", "cpcm", "elprop", "eprnmr",
					   "freq", "geom", "loc", "md", "mdci", "method", "mp2", "mrci", "mrcc", "numgrad", "nbo",
					   "output", "pal", "paras", "rel", "plots", "rocis", "rr", "scf", "vpt2")

	def __init__(self,file_content):
		self.list = file_content
		self.list_l = [a.split() for a in file_content]
		self.str_l = [a.replace(" ","") for a in self.list]
		self.name = self.str_l[0]
		if len(self.name) == 0:	raise Exception("Object has no name")
		self.comt_ls = [i for i,a in enumerate(self.str_l) if a.startswith("#")]
		self.keys_ls = [i for i,a in enumerate(self.str_l) if a.startswith("!")]
		self.ljob_ls = [i for i,a in enumerate(self.str_l) if a.startswith("%base")]
		self.njob_ls = [i for i,a in enumerate(self.str_l) if a.startswith("$new_job")]
		self.file_ls = [i for i,a in enumerate(self.str_l) if a.startswith("*xyzfile")]
		self.astk_ls = [i for i,a in enumerate(self.str_l) if a == "*"]
		self.cord_ls = [i for i,a in enumerate(self.str_l) if a.startswith("*xyz") and not i in self.file_ls]
		self.perc_ls = [i for i,a in enumerate(self.str_l) if a.startswith("%") and not i in self.ljob_ls]
		self.ends_ls = [i for i,a in enumerate(self.str_l) if a == "end"]
		self.start_xyz = min(self.cord_ls)
		self.end_xyz = min(a for a in self.astk_ls if a > self.start_xyz)
		self.keys = " ".join([self.list[a].replace("!","",1) for a in self.keys_ls])
		self.atoms = [a[0] for a in self.list_l[self.start_xyz+1:self.end_xyz]] if self.cord_ls else None
		self.unique_atoms = list(dict.fromkeys(self.atoms)) if self.cord_ls else None
		self.n_atoms = len(self.atoms) if self.cord_ls else None
		self.block_keys = self._block_keys()
		self.n_proc = self._n_proc()
		self.charge = self._charge_multiplicity()[0]
		self.mult = self._charge_multiplicity()[1]
		self.n_electrons = sum(elements.index(e) for e in self.atoms) - self.charge
		self.c_m_validate = self.n_electrons % 2 == self.mult % 2  # True if valid, otherwise False
		self.c_m_validate_txt = "Yes" if self.c_m_validate else "--NO!--"
		self.return_print = "\n".join(self.list[1:])

	def __str__(self):
		for line in self.list[1:]: print(line)
	def _charge_multiplicity(self):
		try:
			charge = int(self.list[self.start_xyz].split()[-2])
			mult = int(self.list[self.start_xyz].split()[-1])
		except IndexError:
			charge = None
			mult = None
			print("Did you provide charge and multiplicity data?")
			print("Line {}:\n'{}'".format(self.start_xyz+1,self.list[self.start_xyz]))
		except ValueError:
			charge = None
			mult = None
			print("Did you properly provide charge and multiplicity data?")
			print("Line {}:\n'{}'".format(self.start_xyz+1,self.list[self.start_xyz]))
		return [charge,mult]
	def _n_proc(self):
		n_proc = None
		for a in self.keys.split():
			if not a.lower().startswith("pal"):continue
			if all(b.isdigit() for b in a.lower()[3:]): n_proc = a.lower()[3:]
		while True:
			if "pal" not in self.block_keys: break
			if "nprocs" not in self.block_keys["pal"]: break
			print("WARNING: Proccessor count seems to have been provided twice: {}".format(self.name))
			if self.block_keys["pal"].index("nprocs") + 1  == len(self.block_keys["pal"]): break
			value = self.block_keys[self.block_keys["pal"].index("nprocs") + 1]
			if all(a.isdigit() for a in value): n_proc = value
			break
		return n_proc
	def _block_keys(self,possible_blocks=possible_blocks):
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
	def xyz_obj(self):
		return XyzFile([self.name,self.n_atoms," ",*self.list[self.start_xyz+1:self.end_xyz]])
	def _hunt_end(self,start):
		for idx,line in enumerate(self.list[start:]):
			if "end" in line.lower().split(): return idx + start
		return None
	def replace_cord(self,xyz_obj):
		new = []
		for line in self.list[0:self.start_xyz+1]: new.append(line)
		for line in xyz_obj.form_cord_block(): new.append(line)
		for line in self.list[self.end_xyz:]: new.append(line)
		return InpFile(new)
	def __print__(self):
		print("name", getattr(self, "name"))
		print("charge", getattr(self, "charge"))
		print("keys", getattr(self, "keys"))
		print("mult", getattr(self, "mult"))
		print("atoms", getattr(self, "atoms"))
		print("unique_atoms",getattr(self,"unique_atoms"))
		print("start_xyz",getattr(self,"start_xyz"))
		print("end_xyz", getattr(self, "end_xyz"))
		print("n_proc", getattr(self, "n_proc"))




