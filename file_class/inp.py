from raachem.util.constants import elements
class InpFile:
	def __init__(self,file_content):
		self.list = file_content
		self.str_l = [a.replace(" ","") for a in self.list]
		self.name = self.str_l[0]
		if len(self.name) == 0:	raise Exception(".inp Object has no name")
		self.comt_ls = [i for i,a in enumerate(self.str_l) if a.startswith("#")]
		self.keys_ls = [i for i,a in enumerate(self.str_l) if a.startswith("!")]
		self.ljob_ls = [i for i,a in enumerate(self.str_l) if a.startswith("%base")]
		self.njob_ls = [i for i,a in enumerate(self.str_l) if a.startswith("$new_job")]
		self.file_ls = [i for i,a in enumerate(self.str_l) if a.startswith("*xyzfile")]
		self.astk_ls = [i for i,a in enumerate(self.str_l) if a == "*"]
		self.cord_ls = [i for i,a in enumerate(self.str_l) if a.startswith("*xyz") and not i in self.file_ls]
		self.perc_ls = [i for i,a in enumerate(self.str_l) if a.startswith("%") and not i in self.ljob_ls]
		self.ends_ls = [i for i,a in enumerate(self.str_l) if a == "end"]
		self.first_xyz = min([*self.file_ls,*self.cord_ls])
		try:
			self.charge = int(self.list[self.first_xyz].split()[-2])
			self.mult = int(self.list[self.first_xyz].split()[-1])
		except IndexError:
			print("Did you provide charge and multiplicity data?")
			print("Line {}:\n'{}'".format(self.first_xyz+1,self.list[self.first_xyz]))
		except ValueError:
			print("Did you properly provide charge and multiplicity data?")
			print("Line {}:\n'{}'".format(self.first_xyz+1,self.list[self.first_xyz]))
		if self.cord_ls:
			self.atoms = []
			self.unique_atoms = []
