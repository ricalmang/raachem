from raachem.util.constants import elements

class GjfFile:
	def __init__(self,file_content):
		self.list = file_content
	def name(self):
		if len(self.list[0]) == 0: raise Exception(".gjf or .com object has no name")
		return self.list[0]
	def empty_line_idxs(self):
		return [idx for idx,line in enumerate(self.list) if line.split() == []]
	def asterisk_line_idxs(self):
		return [idx for idx,line in enumerate(self.list) if line.split() == ["****"]]
	def route_idx(self):
		for idx,line in enumerate(self.list):
			if line.strip().startswith("#"):return idx
		raise Exception("A route section (#) should be specified for .gjf or .com files")
	def route_text(self):
		return " ".join(self.list[self.route_idx():self.title_idx()])
	def title_idx(self):
		for idx,line in enumerate(self.list):
			if idx > self.route_idx() and line.split() == []: return idx+1
	def c_m_idx(self):
		if len(self.list[self.title_idx()+2].split()) < 2:
			raise Exception("Did you provide charge and multiplicity data at line {} of file {}?".format(self.title_idx()+1,self.name()))
		return self.title_idx()+2
	def end_cord_idx(self):
		for idx,line in enumerate(self.list):
			if idx < self.c_m_idx(): continue
			if line.split() == []: return idx+1
	def charge(self):
		return int(self.list[self.c_m_idx()].split()[0])
	def multiplicity(self):
		return int(self.list[self.c_m_idx()].split()[1])
	def cord_block(self):
		cordinates = []
		for line in self.list[self.c_m_idx()+1:]:
			line = line.split()
			if len(line) == 0: break
			if len(line) != 4: continue
			if line[0] in elements:	cordinates.append(line)
			else: cordinates.append([elements[int(line[0])],*line[0:]])
		return cordinates
	def form_cord_block(self):
		return ["{:<5}{:>20f}{:>20f}{:>20f}".format(x[0],*[round(float(x[a]),6) for a in [1,2,3]]) for x in self.cord_block()]
	def cord_strip(self):
		return [line[1:] for line in self.cord_block()]
	def all_elements(self):
		return [line[0] for line in self.cord_block()]
	def elements(self):
		return list(dict.fromkeys(self.all_elements()))
	def n_electrons(self):
		return sum(elements.index(e) for e in self.all_elements()) - self.charge()
	def c_m_validate(self):
		return not self.n_electrons()%2 == self.multiplicity()%2
	def c_m_validate_txt(self):
		return "Yes" if self.c_m_validate() else "--NO!--"
	def gen_basis(self):
		return any(i in self.route_text().lower() for i in ["/gen", "gen ","genecp"])
	def declared_basis_lines(self):
		if not self.gen_basis(): return None
		idxs = [i+1 for idx,i in enumerate(self.asterisk_line_idxs()) if i < self.asterisk_line_idxs()[-1]]
		idxs.insert(0,max(i+1 for i in self.empty_line_idxs() if  i < self.asterisk_line_idxs()[-1]))
		return idxs
	def declared_basis(self):
		e_w_b = [self.list[i].split()[:-1] for i in self.declared_basis_lines()]
		return [j for i in e_w_b for j in i]
	def basis_errors(self):
		if not self.gen_basis(): return []
		#errors
		zero_last = any(self.list[i].split()[-1] == "0" for i in self.declared_basis_lines())
		miss_basis = [a for a in self.elements() if a not in self.declared_basis()]
		surpl_basis = [a for a in self.declared_basis() if a not in self.elements()]
		rep_basis = list(dict.fromkeys([a for a in self.declared_basis() if self.declared_basis().count(a) > 1]))
		#statements
		errors = []
		if not zero_last:errors.append("Missing zero at the end of basis set specification?")
		if miss_basis:errors.append("Missing basis for: {} ?".format(" ".join(miss_basis)))
		if surpl_basis:errors.append("Surplous basis for: {} ?".format(" ".join(surpl_basis)))
		if rep_basis:errors.append("Repeated basis for: {} ?".format(" ".join(rep_basis)))
		return errors
	def gen_ecp(self):
		return any(i in self.route_text().lower() for i in ["pseudo", "genecp"])
	def declared_ecp_lines(self):
		line_idx = []
		if not self.gen_ecp(): return None
		if self.gen_basis(): start_idx = self.declared_basis_lines()[-1] + 1
		else:start_idx = self.end_cord_idx()
		for idx,line in enumerate(self.list):
			if idx < start_idx: continue
			if len(line.split()) <= 1: continue
			if line.split()[-1] != "0": continue
			if all(True if a in elements else False for a in line.split()[:-1]): line_idx.append(idx)
		return line_idx
	def declared_ecp(self):
		ecps = [self.list[i].split()[:-1] for i in self.declared_ecp_lines()]
		return [j for i in ecps for j in i]
	def ecp_errors(self,heavy_e = 36):
		if not self.gen_ecp(): return []
		#errors
		zero_last = any(self.list[i].split()[-1] == "0" for i in self.declared_ecp_lines())
		miss_ecp = [a for a in self.elements() if a not in self.declared_ecp() and elements.index(a) > heavy_e]
		surpl_ecp = [a for a in self.declared_ecp() if a not in self.elements()]
		rep_ecp = list(dict.fromkeys([a for a in self.declared_ecp() if self.declared_ecp().count(a) > 1]))
		#statements
		errors = []
		if not zero_last:errors.append("Missing zero at the end of ecp set specification?")
		if miss_ecp:errors.append("Missing ecp for: {} ?".format(" ".join(miss_ecp)))
		if surpl_ecp:errors.append("Surplous ecp for: {} ?".format(" ".join(surpl_ecp)))
		if rep_ecp:errors.append("Repeated ecp for: {} ?".format(" ".join(rep_ecp)))
		return errors
	def n_proc(self):
		for line in self.list:
			line = line.lower().replace(" ","")
			if "%nprocshared=" in line:	return int(line.replace("%nprocshared=",""))
			elif "%nproc=" in line:	return int(line.replace("%nproc=",""))
		return None
	def mem(self):
		for line in self.list:
			line = line.lower().replace(" ","")
			if "%mem=" in line:
				line = line.replace("%mem=","")
				if "mb" in line: return int(line.replace("mb",""))
				elif "gb" in line: return 1000*int(line.replace("gb",""))
		return None
	def link_one_idxs(self):
		return [idx for idx,l in enumerate(self.list) if "--link1--" in l.lower()]