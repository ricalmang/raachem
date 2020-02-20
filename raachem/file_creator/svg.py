import os, math
from raachem.util.constants import cf
from raachem.util.gen_purp import read_item, is_str_float

class SvgGen:
	def __init__(self,text,wide=10,spam_request=True):
		self.text = text  # ".txt" input
		self.wide = [wide, 60 - wide]  # Inverse widness of platform
		self.spam_request = spam_request  # Did user requested spam
		self.name = self.text[0].replace(".txt","")  # filename
		self.svg_code = ['<?xml version="1.0" encoding="UTF-8" ?>']  # svg code as a list
		self.dg = None  # reaction free energy
		self.spam_worthy = True  # Is the input worthy of spam calculation
		self.title = ["Main_title", "Energy_level", "Reaction_cordinate"]  # Axis and table
		# The "self.graphs" is a quadruply nested list
		# [[g_1],...,[g_3]][0] => [cycle,color,[dots]][-1] => [[a],...,[b]][-1] => [5,"TS_3",7.32,height]
		self.graphs = []
		self.conectors = []
		self.complete_init()
		assert len(self.graphs) != 0, "Could not obtain energies from file!"
		if len(self.graphs) == 1 and len(self.graphs[0][-1]) == 1: self.spam_worthy = False
		self.set_height()
		self.graph_frame()
		self.graph_grid()
		self.graph_connectors()
		self.graph_crt_points()
		self.spam_dg()
		self.graph_spam()
	def complete_init(self):
		dots = []  # [lane,name,Energy,0]
		two_colum_idx = 0
		for line in [a.split() for a in self.text[1:]]:
			if any([len(line) > 6,*[len(line) == a for a in [0,1,5]]]):
				continue
			elif line[0].startswith("#tit") and len(line) == 4:
				self.title = line[1:]
			elif line[0].startswith("#con"):
				self.spam_worthy = False
				self.conectors.append(line)
			elif line[0].startswith("#end"):
				cycle = int(line[1])
				color = line[2]
				if dots:
					dots.sort(key=lambda x: x[0])
					self.graphs.append([cycle, color, dots])
					dots = []
					two_colum_idx = 0
				else:
					print("#end should only be placed at block ends!")
			elif len(line) == 2:
				if is_str_float(line[1]):
					two_colum_idx += 1
					dots.append([two_colum_idx, line[0], float(line[1].replace(",", ".")), 0])
			elif len(line) == 3:
				if is_str_float(line[2]) and line[0].isdigit():
					dots.append([int(line[0]), line[1], float(line[2].replace(",", ".")), 0])
			else:
				pass
		# misuse and option handling
		if len(self.graphs) == 0:
			cycle, color = 1, "black"
			if len(dots) == 0:
				print("No data found for: {}".format(self.name))
				return
			else:
				dots.sort(key=lambda x: x[0])
				self.graphs.append([cycle, color, dots])
		for idx, i in enumerate(self.title):
			if i.lower() == "none":
				self.title[idx] = " "
		if len(self.graphs) == 1:
			self.dg = float(self.graphs[0][-1][-1][2]) - float(self.graphs[0][-1][0][2])
			if self.dg > 0:
				self.spam_worthy = False
				print("The reaction is endergonic!")
				print("No spam will be calculated!")
		else:
			self.spam_worthy = False
	def all_e(self):
		return [x for y in [i[2] for i in self.graphs] for x in y]
	def n_col(self):
		return max(i[0] for i in self.all_e())
	def min_g(self):
		if len(self.graphs) == 1 and len(self.graphs[0][-1]) == 1: return min(i[2] for i in self.all_e() if len(i) == 4) - 10
		return min(i[2] for i in self.all_e() if len(i) == 4)
	def max_g(self):
		if len(self.graphs) == 1 and len(self.graphs[0][-1]) == 1: return max(i[2] for i in self.all_e() if len(i) == 4) + 10
		return max(i[2] for i in self.all_e() if len(i) == 4)
	def range_dg(self):
		return self.max_g() - self.min_g()
	def set_height(self):
		for idx_a,a in enumerate(self.graphs):  #for every block
			for idx_b,b in enumerate(a[2]):  #For every structure
				height = int(round(abs(400 - (b[2] - self.min_g()) * 400 / self.range_dg())))
				self.graphs[idx_a][2][idx_b][3] = height
	def spam_dg(self):
		if self.spam_worthy:
			two_cycles = []
			for idx, a in enumerate(self.graphs[0][-1] * 2):
				if len(self.graphs[0][-1]) > idx:
					two_cycles.append(a)
				elif len(self.graphs[0][-1]) == idx:
					continue
				else:
					two_cycles.append([a[0], a[1], a[2] + self.dg, a[3]])
			ind_spam = []
			for idx, a in enumerate(self.graphs[0][-1]):
				all_it = [[a[0], a[1], b[2] - a[2], b[0]] for idx_b, b in enumerate(two_cycles) if idx_b > idx]
				ind_spam.append(max(all_it,key=lambda c: c[2]))
			return max(ind_spam, key=lambda a: a[2])
	def graph_frame(self):
		a=[
			'<svg width="{0}" viewBox="30 0 {0} 500" height="500" xmlns="http://www.w3.org/2000/svg">',
			'    <line x1="100" y1="25" x2="100" y2="475" stroke="black" stroke-width="2"/>',
			'    <line x1="100" y1="475" x2="{}" y2="475" stroke="black" stroke-width="2"/>',
			'    <text x="{}" y="20" font-size="22" text-anchor="middle" fill="black">{}</text>',
			'    <text x="-250" y="55" font-size="22" {} text-anchor="middle" fill="black">{}</text>',
			'    <text x="{}" y="495" font-size="22" text-anchor="middle" fill="black">{}</text>']
		a[0]=a[0].format((self.n_col() + 1) * 80 + 100)
		a[2]=a[2].format((self.n_col() + 1) * 80 + 75)
		a[3]=a[3].format(int(self.n_col() * 40 + 80), self.title[0].replace("_", " "))
		a[4]=a[4].format('transform="rotate(-90)"',self.title[1].replace("_", " "))
		a[5]=a[5].format(int(self.n_col() * 40 + 80), self.title[2].replace("_", " "))
		self.svg_code.extend(a)
	def graph_grid(self):
		step_size = max(a if 10*a < abs(self.range_dg()) else 0.05 for a in [0.1,0.5,1,2,5,10,100])
		max_e = round((math.ceil(self.max_g())/10))*10
		steps = [max_e]
		while True:
			next_value = steps[-1]-step_size
			steps.append(next_value)
			if next_value < self.min_g(): break
		for item in steps:
			value=int(round(450-400*((item-self.min_g())/(self.range_dg()))))
			if 24 < value < 476:
				b = [
					'    <line x1="100" y1="{0}" x2="105" y2="{0}" stroke="black" stroke-width="2"/>',
					'    <text x="80" y="{}" text-anchor="middle" fill="black">{:.1f}</text>']
				b[0]=b[0].format(value)
				b[1]=b[1].format(value, item)
				self.svg_code.extend(b)
	def graph_crt_points(self):
		for i in self.graphs:
			if not len([a[0] for a in i[-1]]) == len(set([a[0] for a in i[-1]])):
				print("WARNING: Two or more structures are occupying the same block lane!")
			l_c = [0, 0, 0] # last collumn
			for idx, item in enumerate(i[2]):
				c_p = [int((item[0]+1)*80+self.wide[0]),int(round(item[3]+50)),int((item[0]+1)*80+self.wide[1])]
				# [x1,y1,x2], y1=y2
				a = [
					'    <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="4"/>',
					'    <text x="{}" y="{}" text-anchor="middle" fill="{}">{}</text>',
					'    <text x="{}" y="{}" text-anchor="middle" fill="{}">{}</text>']
				a[0]=a[0].format(c_p[0], c_p[1], c_p[2], c_p[1], i[1])
				a[1]=a[1].format(c_p[0] + 20, c_p[1] - 5, i[1], item[1].replace("_", " "))
				a[2]=a[2].format(c_p[0] + 20, c_p[1] + 15, i[1], item[2])
				self.svg_code.extend(a)
				if not idx == 0:
					b='    <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="2"/>'
					b=b.format(l_c[2], l_c[1], c_p[0], c_p[1], i[1])
					self.svg_code.append(b)
				l_c = c_p
	def graph_spam(self):
		if self.spam_worthy and self.spam_request:
			data = [self.graphs[0][2][self.spam_dg()[i] - 1] for i in [0, 3]]
			data.sort(key=lambda a: a[3], reverse=True)
			# TDI arrow
			p = [(data[0][0] + 1) * 80 + 10, data[0][3] - 40]
			a = [
				'    <text x="{}" y="{}" text-anchor="middle" fill="black">TDI</text>',
				'    <path d=" M {0} {1} L {2} {1} L {2} {3} L {4} {3} L {5} {6} L {7} {3} L {0} {3} Z "/>']
			a[0]=a[0].format(p[0]+20,p[1])
			a[1]=a[1].format(10+p[0],10+p[1],30+p[0],40+p[1],40+p[0],20+p[0],70+p[1],0+p[0])
			self.svg_code.extend(a)
			# TDTS arrow
			p = [(data[1][0] + 1) * 80 + 50, data[1][3] + 140]
			a = [
				'    <text x="{}" y="{}" text-anchor="middle" fill="black">TDTS</text>',
				'    <path d=" M {0} {1} L {2} {1} L {2} {3} L {4} {3} L {5} {6} L {7} {3} L {0} {3} Z "/>']
			a[0]=a[0].format(p[0]-20,p[1]+10)
			a[1]=a[1].format(-10+p[0],-10+p[1],-30+p[0],-40+p[1],-40+p[0],-20+p[0],-70+p[1],0+p[0])
			self.svg_code.extend(a)
			# spam and dg anotations
			a = [
				'    <text x="120" y="450" text-anchor="left" fill="black">Gr = {:.2f}</text>',
				'    <text x="120" y="470" text-anchor="left" fill="black">Es = {:.2f}</text>']
			a[0]=a[0].format(self.dg)
			a[1]=a[1].format(self.spam_dg()[2])
			self.svg_code.extend(a)
	def graph_connectors(self):
		for i in self.conectors:
			i = [int(a) if idx in [1, 2, 3, 4] else a for idx, a in enumerate(i)]
			con = []
			for line in self.graphs[i[1] - 1][2]:
				if line[0] == i[2]:
					con.append(line)
			for line in self.graphs[i[3] - 1][2]:
				if line[0] == i[4]:
					con.append(line)
			if con[0] > con[1]:
				con.reverse()
			if con[0] < con[1]:
				a='    <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="2"/>'
				a=a.format(int((con[0][0]+1)*80+self.wide[1]),con[0][3]+50,int((con[1][0]+1)*80+self.wide[0]),con[1][3]+50,i[5])
				self.svg_code.append(a)
			if con[0] == con[1]:
				print("Cannot conect items on same column")
	
	def save_svg(self):
		with open(os.path.join(cf, "E_profile.svg"), "w") as out_file:
			for line in self.svg_code: out_file.write(line + "\n")
			out_file.write('</svg>')
		print("Take a look at file E_profile.svg!")

def vector_graph():
	while True:
		text = read_item(None, "Which txt item contains the names and energies?", [".txt"])
		if text: break
		else: return
	try:
		SvgGen(text).save_svg()
	except AssertionError:
		return
