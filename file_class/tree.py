class Tree:
	def __init__(self,name=None,data=None, parent=None):
		self.parent = parent
		self.child = []
		self.name = name
		self.data = data

	def add_child(self, name, data):
		new_child = Tree(name, data, parent=self)
		self.child.append(new_child)
		return new_child
	def is_root(self):
		return self.parent is None
	def is_leaf(self):
		return not self.children
	def print_tree(self,text=""):
		text = text + self.name
		if not self.child: print(text)
		else:
			for idx,a in enumerate(self.child):
				if idx == 0: a.print_tree(text + "--")
				else: a.print_tree(" " * (len(text)) + "--")
	def print_con(self,text=""):
		text = text + self.name
		if not self.child: print(text)
		else:
			for idx,a in enumerate(self.child):
				a.print_con(text + "--")

