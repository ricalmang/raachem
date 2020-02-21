from raachem import read_item
from raachem import InpFile, XyzFile, GjfFile, LogFile
import os
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"testes")
inp = InpFile(read_item(file_name="teste.001.inp", promp = False, cf=test_dir))
gjf = GjfFile(read_item(file_name="teste.001.gjf", promp = False, cf=test_dir))
assert isinstance(inp, InpFile)
assert isinstance(gjf, GjfFile)
a = None
def wal(b): global a; a = b; return b
for obj in [inp,gjf]:
	assert wal(obj.name()[:-4]) == "teste.001", print(a)
	assert wal(obj.charge()) == 1, print(a)
	assert wal(obj.multiplicity()) == 1, print(a)
	assert wal(obj.n_atoms()) == 2, print(a)
	assert wal(obj.n_electrons()) == 9, print(a)
	assert wal(type(obj.c_m_idx())) is int, print(a)
	assert all(wal(a == b) for a,b in zip(obj.elements(),["H","F"])), print(a)
	assert wal(obj.c_m_validate() )is False
	assert wal(obj.c_m_validate_txt()) == "--NO!--", print(a)
	assert wal(obj.n_proc()) == 8, print(a)
	assert wal(type(obj.route_text())) is str, print(a)
	assert wal(type(obj.cord_block())) is list, print(a)
	assert wal(type(obj.cord_block()[0])) is list, print(a)
	assert all(wal(type(a) is str) for a in obj.cord_block()[0]), print(a)
	assert wal(isinstance(obj.xyz_obj(),XyzFile)), print(a)
	assert wal(type(obj.end_cord_idx())) is int, print(a)
	assert wal(obj.xyz_obj().n_atoms()) == obj.n_atoms(), print(a)
	assert wal(obj.xyz_obj().name()) == obj.name(), print(a)
	assert wal(obj.xyz_obj().cord_block()) == obj.cord_block(), print(a)
	assert wal(obj.xyz_obj().elements()) == obj.elements(), print(a)
	assert wal(obj.xyz_obj().all_elements()) == obj.all_elements(), print(a)
	assert wal(obj.replace_cord(obj.xyz_obj()).n_electrons()) == obj.n_electrons(), print(a)
	assert wal(type(obj.list)) is list, print(a)
	assert wal(type(obj.list[0])) is str, print(a)
	assert wal(type(obj.list_l)) is list, print(a)
	assert wal(type(obj.list_l[0])) is list, print(a)
assert wal(type(gjf.title_idx())) is int, print(a)

