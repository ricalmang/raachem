import shutil, os
from raachem.util.constants import cf
from raachem.util.gen_purp import file_weeder, w_any, Var

def deploy():
	scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"scripts")
	while True:
		options = file_weeder([".py"], scripts_dir)
		options = [a for a in options if a != "__init__.py"]
		options.insert(0, "To cancel")
		print("Which script do you want to paste in the current folder?")
		for idx, i in enumerate(options):
			print("{:>3} - {}".format(idx, i))
		option = input()
		if option in [str(a) for a in range(len(options))]:
			break
	if option == "0":
		return
	else:
		with open(os.path.join(scripts_dir,option[int(option)])) as file:
			script = file.readlines()
		script = [a.replace("sub_s_name",Var().sub_s_name) for a in script]
		w_any(script,write_mod="w",filename=options[int(option)],folder=cf)
		
