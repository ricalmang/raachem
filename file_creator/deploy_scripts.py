import shutil, os
from raachem.util.constants import cf
from raachem.util.gen_purp import file_weeder

def deploy():
	scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"scripts")
	while True:
		options = file_weeder([".py"], scripts_dir)
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
		shutil.copy(os.path.join(scripts_dir, options[int(option)]), os.path.join(cf,options[int(option)]))
