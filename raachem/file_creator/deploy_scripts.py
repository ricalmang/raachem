import os
from raachem.util.gen_purp import file_weeder, preferences

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
		with open(os.path.join(scripts_dir,options[int(option)])) as file:
			script = file.readlines()
		script = [a.replace("sub_s_name",preferences.sub_s_name) for a in script]
		script = [a.replace(".gjf", preferences.gauss_ext) for a in script]
		with open((os.path.join(os.getcwd(), options[int(option)])), "w", newline="\n") as file_b:
			for line in script:
				file_b.write(str(line))

