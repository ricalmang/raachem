import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="raachem",
    version="0.0.6",
    author="Ricardo Almir Angnes",
    author_email="ricardo_almir@hotmail.com",
    description="Tools I use for computational chemistry",
    long_description=long_description,
	license="MIT",
    url="https://github.com/ricalmang/raachem",
	keywords = ['chemistry'],
	install_requires = ["numpy"],
    packages=setuptools.find_packages(),
	package_data = {"":["*.GAUSSIAN","*.ORCA","*.gbs"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)