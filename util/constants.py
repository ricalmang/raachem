import os
from raachem.util.gen_purp import w_any
cf=os.getcwd()

# Data from: Pyykko, P. and Atsumi, M., Chem. Eur. J. 2009, 15, 186.
element_radii=[
    ["None",None],['H', 32],['He', 46],[ 'Li', 133],[ 'Be', 102],[ 'B', 85],[ 'C', 75],
    ['N', 71],['O', 63],['F', 64],[ 'Ne', 67],[ 'Na', 155],[ 'Mg', 139],
    ['Al', 126],['Si', 116],['P', 111],[ 'S', 103],[ 'Cl', 99],[ 'Ar', 96],
    ['K', 196],['Ca', 171],['Sc', 148],[ 'Ti', 136],[ 'V', 134],[ 'Cr', 122],
    ['Mn', 119],['Fe', 116],['Co', 111],[ 'Ni', 110],[ 'Cu', 112],[ 'Zn', 118],
    ['Ga', 124],['Ge', 121],['As', 121],[ 'Se', 116],[ 'Br', 114],[ 'Kr', 117],
    ['Rb', 210],['Sr', 185],['Y', 163],[ 'Zr', 154],[ 'Nb', 147],[ 'Mo', 138],
    ['Tc', 128],['Ru', 125],['Rh', 125],[ 'Pd', 120],[ 'Ag', 128],[ 'Cd', 136],
    ['In', 142],['Sn', 140],['Sb', 140],[ 'Te', 136],[ 'I', 133],[ 'Xe', 131],
    ['Cs', 232],['Ba', 196],['La', 180],[ 'Ce', 163],[ 'Pr', 176],[ 'Nd', 174],
    ['Pm', 173],['Sm', 172],['Eu', 168],[ 'Gd', 169],[ 'Tb', 168],[ 'Dy', 167],
    ['Ho', 166],['Er', 165],['Tm', 164],[ 'Yb', 170],[ 'Lu', 162],[ 'Hf', 152],
    ['Ta', 146],['W', 137],['Re', 131],[ 'Os', 129],[ 'Ir', 122],[ 'Pt', 123],
    ['Au', 124],['Hg', 133],[ 'Tl', 144],[ 'Pb', 144],[ 'Bi', 151],[ 'Po', 145],
    ['At', 147],['Rn', 142],[ 'Fr', 223],[ 'Ra', 201],[ 'Ac', 186],[ 'Th', 175],
    ['Pa', 169],['U', 170],[ 'Np', 171],[ 'Pu', 172],[ 'Am', 166],[ 'Cm', 166],
    ['Bk', 168],['Cf', 168],[ 'Es', 165],[ 'Fm', 167],[ 'Md', 173],[ 'No', 176],
    ['Lr', 161],['Rf', 157],[ 'Db', 149],[ 'Sg', 143],[ 'Bh', 141],[ 'Hs', 134],
    ['Mt', 129],['Ds', 128],['Rg', 121],[ 'Cn', 122],[ 'Nh', 136],[ 'Fl', 143],
    ['Mc', 162],['Lv', 175],['Ts', 165],['Og', 157]]
elements = tuple(i[0] for i in element_radii)
keywords = \
    ['1-bromo-2-methylpropane', '1-bromooctane', '1-bromopentane', '1-bromopropane', '1-butanol',
    '1-chlorohexane', '1-chloropentane', '1-chloropropane', '1-decanol', '1-fluorooctane', '1-heptanol',
    '1-hexanol', '1-hexene', '1-hexyne', '1-iodobutane', '1-iodohexadecane', '1-iodopentane',
    '1-iodopropane', '1-nitropropane', '1-nonanol', '1-pentanol', '1-pentene', '1-propanol',
    '1-trichloroethane', '2-bromopropane', '2-butanol', '2-chlorobutane', '2-dibromoethane',
    '2-dichloroethene', '2-dimethylcyclohexane', '2-ethanediol', '2-heptanone', '2-hexanone',
    '2-methoxyethanol', '2-methyl-1-propanol', '2-methyl-2-propanol', '2-methylpentane',
    '2-methylpyridine', '2-nitropropane', '2-octanone', '2-pentanone', '2-propanol', '2-propen-1-ol',
    '2-trichloroethane', '2-trifluoroethanol', '3-methylpyridine', '3-pentanone', '4-dimethylpentane',
    '4-dimethylpyridine', '4-dioxane', '4-heptanone', '4-methyl-2-pentanone', '4-methylpyridine',
    '4-trimethylbenzene', '4-trimethylpentane', '5-nonanone', '6-dimethylpyridine', 'a-chlorotoluene',
    'aceticacid', 'acetone', 'acetonitrile', 'acetophenone', 'allcheck', 'aniline', 'anisole', 'apfd',
    'argon', 'b1b95', 'b1lyp', 'b3lyp', 'b3p86', 'b3pw91', 'b971', 'b972', 'b97d', 'b97d3', 'benzaldehyde',
    'benzene', 'benzonitrile', 'benzylalcohol', 'betanatural', 'bhandh', 'bhandhlyp', 'bromobenzene',
    'bromoethane', 'bromoform', 'butanal', 'butanoicacid', 'butanone', 'butanonitrile', 'butylamine',
    'butylethanoate', 'calcall', 'calcfc', 'cam-b3lyp', 'carbondisulfide', 'carbontetrachloride',
    'cartesian', 'checkpoint', 'chkbasis', 'chlorobenzene', 'chloroform', 'cis-1', 'cis-decalin',
    'connectivity', 'counterpoise', 'cyclohexane', 'cyclohexanone', 'cyclopentane', 'cyclopentanol',
    'cyclopentanone', 'd95v', 'decalin-mixture', 'def2qzv', 'def2qzvp', 'def2qzvpp', 'def2sv', 'def2svp',
    'def2svpp', 'def2tzv', 'def2tzvp', 'def2tzvpp', 'density', 'densityfit', 'dibromomethane',
    'dibutylether', 'dichloroethane', 'dichloromethane', 'diethylamine', 'diethylether', 'diethylsulfide',
    'diiodomethane', 'diisopropylether', 'dimethyldisulfide', 'dimethylsulfoxide', 'diphenylether',
    'dipropylamine', 'e-2-pentene', 'empiricaldispersion', 'ethanethiol', 'ethanol', 'ethylbenzene',
    'ethylethanoate', 'ethylmethanoate', 'ethylphenylether', 'extrabasis', 'extradensitybasis', 'finegrid',
    'fluorobenzene', 'formamide', 'formicacid', 'freq', 'full', 'gd3bj', 'genecp', 'geom', 'gfinput',
    'gfprint', 'hcth', 'hcth147', 'hcth407', 'hcth93', 'heptane', 'hexanoicacid', 'hissbpbe', 'hseh1pbe',
    'integral', 'iodobenzene', 'iodoethane', 'iodomethane', 'isopropylbenzene', 'isoquinoline', 'kcis',
    'krypton', 'lanl2dz', 'lanl2mb', 'lc-wpbe', 'loose', 'm-cresol', 'm-xylene', 'm062x', 'm06hf', 'm06l',
    'm11l', 'maxcycles', 'maxstep', 'mesitylene', 'methanol', 'methylbenzoate', 'methylbutanoate',
    'methylcyclohexane', 'methylethanoate', 'methylmethanoate', 'methylpropanoate', 'minimal', 'mn12l',
    'mn12sx', 'modredundant', 'mpw1lyp', 'mpw1pbe', 'mpw1pw91', 'mpw3pbe', 'n-butylbenzene', 'n-decane',
    'n-dimethylacetamide', 'n-dimethylformamide', 'n-dodecane', 'n-hexadecane', 'n-hexane',
    'n-methylaniline', 'n-methylformamide-mixture', 'n-nonane', 'n-octane', 'n-octanol', 'n-pentadecane',
    'n-pentane', 'n-undecane', 'n12sx', 'nitrobenzene', 'nitroethane', 'nitromethane', 'noeigentest',
    'nofreeze', 'noraman', 'nosymm', 'nprocshared', 'o-chlorotoluene', 'o-cresol', 'o-dichlorobenzene',
    'o-nitrotoluene', 'o-xylene', 'o3lyp', 'ohse1pbe', 'ohse2pbe', 'oniom', 'output', 'p-isopropyltoluene',
    'p-xylene', 'pbe1pbe', 'pbeh', 'pbeh1pbe', 'pentanal', 'pentanoicacid', 'pentylamine',
    'pentylethanoate', 'perfluorobenzene', 'pkzb', 'population', 'propanal', 'propanoicacid',
    'propanonitrile', 'propylamine', 'propylethanoate', 'pseudo', 'pw91', 'pyridine', 'qst2', 'qst3',
    'quinoline', 'qzvp', 'rdopt', 'read', 'readfc', 'readfreeze', 'readopt', 'readoptimize', 'regular',
    'restart', 's-dioxide', 'savemixed', 'savemulliken', 'savenbos', 'savenlmos', 'scrf', 'sddall',
    'sec-butylbenzene', 'sogga11', 'sogga11x', 'solvent', 'spinnatural', 'tert-butylbenzene',
    'tetrachloroethene', 'tetrahydrofuran', 'tetrahydrothiophene-s', 'tetralin', 'thcth', 'thcthhyb',
    'thiophene', 'thiophenol', 'tight', 'toluene', 'tpss', 'tpssh', 'trans-decalin', 'tributylphosphate',
    'trichloroethene', 'triethylamine', 'tzvp', 'ultrafine', 'uncharged', 'v5lyp', 'verytight', 'vp86',
    'vsxc', 'vwn5', 'water', 'wb97', 'wb97x', 'wb97xd', 'wpbeh', 'x3lyp', 'xalpha', 'xenon',
    'xylene-mixture']
subshells = ("1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p",
             "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p", "8s", "5g", "6f",
             "7d", "8p", "9s")
subshell_size = (2,2,6,2,6,2,10,6,2,10,6,2,14,10,6,2,14,10,6,2,18,14,10,6,2)
