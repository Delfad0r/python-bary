from distutils.core import setup
import sys

if sys.version_info[0] <= 2:
    print("python-bary requires Python3")
    sys.exit(-1)

setup(
	name 			= 'bary',
	packages		= ['bary'],
	version			= '1.0.0',
	description		= 'Barycentric coordinates computing system',
	author			= 'Delfad0r',
	author_email	= '',
	url				= 'https://github.com/Delfad0r/python-bary',
)
