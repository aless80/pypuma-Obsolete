from setuptools import setup, find_packages

setup(name='pypuma',
    version='0.1',
    description='bla.',
    url='https://github.com/aless80/pypuma',
    author='Alessandro Marin',
    author_email='AlessandroMarin80@gmail.com',
    license='MIT',
    packages=['pypuma'],
    install_requires=['pandas',
    'numpy',
    'networkx',
    'matplotlib'
    ],
    scripts=['bin/pypuma'],
    zip_safe=False)
