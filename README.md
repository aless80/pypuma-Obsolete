## TODO
s1,s2,t1,t2 are hardcoded  
rename scripts from panda to puma and rebuild the setup
rename all 'panda' to 'puma' int his README and in all scripts
rewrite this README. Anaconda nice to mention

The main script is still called pyp but if you after the install see s1,s2, t1,t2 in the command options for this:  
pypanda --help  
you can make sure this is pypuma


## PyPanda (Python Panda)
Python implementation of PANDA (Passing Attributes between Networks for Data Assimilation)  

_Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions, PLoS One, 2013 May 31;8(5):e64832_

### Table of Contents
* [Panda implementation](#panda-algorithm)  
* [Installation](#installation)  
* [Usage](#usage)  
  * [python](#run-from-python)
  * [Terminal](#run-from-the-terminal)  
* [Results] (#results)

### Panda algorithm
To find agreement between the three input networks first the responsibility (R) is calculated.  

<img src="https://github.com/aless80/pypuma/raw/develop/img/responsibility.png" height="30">  

Thereafter availability (A) is calculated.  

<img src="https://github.com/aless80/pypuma/raw/develop/img/availability.png" height="30">  

Availability and responsibility are combined with the following formula.  

<img src="https://github.com/aless80/pypuma/raw/develop/img/combine.png" height="30">  

Protein cooperativity and gene co-regulatory networks are updated.  

<img src="https://github.com/aless80/pypuma/raw/develop/img/cooperativity.png" height="30">  
<img src="https://github.com/aless80/pypuma/raw/develop/img/co-regulatory.png" height="30">  

P and C are updated to satisfy convergence.  

<img src="https://github.com/aless80/pypuma/raw/develop/img/p.png" height="30">  
<img src="https://github.com/aless80/pypuma/raw/develop/img/c.png" height="30">  

Hamming distance is calculated every iteration.  

<img src="https://github.com/aless80/pypuma/raw/develop/img/hamming.png" height="40">  


### Installation
PyPuma requires Python 2.7. We recommand the following commands to install PyPuma (on Ubuntu and Debian derived systems, also works on OSX):
#### Using a virtual environment
Using [python virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/) is the cleanest installation method. 

Cloning git and setting up the [python virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/):
```no-highlight
pip install --user pipenv   #Make sure you have Pipenv
git clone https://github.com/aless80/pypuma.git
cd pypuma
virtualenv pypumaenv #Create a folder for the virtual environment inside the cloned git folder 
source pypumaenv/bin/activate
```
Installing and running pypuma:
```no-highlight
(pypumaenv)$ pip install -r requirements.txt
(pypumaenv)$ python setup.py install #--user
```

Complete uninstall:
```no-highlight
(pypuma)$ deactivate	#Exit
rm -rf ....
```

#### Using pip on the user's install directory
```no-highlight
git clone https://github.com/aless80/pypuma.git
cd pypuma
python setup.py install --user
#to run from the command line you will need to make pypuma executable and add the bin directory to your PATH:
cd bin
chmod +x pypuma
echo "$(pwd):PATH" >> ~/.bashrc
source ~/.bashrc
```
To run PyPuma from Windows (tested on Windows 10) install Git (https://git-scm.com/downloads) and Anaconda Python2.7 (https://www.continuum.io/downloads) and from the Anaconda Prompt run:
```no-highlight
git clone https://github.com/aless80/pypuma.git
cd pypuma
python setup.py install
```
### Usage
#### Run from the terminal
PyPuma can be run directly from the terminal with the following options:
```
-h help
-e (required) expression values
-m (optional) pair file of motif edges, when not provided analysis continues with Pearson correlation matrix
-p (optional) pair file of PPI edges
-f (optional) remove missing values (default is Fales)
-o (required) output file
-q (optional) output lioness single sample network
```
To run PyPuma on the example data:
```
$ pypuma -e ToyData/ToyExpressionData.txt -m ToyData/ToyMotifData.txt -p ToyData/ToyPPIData.txt -f True -o test_panda.txt -q test_lioness.txt
```
To reconstruct a single sample Lioness Pearson correlation network:
```
$ pypuma -e ToyData/ToyExpressionData.txt -o test_panda_pearson.txt -q test_lioness_pearson.txt
```
#### Run from python
Fire up your python shell or ipython notebook. 
Import PyPuma library:
```python
from pypuma import Puma
from pypuma import Lioness
import pandas as pd
```
Run Panda algorithm, leave out motif and PPI data to use Pearson correlation network:
```python
p = Panda('ToyData/ToyExpressionData.txt', 'ToyData/ToyMotifData.txt', 'ToyData/ToyPPIData.txt', remove_missing=False)
```
Save the results:
```python
p.save_panda_results(file = 'Toy_Panda.pairs')
```
Return a network plot:
```python
plot = AnalyzePanda(p)
plot.top_network_plot(top=100, file='top_100_genes.png')
```
Calculate indegrees for further analysis:
```python
indegree = p.return_panda_indegree()
```
Calculate outdegrees for further analysis:
```python
outdegree = p.return_panda_outdegree()
```
Run the Lioness algorithm for single sample networks:
```python
l = Lioness(p)
```
Save Lioness results:
```python
l.save_lioness_results(file = 'Toy_Lioness.txt')
```
Return a network plot for one of the Lioness single sample networks:
```python
plot = AnalyzeLioness(l)
plot.top_network_plot(column= 0, top=100, file='top_100_genes.png')
```
### Results
```
Example Panda output:
TF  Gene  Motif Force
---------------------
CEBPA	AACSL	0.0	-0.951416589143
CREB1	AACSL	0.0	-0.904241609324
DDIT3	AACSL	0.0	-0.956471642313
E2F1	AACSL	1.0	3.6853160511
EGR1	AACSL	0.0	-0.695698519643

Example lioness output:
Sample1 Sample2 Sample3 Sample4
-------------------------------
-0.667452814003	-1.70433776179	-0.158129613892	-0.655795512803
-0.843366539284	-0.733709815256	-0.84849895139	-0.915217389738
3.23445386464	2.68888472802	3.35809757371	3.05297381396
2.39500370135	1.84608635425	2.80179804094	2.67540878165
-0.117475863987	0.494923925853	0.0518448588965	-0.0584810456421

TF, Gene and Motif order is identical to the panda output file.
```
