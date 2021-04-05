## CBB 752 Assignment 1

I'm implementing the Smith-Waterman algorithm from scratch. You can see the assignment specs [here](http://cbb752b21.gersteinlab.org/assignments). 

### System requirements

* Python 3
* git [if you don't have git try running ```pip install git```]


If you want to install and run this, run the following: 

```
git clone git://github.com/a-ram-1/cbb752-pset1.git
```

Enclosed are the following files: 

* **blosum62.txt** -- the blosum62 scoring matrix
* **input.txt** -- the input file
* **sample-input1.txt** -- a sample input file
* **sample-input2.txt** -- another sample input file
* **sample-output1.txt** -- what the output should look like if you run the algorithm with **sample-input1.txt**
* **sample-output2.txt** -- what the output should look like if you run the algorithm with **sample-input2.txt**
* **pset1.py** -- code for pset 1
* **output.txt** -- output when running my script on the input file

Here's the command you'll run: 

```
python pset1.py -i [inputfile].txt -s blosum62.txt 
```

There are two additional optional parameters, -o and -e which are gap opening and extension penalties respectively. By default they are set to -2 and -1. If you wanted to set them, add arguments like "-o -3 -e -4".

