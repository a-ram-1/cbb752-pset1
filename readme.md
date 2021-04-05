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


### Example input 

Let's contextualize this with an example. Say we use these input sequences: 

```
ASDASDFLWE
ALSKERASDASDLWERI
```

I'll put this in a file called example-input.txt. Now if I run the following: 

```
python pset1.py -i example-input.txt -s blosum62.txt
```

I'll get the following output: 

```
-----------
|Sequences|
-----------
sequence1
ASDASDFLWE
sequence2
ALSKERASDASDLWERI
--------------
|Score Matrix|
--------------
		A	S	D	A	S	D	F	L	W	E	
	0	0	0	0	0	0	0	0	0	0	0	
A	0	4	2	1	4	2	1	0	0	0	0	
L	0	2	2	0	2	2	0	1	4	2	1	
S	0	1	6	4	3	6	4	3	2	1	2	
K	0	0	4	5	3	4	5	3	2	1	2	
E	0	0	3	6	4	3	6	4	3	2	6	
R	0	0	2	4	5	3	4	3	2	0	4	
A	0	4	2	3	8	6	5	4	3	2	3	
S	0	2	8	6	6	12	10	9	8	7	6	
D	0	1	6	14	12	11	18	16	15	14	13	
A	0	4	5	12	18	16	16	16	15	13	13	
S	0	2	8	11	16	22	20	19	18	17	16	
D	0	1	6	14	15	20	28	26	25	24	23	
L	0	0	5	12	14	19	26	28	30	28	27	
W	0	0	4	11	13	18	25	27	28	41	39	
E	0	0	3	10	12	17	24	25	27	39	46	
R	0	0	2	9	11	16	23	24	26	38	44	
I	0	0	1	8	10	15	22	23	26	37	43	
----------------------
|Best Local Alignment|
----------------------
Alignment Score:46
Alignment Results:
      (ASDASDFLWE)  
       |||||| |||   
ALSKER(ASDASD-LWE)RI

```


