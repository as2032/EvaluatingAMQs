#Evaluating different AMQ's

Alexander Straub
Hw3

##Input
The amqs.py script requires two text file's as input. One text file represents K, a set of k-mers of any type, separated by '\n', which will be the input set that will be used to initialize the datastructures. The second text file represents K', another set of k-mers of any type separated by '\n', which will be the query keys to query the datastructures built on K.
An example input file can be found in ```input_files```.

##Running the code

Once all requirements are installed (those in requirements.txt) running the program requires various command line arguments.

Arg1 (string):
filename for key set K
Arg2 (string):
filename for key set K'
Arg3 (float):
error_rate (between 0 and 1) for bloom filter
Arg4 (int):
Number of bits to store for fingerprint array, one of (7,8,9)
Arg5 (int):
Number of true positives in query key set K'
Arg6 (int):
Number of true negatives in query key set K'

Example:
```python3 amqs.py 15k_31_kmerT.txt 15k_31_kmerFP_t4000f2000.txt 0.0078125 7 4000 2000```
