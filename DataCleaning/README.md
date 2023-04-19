This code will clean your Chandra observations -- it has been tailored for galaxy clusters. 

To run this code, you need to first create an input file (see MS0735.i or Perseus.i for the structure).

Once you have your input file, you can enter the terminal, activate your ciao installation, and run the following
`python DataCleaningPipeline.py InputFileName.i` where InputFileName is the name of the input file you just made.

This will create a folder called `repro` in the ObsID folder.