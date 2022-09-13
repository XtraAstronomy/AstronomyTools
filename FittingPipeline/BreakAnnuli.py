"""
Python routine to break up .reg file for  multiple annuli into individual .reg files for all annuli
"""

# ----- Inputs ----- # 
num_annuli = 10  # Number of annuli
ann_file = '/mnt/carterrhea/carterrhea/Perseus/Test/annuli.reg'  # Full path to annuli files
output_dir = '/mnt/carterrhea/carterrhea/Perseus/Test/regions/'
# ---- Code ---- #
def break_annuli(ann_file, num_annuli, output_dir):
    """
    Function to break a single .reg file containing multiple annuli into individual .reg files for all annuli

    Args:
        ann_file: Full path to annuli file
        num_annuli: Number of annuli
        output_dir: Full path to output directory
    """
    x_coord = None; y_coord = None; inner = []; outer = [] 
    preamble = ''  # Preamble for .reg
    with open(ann_file, 'r') as file:
        #print(file)
        #preamble = file.readlines() # Get preamble
        lines = file.readlines()  # skip preamble
        ct = 0
        for line in lines: 
            if ct < 3:
                preamble += line
                ct += 1
            else:
                line = line.rstrip().replace('annulus(','').replace(')','').split(',')
                x_coord = line[0]
                y_coord = line[1]
                inner = [val for val in line[2:-1]]
                outer = [val for val in line[3:]]
    # Write to new annulus files
    for new_ann in range(num_annuli):
        with open(output_dir+'annulus_%i.reg'%new_ann, 'w+') as new_file:
            new_file.write(preamble)
            new_file.write('annulus(%s,%s,%s,%s)'%(x_coord, y_coord, inner[new_ann], outer[new_ann]))

break_annuli(ann_file, num_annuli, output_dir)