#Read input file
#   parameters:#       input file - .i input file
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
def read_input_file(input_file):
    inputs = {}
    with open(input_file) as f:
        for line in f:  # Each line in file
            if '=' in line and '#' not in line:  # Set file
                inputs[line.split("=")[0].strip()] = line.split("=")[1].strip()
            else: pass
        # Update numbers as floats
        for key,val in inputs.items():
            if is_number(val) == True:
                inputs[key] = float(val)
            if key == 'sigmas':
                sigmas = [int(sig) for sig in val.split(',')]
                inputs[key] = sigmas
    return inputs
