from intron_comparison import process_file

path = '/home/julia/Documents/licencjat/'

scaffolds_dict = {}
with open(path + 'longa_rascaf.fasta') as f:
    for line in f.readlines():
        if line[0] == '>':
            scaff = line[1:].strip()
        else:
            scaffolds_dict[scaff] = line.strip()


