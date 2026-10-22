from pathlib import Path, PurePath
import subprocess

def extract_region(bam_file, region, output_file):
    # Construct the samtools command
    command = ['samtools', 'view', '-h', bam_file, region, '-o', output_file]

    # Execute the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check if the command was successful
    if result.returncode == 0:
        print("Region extracted successfully.")
    else:
        print("Error in running samtools:", result.stderr)

def create_filename(line):
    # Split the line by ';' to isolate each attribute
    attributes = line.split(';')
    
    # Initialize empty dictionary to store key-value pairs
    attr_dict = {}
    #print(attributes)
    # Process each attribute
    for attr in attributes:
        # Strip any leading/trailing whitespace
        attr = attr.strip().replace('\"','')
        if attr:
            # Split the attribute into key and value using the first space found
            parts = attr.split(' ', 1)
            if len(parts) == 2:  # Make sure there are exactly two parts
                key, value = parts
                # Remove extra double quotes from value and store in dictionary
                attr_dict[key.strip()] = value.strip('\"')
    
    # Extract gene_name and gene_id from the dictionary
    gene_name = attr_dict.get('gene_name', 'default_gene_name')
    gene_id = attr_dict.get('gene_id', 'default_gene_id')
    
    # Create the filename
    filename = f"{gene_name}.{gene_id}.sam"
    return filename

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE
gtf_filename = 'gencode.v40.basic.annotation.CancerGeneCensus.gtf'
gtf_file_path =Path(PurePath(base_path, gtf_filename))
bam_file = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/EvRo_ethnic_files/test_bamfile/LP6005443-DNA_D06.srt.aln.bam"

# Open the GTF file and read line by line
with open(gtf_file_path, 'r') as file:
    for line in file:
        # Skip header lines if any
        if line.startswith('#'):
            continue
        
        # Split the line into fields
        fields = line.strip().split('\t')
        
        # Check if the third field is 'gene'
        if fields[2] == 'gene':
            region = f"{fields[0].lstrip('chr')}:{int(fields[3])-500}-{int(fields[4])+500}"
            output_file = Path(PurePath(base_path, "LP6005443-DNA_D06", create_filename(fields[8])))
            print((int(fields[4])+500)-(int(fields[3])-500))
            #print(output_file)
            #extract_region(bam_file, region, output_file)

