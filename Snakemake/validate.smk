'''
Configuration & sample sheet validation
'''

################################################
## Sample information Validation

## Id / Name column must be unique
if not samples.Id.is_unique:
	print( "Error: Id column in sample.tsv is not unique" )
	sys.exit(1)
if not samples.Name.is_unique:
	print( "Error: Name column in sample.tsv is not unique" )
	sys.exit(1)

## Name and Group columns should only contain alphanumeric, underscore, dot (no dash)
invalid_elem = []
for col in ['Name', 'Group']:
    is_valid = samples[col].str.match(r'^[a-zA-Z0-9][a-zA-Z0-9_.]*$')
    index_invalid = ~is_valid.fillna(True)
    if index_invalid.any():
        invalid_elem = invalid_elem + samples[col][index_invalid].tolist()

## all other columns should only contain alphanumeric, dash, underscore, dot
for col in samples.columns:
    if col not in ['Name', 'Group']:
        is_valid = samples[col].str.match(r'^[a-zA-Z0-9][a-zA-Z0-9_.\-]*$')
        index_invalid = ~is_valid.fillna(True)
        if index_invalid.any():
            invalid_elem = invalid_elem + samples[col][index_invalid].tolist()

if len(invalid_elem) > 0:
    print( "Error: Must be at least two characters; Only alphanumeric, dash (-), underbar (_) and dots (.) are allowed in the sample sheet. Dashes are not allowed in the Name or Group column.")
    print( "Invalid values:" )
    for elem in invalid_elem: print( "  - %s" % elem )
    sys.exit(1)

###########################################
## File / directory check
## "NULL" or must exist
filelist = [ chrom_size, star_index, src_sampleInfo, cluster_yml ]
for f in filelist:
    if f.upper() != "NULL":
        if not os.path.exists(f):
            print( "Error: %s does not exist" % f )
            sys.exit(1)


dirList = [fastqDir, trimDir, dedupDir, qcDir, sampleDir]
dirNames = ["fastqDir", "trimDir", "filteredDir", "dedupDir", "qcDir", "sampleDir"]

for name, value in zip(dirNames, dirList):
    if not value:
        print(f"{name} is undefined or empty.")
        sys.exit(1)

import os
if os.path.exists("config.yml"):
    if pooling not in ["replicate", "pooled"]:
        print(f"html_report in the config file can only be defined as replicate or pooled.")
        sys.exit(1)
    