'''
Configuration & sample sheet validation
'''

################################################
## Sample information Validation

## if column numbers are the same in sample.tsv
def validate_tsv_with_pandas(file_path):
    expected_cols = None
    errors = []

    # Read file manually (so we control filtering)
    with open(file_path, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, start=1):
            stripped = line.strip()

            # Skip empty and comment lines
            if not stripped or stripped.startswith("#"):
                continue

            # Split columns
            cols = stripped.split("\t")
            col_count = len(cols)

            if expected_cols is None:
                expected_cols = col_count
            elif col_count != expected_cols:
                errors.append(
                    f"Line {line_num}: expected {expected_cols}, found {col_count}"
                )

    # Now safely load into pandas (knowing it's consistent or not)
    df = pd.read_csv(file_path, sep="\t", comment="#",skip_blank_lines=True, dtype=str)

    return df, errors


df, errors = validate_tsv_with_pandas(src_sampleInfo)
if errors:
    print("Column mismatch detected:")
    for e in errors:
        print(e)
        sys.exit(1)


## Check if required columns exist in sample.tsv
req_columns = ['Id', 'Name', 'Group', 'Fq1', 'Fq2', 'Ctrl', 'PeakMode']
missing_columns = [col for col in req_columns if col not in samples.columns]
if missing_columns:
    print(f"Error: Missing required columns: {', '.join(missing_columns)}", file=sys.stderr)
    sys.exit(1)

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
for col in ['Id', 'Fq1', 'Fq2', 'Ctrl', 'PeakMode']:
    is_valid = samples[col].str.match(r'^[a-zA-Z0-9][a-zA-Z0-9_.\-]*$')
    index_invalid = ~is_valid.fillna(True)
    if index_invalid.any():
        invalid_elem = invalid_elem + samples[col][index_invalid].tolist()

if len(invalid_elem) > 0:
    print( "Error: Must be at least two characters; Only alphanumeric, dash (-), underbar (_), dots (.) and NULL are allowed in the sample sheet. Dashes are not allowed in the Name or Group column.")
    print( "Invalid values:" )
    for elem in invalid_elem: print( "  - \"%s\"" % elem )
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

