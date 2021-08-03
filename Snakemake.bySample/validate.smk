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

## Group column must not overlap with Name column
#if not len(set(samples.Group).intersection(set(samples.Name)))==0:
#	print( "Error: Sample Group name must not overlap with Name" )
#	sys.exit(1)

## Only alphanumeric / dash / underbar / dot in sample sheet
## Must start with alphanumeric only
invalid_elem=[]
for col in samples:
    #tmp = samples[col].str.count(r'(^[a-zA-Z0-9][a-zA-Z0-9-_\.]+$)')
    tmp = samples[col].str.count(r'(^[a-zA-Z0-9][a-zA-Z0-9-_\.]+$)')
    index_invalid = (tmp == 0)
    if index_invalid.any():
        invalid_elem = invalid_elem + samples[col][index_invalid].tolist()

if len(invalid_elem) > 0:
    print( "Error: Must be at least two character; Only alphanumeric, dash (-) and underbar (_) are allowed in a sample sheet")
    print( "Invalid values:" )
    for elem in invalid_elem: print( "  - %s" % elem )
    sys.exit(1)



## To add
## - file check
if peak_mask.upper() != "NULL":
    assert os.path.exists(peak_mask)
