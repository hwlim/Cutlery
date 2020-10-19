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

## Only alphanumeric / dash / underbar in sample sheet
## NOTE: should be upgraded to regular expression of allowed pattern (not only checking a space)
## Series.str.count(r'(^[a-zA-Z0-9][a-zA-Z0-9-_]+$)').sum()

invalid_elem=[]
for col in samples:
    flag = samples[col].str.contains(" ")
    if flag.any:
        invalid_elem = invalid_elem + samples[col][ flag.tolist() ].tolist()

if len(invalid_elem) > 0:
    print( "Error: Only alphanumeric, dash (-) and underbar (_) are allowed in a sample sheet" )
    sys.exit(1)

## To add
## - file check
