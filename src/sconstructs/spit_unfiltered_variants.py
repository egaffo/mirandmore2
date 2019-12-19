import os 

Import('*')

# outputs a table of variants not filtered

try:
    env = env_spit_unfiltered_variants.Clone()
except NameError as ne:
    print('mirandmore_spit_unfiltered_variants: command line mode.')
    
    vars = Variables('vars.py')
    vars.Add('SAMPLE', '', '')
    vars.Add('PRE_BLOB', '', '')

    env = Environment(variables = vars, ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    
target = env['SAMPLE'] + "_unfiltered_variants.txt"

unfiltered_variants_cmd = '''filtered_to_variants_table.py -i ${SOURCE} -o ${TARGET}'''

unfiltered_variants =  env.Command(target, 
                                   env['PRE_BLOB'], 
                                   unfiltered_variants_cmd)

Return('unfiltered_variants')
