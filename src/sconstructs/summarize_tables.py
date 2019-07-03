
Import("*")

try:
    env = env_summarize_tables.Clone()
except NameError as ne:
    vars = Variables('vars.py')
    vars.Add('MIR_TABLES', 'comma separated list of miR tables', '')
    vars.Add('UNFILTERED_VARIANT_TABLES', 'comma separated list of unfiltered variants', '')
    vars.Add('META', '', '')

    env = Environment(variables = vars, ENV = os.environ)
 
    Help(vars.GenerateHelpText(env))

    env['MIR_TABLES'] = env['MIR_TABLES'].split(',')
    env['UNFILTERED_VARIANT_TABLES'] = env['UNFILTERED_VARIANT_TABLES'].split(',')

targets = ["raw_data.txt", #"normalized_data.txt", 
           "raw_variants.txt"]#, "normalized_variants.txt"]
           
sources = [env['MIR_TABLES'], 
	   env['UNFILTERED_VARIANT_TABLES']]

mir_table_num = len(env['MIR_TABLES'])
unfiltered_variants_tables_num = len(env['UNFILTERED_VARIANT_TABLES'])

cmd = "vacuum_cleaner.R -m $META -i \"${SOURCES[0: " + str(mir_table_num) +\
      "]}\" -d $TARGET.dir -u \"${SOURCES[" + str(mir_table_num) +\
      ":" + str(unfiltered_variants_tables_num + mir_table_num + 1) + "]}\""

collective_tables = env.Command(targets, 
                                sources, 
                                cmd)

Return('collective_tables')
