Import("*")

try:
    env = env_summarize_depth.Clone()
except NameError as ne:
    vars = Variables('vars.py')
    vars.Add('TRIMMING_REPORT_FILES', '', '')

    env = Environment(variables = vars, ENV = os.environ)
 
    Help(vars.GenerateHelpText(env))

    env['TRIMMING_REPORT_FILES'] = env['TRIMMING_REPORT_FILES'].split(',')


depth_csv = env.Command('depth.csv', 
			env['TRIMMING_REPORT_FILES'], 
			"collect_depth.py -i \"$SOURCES\" -o $TARGET")
 
Return('depth_csv')

