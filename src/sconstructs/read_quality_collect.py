import os

Import('*')

try:
    env = env_read_quality_stats.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('SAMPLES', '', '')
    vars.Add('FASTQC_HTMLS', '', '')
    env = Environment(variables = vars, ENV = os.environ)
 
    Help(vars.GenerateHelpText(env))
 
    SAMPLES = env['SAMPLES']
    fastqc_html_files = env['FASTQC_HTML_FILES'].split(',')
 
copied_fastqc = []
#for source in fastqc_html_files:
for source in env['FASTQC_HTMLS']:
    for sample in env['SAMPLES']:
        if sample + os.path.sep in source.path:
            copied = env.Command(sample + '_$SOURCE.file', source, "cp $SOURCE $TARGET")
            copied_fastqc.append(copied)

Return('copied_fastqc')

