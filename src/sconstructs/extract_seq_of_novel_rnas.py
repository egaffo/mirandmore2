Import('*')

try:
    env = env_extract_seq_of_novel_rnas.Clone()
except NameError as ne:
    print 'mirandmore_extract_seq_of_novel_rnas: command line mode.'

    vars = Variables('vars.py')
    vars.Add('SAMPLE', '', '')
    vars.Add('MIN_COUNT', '', '')
    vars.Add('HAIRPIN_EXTENDED_TO_FOLD', '', '')
    vars.Add('MIR_TABLE', '', '')

    env = Environment(variables = vars,ENV = os.environ)
    
    Help(vars.GenerateHelpText(env))

env.Replace(MIN_COUNT = str(int(env['MIN_COUNT'])-1))

targets = [env['SAMPLE'] + "_novel_mirs_exp_seq.txt", 
           env['SAMPLE'] + "_mor_exp_seq.txt"]

novel_seqs_cmd =  '''extract_seq_of_novel_rnas.py -i ${SOURCE} -t $MIN_COUNT '''\
                  '''--hairpin_extended_to_fold $HAIRPIN_EXTENDED_TO_FOLD '''\
                  '''--mirna ${TARGETS[0]} --morna ${TARGETS[1]}'''

novel_seqs =  env.Command(targets, 
                          env['MIR_TABLE'], 
                          novel_seqs_cmd)

Return('novel_seqs')
