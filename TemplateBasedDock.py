import os,re,sys,datetime,itertools
import traceback
import pandas as pd
from subprocess import check_output
import InterPred as IP
import csvLoader


def run_tm_align_parse_output(template, target, target_chain="_"):
    # Runs TM align between the target and a template, parses the output and puts it in a dictionary

    try:

        TMC = os.environ["TM"]

    except KeyError:
        print("Environment variable not found. Please make sure that the InterPred.run script has been correctly setup")


    tm_align_output = check_output([TMC, '-B', target, '-chB', target_chain, '-A', template[0], '-chA', template[1]])


    pattern=re.compile((
    '(.*\.[pe][dn][bt])\schain:\s(\w).*\n'
    '(.*\.[pe][dn][bt])\schain:\s(\w).*\n' #\_(\w)
    'Length of Chain_1:\s+(\d+).+\n'
    'Length of Chain_2:\s+(\d+).+\n\n'
    'Aligned length=\s*(\d+), RMSD=\s*([\d+\.]+), Seq_ID=n_identical\/n_aligned=\s*([\d+\.]+).*\n'
    'TM-score=\s*([\d+\.]+).*\n'
    'TM-score=\s*([\d+\.]+).*\n\n'
    '([A-Z-]+)\n'
    '^(.+)$\n'
    '([A-Z-]+)\n'
    ), re.MULTILINE)

    match = pattern.findall(tm_align_output)
    keys = ['pdb1', 'ch1', 'pdb2', 'ch2', 'len1', 'len2', 'ali_len', 'rmsd', 'seqid', 'tm1', 'tm2', 'ali1', 'dot_colon_str', 'ali2']
    for values in match:

        parsed_line = dict(zip(keys, values))
        return parsed_line
class TemplateBasedDock(IP.InterPred):
    PDBrepo = '/home/limin/limin/InterPred/interpred/homology_models'
    TEMPLATErepo = '/home/limin/limin/InterPred/task0110/out'
    def __init__(self):
        self.context = os.getcwd()
        super(TemplateBasedDock, self)
    def make_model(self, pdb1, pdb2, pdbt ):
# pdb1,2 in forms of absolute pdb path
# pdbt in forms of absolute pdb path and 2 selected chain code
        self._env_set()
        import make_coarse_models
        import make_structural_alignments
        def structural_alignment(template, target, chain='_'):
            id = os.path.splitext(os.path.basename(target))[0] + '_' + chain
            out_file = os.environ['TMPDIR'] + '/' + id + '.structural_alignments' 
            hit =run_tm_align_parse_output(
                template,
                target, chain)
            if not hit:
                raise Exception('hit return None')
            hit = [hit]
            df=pd.DataFrame(hit)
            df.to_csv(out_file)

            return out_file
            
        #structural alignment
        #this method will return a absolute path of the *.structural_alignment
        sa1 = structural_alignment([pdbt[0], pdbt[1]], pdb1) 
        sa2 = structural_alignment([pdbt[0], pdbt[2]], pdb2) 
        #assert that pdb1 and pdb2 has only one molecule and the chain id is
        #empty
        outfiles = make_coarse_models.main((sa1, pdb1, '_', sa2, pdb2, '_'))
        self._env_restore()
        return sa1, sa2, outfiles
    def grabPDBpath(self,target):
        files = os.listdir(self.PDBrepo)
        files = filter(
            lambda f: re.match('{}_\d.pdb'.format(target),f),
            files
        )
        assert(files)
        return map(lambda f: os.path.join(self.PDBrepo, f),files)
    def grabTEMPLATEpath(self,target):
        assert(os.path.isfile(target))
        return os.path.abspath(target)
    def errlog(self,msg ,filename='./errlog', stdout=False):
        f = open(self.context + '/' + filename, 'a')
        f.write('#'*10 + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        + '#'*10 + '\n')
        f.write(msg + '\n')
        f.close()
        if stdout:
           sys.stderr.write(msg+'\n' )
    def run(self, config):
        def preprocess(s):
           return os.path.splitext(os.path.basename(s))[0]
        f = open(config)
        cl = csvLoader._csv_r(f)
        ret = []
        for line in cl.reader:
            try:
                pdb1 = self.grabPDBpath(preprocess(line['SegA'])) #
                pdb2 = self.grabPDBpath(preprocess(line['SegB'])) #
                pdbt =[
                self.grabTEMPLATEpath(line['path']),line['codeA'],line['codeB']] #
            except Exception as e:
                _, _, tb = sys.exc_info()
                filename, line, func, text = traceback.extract_tb(tb)[-1]
                print('An error occurred on line {} in statement {}'.format(line, text))
                self.errlog(str(e), stdout=True)
                traceback.print_tb(tb)
                continue
            for (p1, p2) in itertools.product(pdb1, pdb2):
                print '{} and {} align to {}'.format(p1,p2,pdbt[0])
                try:
                    ret.append(self.make_model(p1, p2, pdbt))
                except Exception as e:
                    _, _, tb = sys.exc_info()
                    filename, line, func, text = traceback.extract_tb(tb)[-1]
                    print('An error occurred on line {} in statement {}'.format(line, text))
                    traceback.print_tb(tb)
                    self._env_restore()
                    raise
        return ret
if __name__ == '__main__':
    import pprint,sys,json
    so = TemplateBasedDock()
#    res = so.make_model(
#        pdb1 = '/home/limin/limin/InterPred/interpred/homology_models/CSF1R_HUMAN_seg0_0.pdb',
#        pdb2 = '/home/limin/limin/InterPred/interpred/homology_models/CSF1_HUMAN_seg0_0.pdb',
#        pdbt = [
#            '/home/limin/limin/InterPred/task0110/out/4WRL_20180104174819744000.pdb',
#            'A','B']
#    )
#    pp = pprint.PrettyPrinter(indent=4)
#    pp.pprint(res)
    ret = so.run(sys.argv[1])
    json.dump(ret, open('log.json','w'), indent=4)
