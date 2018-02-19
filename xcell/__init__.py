"""Execute xCell transformation of gene expression.

This python package gives both a CLI interface and a python module to work with xCell in Python Pandas DataFrames.

Find the official R package here:

https://github.com/dviraran/xCell

And if you find this useful, please cite the authors' publication:

Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol. 2017 Nov 15;18(1):220. doi: 10.1186/s13059-017-1349-1. PubMed PMID: 29141660; PubMed Central PMCID: PMC5688663.

"""
import argparse, sys, os
import pandas as pd 
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

def xCell(expression_df,
         rnaseq=True,
         scale=True,
         alpha=0.5,
         nperm=250,
         parallel_sz=0,
         verbose=False,
         tempdir= None,
         beta_pval=False,
         perm_pval=False,
         matrix=False
         ):
    """xCell function for use with pandas DataFrame objects

    :param expression_df: REQUIRED: Expression data indexed on gene names column labels as sample ids
    :type expression_df: pandas.DataFrame
    :param parallel_sz: Number of processors to use when doing the calculations in parallel. This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel_sz=0) then it will use all available core processors unless we set this argument with a smaller number. If snow is loaded then we must set this argument to a positive integer number that specifies the number of processors to employ in the parallel calculation.
    :type parallel_sz: int Default: 0
    :param verbose: Gives information about each calculation step.
    :type verbose: bool Default: False
    :param tempdir: Location to write temporary files
    :type tempdir: string Default: System Default
    :returns: pandas.DataFrame
    """
    if matrix and (beta_pval or perm_pval): raise ValueError("can't return pvalues as a matrix")
    df = expression_df

    if not tempdir:
        tempdir =  mkdtemp(prefix="weirathe.",dir=gettempdir().rstrip('/'))
    if verbose:
        sys.stderr.write("Caching to "+tempdir+"\n")
    ## Remove genes from the genesets that do not occur in the dataset
    #members = gmt_df['member'].unique()
    #missing = set(members)-set(df.index)
    #original = df.index
    #if len(missing) > 0:
    #    if verbose: sys.stderr.write("WARNING removing "+str(len(missing))+\
    #      " genes from gene sets that don't exist in the data\n"+\
    #      ",".join(sorted(list(missing)))+"\n")
    #gmt_df = gmt_df[~gmt_df['member'].isin(list(missing))]
    ## Write our gene sets
    #gmt_df = gmt_df.groupby(['name']).\
    #    apply(lambda x: "\t".join(sorted(list(x['member'])))).reset_index().rename(columns={0:'members'})
    #of = open(os.path.join(tempdir,"gs.gmt"),'w')
    #for row in gmt_df.itertuples():
    #    name = row.name
    #    description = 'description'
    #    fields = row.members
    #    of.write(name+"\t"+description+"\t"+fields+"\n")
    #of.close()
    df.to_csv(os.path.join(tempdir,"expr.csv"))
    cur = os.path.dirname(os.path.realpath(__file__))
    rscript = os.path.join(cur,"xcell.r")
    cmd = ["Rscript",rscript]+[str(x) for x in \
        [rnaseq,scale,alpha,nperm,parallel_sz,verbose,tempdir,beta_pval,perm_pval]]
    if verbose: sys.stderr.write(" ".join(cmd)+"\n")
    sp = Popen(cmd,stdout=PIPE,stderr=PIPE)
    if not verbose: sp.communicate()
    else:
        for line in sp.stderr: sys.stderr.write(line)
    if verbose: sys.stderr.write("finished R script\n")
    output1 = pd.read_csv(os.path.join(tempdir,"pathways.csv"),index_col=0)
    output1.index.name = 'name'
    if matrix: return output1
    df = output1.unstack().reset_index().rename(columns={0:'score','level_0':'sample'})

    if beta_pval:
        output2 = pd.read_csv(os.path.join(tempdir,"beta.csv"),index_col=0)
        output2.index.name = 'name'
        output2.columns = output1.columns
        d2 = output2.unstack().reset_index().rename(columns={0:'beta_pval','level_0':'sample'})
        df = df.merge(d2,on=['sample','name'])

    if perm_pval:
        output3 = pd.read_csv(os.path.join(tempdir,"randomP.csv"),index_col=0)
        output3.index.name = 'name'
        output3.columns = output1.columns
        output4 = pd.read_csv(os.path.join(tempdir,"randomD.csv"),index_col=0)
        output4.index.name = 'name'
        d3 = output3.unstack().reset_index().rename(columns={0:'permute_pval','level_0':'sample'})
        df = df.merge(d3,on=['sample','name'])

    return df

def __cli():
    args = __do_inputs()
    # Now read in the input files for purposes of standardizing inputs
    df = None
    if args.tsv_in:
        df = pd.read_csv(args.input,sep="\t",index_col=0)
    else:
        df = pd.read_csv(args.input,index_col=0)
    result = xCell(df,
                  rnaseq=args.rnaseq,
                  scale=args.scale,
                  alpha=args.alpha,
                  nperm=args.nperm,
                  parallel_sz=args.parallel_sz,
                  verbose=args.verbose,
                  tempdir=args.tempdir,
                  beta_pval=args.beta_pval,
                  perm_pval=args.perm_pval,
                  matrix=args.matrix
                 )
    sep = ','
    use_index = False
    if args.matrix: use_index = True
    if args.tsv_out: sep = "\t"
    if args.output:
        result.to_csv(args.output,sep=sep,index=use_index)
    else:
        result.to_csv(os.path.join(args.tempdir,'final.csv'),sep=sep,index=use_index)
        with open(os.path.join(args.tempdir,'final.csv')) as inf:
            for line in inf:
                sys.stdout.write(line)

def __do_inputs():
    # Setup command line inputs
    parser=argparse.ArgumentParser(description="Execute R xCell",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group0 = parser.add_argument_group('Input options')
    group0.add_argument('input',help="Use - for STDIN")
    group0.add_argument('--tsv_in',action='store_true',help="Exepct CSV by default, this overrides to tab")


    group2 = parser.add_argument_group('Output options')
    group2.add_argument('--tsv_out',action='store_true',help="Override the default CSV and output TSV")
    group2.add_argument('--output','-o',help="Specifiy path to write transformed data")
    group2.add_argument('--meta_output',help="Speciify path to output additional run information")
    group2.add_argument('--matrix',action='store_true',help="Output results as a matrix")

    group4 = parser.add_argument_group("Add pvalue calculations. Cannot output as matrix. These will be added as columns to the DataFrame")
    group4.add_argument('--beta_pval',action='store_true',help="output the beta pvalue")
    group4.add_argument('--perm_pval',action='store_true',help="output the random permutation pvalue")

    group1 = parser.add_argument_group('command options')
    parallel_sz_str = '''
Number of processors to use when doing the calculations in parallel. This requires
to previously load either the parallel or the snow library. If parallel is
loaded and this argument is left with its default value (parallel_sz=0) then it
will use all available core processors unless we set this argument with a smaller
number. If snow is loaded then we must set this argument to a positive integer
number that specifies the number of processors to employ in the parallel
calculation.
    '''
    group1.add_argument('--parallel_sz',type=int,default=0,help=parallel_sz_str)
    verbose_str = '''
Gives information about each calculation step.
    '''
    group1.add_argument('--verbose',action='store_true',help=verbose_str)
    rnaseq_str = '''
Inputs are RNAseq.
    '''
    group1.add_argument('--rnaseq',type=bool, default=True,help=rnaseq_str)
    scale_str = '''
Scaling transforms with fit.vals.
    '''
    group1.add_argument('--scale',type=bool, default=True,help=scale_str)
    alpha_str = '''
Value to override spillover alpha parameter.
    '''
    group1.add_argument('--alpha',type=float, default=0.5,help=alpha_str)
    nperm_str = '''
Number of random resamplings.
    '''
    group1.add_argument('--nperm',type=int, default=250,help=nperm_str)

    # Temporary working directory step 1 of 3 - Definition
    label4 = parser.add_argument_group(title="Temporary folder parameters")
    group3 = label4.add_mutually_exclusive_group()
    group3.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
    group3.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")


    args = parser.parse_args()
    if args.matrix and (args.beta_pval or args.perm_pval): raise ValueError("can't return pvalues in a matrix.")
    setup_tempdir(args)
    return args  

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return

if __name__=="__main__":
  __cli()

