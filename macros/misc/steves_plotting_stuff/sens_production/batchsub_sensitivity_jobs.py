#!/usr/bin/env python

# Submits a bunch of sensitivity jobs

import math
import os
import subprocess
import distutils.spawn
import sys

try:
  os.environ["T2KOSCANA"]
except:
  print "Environment not set up, exiting"
  exit(-1)


queue="long" 
walltime="24:00"
#~ mem_limit=4194304
mem_limit=524288
vmem_limit=8388608

n_submitted = 0

script_location = os.path.dirname(os.path.realpath(__file__))

try:
  working_dir=os.environ["PBS_O_WORKDIR"]
except KeyError:
  working_dir=os.environ["PWD"]

output_dir = ""
for possible_outdir in ["OUTFOLDER","PWD"]:
  try:
    output_dir = os.path.abspath( os.environ[possible_outdir] )
    break
  except:
    continue

if not output_dir: output_dir = working_dir
  
prog    =  os.path.abspath( script_location+"/run_sens_job.py" )
log_dir =  os.path.abspath( output_dir+"/log/" )

if not os.path.isdir(output_dir):
  print "Making output folder ",output_dir
  os.makedirs(output_dir)
if not os.path.isdir(log_dir):
  print "Making log folder ",log_dir
  os.makedirs(log_dir)

extra_qsub = ["-o", log_dir, "-e", log_dir]

def main(include1d=True,include2d=False,
         min_seed = 1, max_seed=10,
         fit_param = None, bins_per_job=21,
         hierarchy="normal", tag="valor_sens",
         config="analysis_configuration.xml",
         dataset_id = 0):
  # Warwick cluster
  if distutils.spawn.find_executable("bsub"):
    base_command="bsub"
    base_args=[base_command,"-M",mem_limit,"-v",vmem_limit, "-q", queue]
    job_name_switch="-J"
  
  # torgue w/ bqsub (ral cluster)
  elif distutils.spawn.find_executable("bqsub"):
    base_command="bqsub"
    base_args=[base_command,"-c",walltime]+extra_qsub
    job_name_switch="-N"
  
  # torque
  elif distutils.spawn.find_executable("qsub"):
    base_command="qsub"
    resources = "walltime="+`walltime`+",mem="+`mem_limit`+"kb"
    base_args=[base_command,"-l",resources]+extra_qsub
    job_name_switch="-N"
  
  else:
    print >>sys.stderr,"No supported queue submission found"
    exit(-1)
  
  var_bins = dict( [ ("mixing23", 40),
                     ("dmsq"  , 20)
                  ])
  
  seeds = range(min_seed,max_seed+1)
  if not min_seed:
    import random
    seeds = [random.randint(1,2**31) for s in seeds]
  
  for seed in seeds:
    #1d
    if include1d:
      for surface in var_bins.iterkeys():
        var1=surface
        nbins1=var_bins[var1]
        njobs= int(math.ceil(float(nbins1)/bins_per_job))
        
        for ijob in xrange(njobs):
          min_bin1 = bins_per_job*ijob+1
          max_bin1 = min( nbins1, min_bin1+bins_per_job )
          
          job_info = [tag,seed,surface,min_bin1,max_bin1]
          job_name = "_".join(map(str,job_info))
          
          submit_args = base_args+[job_name_switch, job_name]
          my_args   = [prog,
                       "-r",seed,
                       "--p1",surface,
                       "--p2","none",
                       "-t",tag,
                       "-m",hierarchy,
                       "-c",config,
                       "--min-param1-bin", min_bin1,
                       "--max-param1-bin", max_bin1,
                       "-d",dataset_id,
                      ]
          
          if fit_param is not None: my_args+=["-f", fit_param]
          
          submit(submit_args,my_args)
          
    #2d
    if include2d:
      for surface in ["mixing23_dmsq"]:
        var1,var2=surface.split("_")
        nbins1=var_bins[var1]
        nbins2=var_bins[var2]
        njobs= int(math.ceil(float(nbins2)/bins_per_job))
        
        for bin1 in xrange(nbins1):
          min_bin1 = bin1
          max_bin1 = bin1
          
          for ijob in xrange(njobs):
            min_bin2 = bins_per_job*ijob
            max_bin2 = min( nbins2, min_bin2+bins_per_job ) - 1

            job_info = [tag,seed,surface,min_bin1,max_bin1,min_bin2,max_bin2]
            job_name = "_".join(map(str,job_info))
            
            submit_args = base_args+[job_name_switch, job_name]
            my_args   = [prog,
                         "-r",seed,
                         "--p1",var1,
                         "--p2",var2,
                         "-t",tag,
                         "-m",hierarchy,
                         "-c",config,
                         "--min-param1_bin", min_bin1,
                         "--max-param1_bin", max_bin1,
                         "--min-param2_bin", min_bin2,
                         "--max-param2_bin", max_bin2,
                         "-d",dataset_id,
                        ]
            
            if fit_param is not None: my_args+=["-f", fit_param]
            
            submit(submit_args,my_args)
            
  print "Submitted ",n_submitted," jobs"
  
def submit(submit_args,my_args):
  base_command=submit_args[0]
  
  if base_command=="bsub":
    full_call = map(str,submit_args+my_args)
    print " ".join(full_call)
    subprocess.call(full_call)
    
  elif base_command=="bqsub":
    full_call = map(str,submit_args+my_args)
    print " ".join(full_call)
    subprocess.call(full_call)
    
  elif base_command=="qsub":
    my_call=" ".join(map(str,my_args))
    
    sub_call = " ".join(map(str,submit_args))
    
    full_call = "echo "+my_call+ " | "+sub_call
    print full_call
    subprocess.call(full_call,shell=True)
  
  global n_submitted
  n_submitted+=1

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Submit many VALOR sensitivity jobs")
  parser.add_option('-1','--with-1d', action="store_true", dest='submit1d',default=False,help='Submit jobs for 1d plots')
  parser.add_option('-2','--with-2d', action="store_true", dest='submit2d',default=False,help='Submit jobs for 2d plots')
  
  parser.add_option('--min-seed'  , dest='min_seed'   , default=1, help="Min random seed, 0 sets all random" , type='int')
  parser.add_option('--max-seed'  , dest='max_seed'   , default=3, help="Max random seed" , type='int')

  parser.add_option('-f','--fit-param' , dest='fit_param'  , default=None, help="Fit param list")
  
  parser.add_option('-m', dest="hierarchy", default="normal", help='Mass hierarchy: "normal" or "inverted"')
  
  parser.add_option('-b', '--bins-per-job', dest='bins_per_job', default=21, help="Number of fits to do per job",type='int')
  parser.add_option('-t', '--tag', dest='tag', default="valor_sens", help="Tag for output files")
  
  parser.add_option('-c', '--config', dest='config', default="analysis_configuration.xml", help="Configuration file")
  
  parser.add_option('-d', '--dataset', dest='dataset', default=0, help="Which dataset in the XML file to use", type='int')
  
  options,other = parser.parse_args()
  
  if not (options.submit1d or options.submit2d):
    print "No plots to generate selected"
    sys.exit(-1)
  else:
    include1d = options.submit1d
    include2d = options.submit2d
  
  if other:
    print "Unrecognised arguments: "+other
    sys.exit(-1)
    
  sys.argv=[]

  main(include1d, include2d,
       options.min_seed, options.max_seed,
       options.fit_param, options.bins_per_job,
       options.hierarchy, options.tag,
       options.config, options.dataset)
