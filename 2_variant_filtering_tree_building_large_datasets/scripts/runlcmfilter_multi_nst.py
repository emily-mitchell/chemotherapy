Last login: Fri Mar  8 08:46:15 on ttys000

The default interactive shell is now zsh.
To update your account to use zsh, please run `chsh -s /bin/zsh`.
For more details, please visit https://support.apple.com/kb/HT208050.
mib115090i:~ em16$ ssh farm5-login
ssh: Could not resolve hostname farm5-login: nodename nor servname provided, or not known
mib115090i:~ em16$ ssh farm5-login
em16@farm5-login's password: 
Welcome to Ubuntu 18.04 (bionic) (4.15.0-219-generic x86_64)

* Managed by WSI Ansible.

* For information on the farm 22.04 migration please visit:
* https://ssg-confluence.internal.sanger.ac.uk/display/FARM/Ubuntu+22.04

* Auto-update: Next check 23-03-2024; Auto-reboot: Disabled

  System information as of Thu Mar 21 11:54:15 GMT 2024

  System load:  21.52               Processes:            4262
  Usage of /:   28.9% of 137.20GB   Users logged in:      78
  Memory usage: 31%                 IP address for bond0: 10.160.12.32
  Swap usage:   3%

  => There are 44 zombie processes.
Your Hardware Enablement Stack (HWE) is supported until April 2023.
*** System restart required ***
*** Livepatch has fixed kernel vulnerabilities. System restart recommended on the closest maintenance window ***





Plan your installation, and FAI installs your plan.

Last login: Thu Mar  7 09:30:10 2024 from 10.80.240.51
Loading setup for ubuntu:bionic
em16@farm5-head2:~$ cd /lustre/scratch119/realdata/mdt1/team154/ms44/scripts/
-bash: cd: /lustre/scratch119/realdata/mdt1/team154/ms44/scripts/: No such file or directory
em16@farm5-head2:~$ cd /lustre/scratch126/casm/team273jn/share/
em16@farm5-head2:/lustre/scratch126/casm/team273jn/share$ ls
jobs.py  pileups  __pycache__  README.txt  runlcmfilter_multi_nst.py  runlcmfilter.sh  test
em16@farm5-head2:/lustre/scratch126/casm/team273jn/share$ nano runlcmfilter_multi_nst.py

  GNU nano 2.9.3                            runlcmfilter_multi_nst.py                                       

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import jobs
import argparse
    

LCMFILTER_SCRIPT="/lustre/scratch126/casm/team273jn/share/runlcmfilter.sh"
NST="/nfs/cancer_ref01/nst_links/live"

def runfilters(samples,inputdir,outputdir,build,ft,vcfext,bamext):
    cmds=[]
    os.makedirs(outputdir,exist_ok=True)
    for sample in samples:
        cmds=cmds+[f"{LCMFILTER_SCRIPT} {inputdir}/{sample}/{sample}.{vcfext} {inputdir}/{sample}/{sample}.$
    lsf_dir=f"{outputdir}/lsf"
    os.makedirs(lsf_dir,exist_ok=False)
    jobs.submit_as_job_array(cmds,samples[0],lsf_dir,queue="normal",mem_gb=8,max_chunks=20)
  
if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Runs hairpin filters")
    parser.add_argument("--project",help="File containing commands",required=True)
    parser.add_argument("--samples",type=str,required=True,help="Comma delimited list of samples")
    parser.add_argument("--outdir",type=str,required=True,help="Output directory")
    parser.add_argument("--build",type=str,required=True,help="Genome build: hg19,hg38 or mm10")
    parser.add_argument("--fragment_threshold",type=int,default=4)
    parser.add_argument("--vcf_suffix",type=str,default="caveman_c.annot.vcf.gz")
    parser.add_argument("--bam_suffix",type=str,default="sample.dupmarked.bam")
    args=parser.parse_args()
    samplelist=args.samples
    samples=[x for x in samplelist.strip().split(",")]
    inputdir=f"{NST}/{args.project}"
    runfilters(samples,inputdir,args.outdir,args.build,args.fragment_threshold,args.vcf_suffix,args.bam_suf$
                             [ File 'runlcmfilter_multi_nst.py' is unwritable ]
^G Get Help    ^O Write Out   ^W Where Is    ^K Cut Text    ^J Justify     ^C Cur Pos     M-U Undo
^X Exit        ^R Read File   ^\ Replace     ^U Uncut Text  ^T To Linter   ^_ Go To Line  M-E Redo
