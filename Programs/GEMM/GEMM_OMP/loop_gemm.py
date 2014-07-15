import os
import subprocess as sp

#Target number of flops per run
total_flops = 1E9


#Set thread counts to loop over
thread_list = [1, 4, 8]

#Set kernal list to loop over
kernal_list = ['blas','block32','block64','block128']
#kernal_list = ['plain','block16','blas']
#kernal_list = ['block'+str(x) for x in [32, 64, 128]]

#Set k size to loop over
k_list = range(100, 1100, 200)

output = []
for threads in thread_list:
    os.environ['OMP_NUM_THREADS'] = str(threads)

    print '\nCurrent running with %d threads.' % threads
    print '_'*80
    for kernal in kernal_list:
        print '\nStarting %s kernal:' % kernal
    
        for k in k_list:

            iter_total_flops = int(total_flops/k**3)
     
            # If it is too big, just continue 
            if iter_total_flops == 0:
                    continue
            
            # Run gemm and record the timing 
            cmd = './gemm '+str(k)+' '+str(iter_total_flops)+' '+kernal
            proc = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT, close_fds=True)
            data = proc.stdout.read()
            print data.strip()
