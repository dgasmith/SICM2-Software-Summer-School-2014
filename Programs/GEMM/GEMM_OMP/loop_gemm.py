import os
import subprocess as sp

# Target number of flops per run
total_flops = 1E10

# Set thread counts to loop over
thread_list = [1, 8]

# Set kernal list to loop over
kernal_list = ['blas','block8']
# kernal_list = ['block16','block32','block48','block64','block80','block96']
# kernal_list = ['block'+str(x) for x in [16,32,64,128]]

# Set k size to loop over
k_list = range(200, 3100, 400)

output = []
for threads in thread_list:
    os.environ['OMP_NUM_THREADS'] = str(threads)

    print '\nCurrent running with %d threads.' % threads
    print '_'*80
    for kernal in kernal_list:
        print '\nStarting %s kernal:' % kernal
    
        for k in k_list:
            iter_total_flops = int(total_flops/(2*k**3))
     
            # If it is too big, just continue 
            if iter_total_flops == 0:
                    continue
            
            # Run gemm and record the timing 
            cmd = './gemm '+str(k)+' '+str(iter_total_flops)+' '+kernal
            proc = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT, close_fds=True)
            data = proc.stdout.read()
            print data.strip()
