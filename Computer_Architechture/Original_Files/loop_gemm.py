import subprocess as sp

#Target number of flops per run
total_flops = 5E9

output = []
for kernal in ['blas']:
#for kernal in ['plain','block16','blas']:
    print '\nStarting %s kernal:' % kernal

    for k in range(4, 15, 2):
 
        k = 2**k
        iter_total_flops = int(total_flops/k**3)
 
        #If it is too big, just continue 
        if iter_total_flops == 0:
                continue
 
        cmd = './gemm '+str(k)+' '+str(iter_total_flops)+' '+kernal
        proc = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT, close_fds=True)
        data = proc.stdout.read()
        print data.strip()
