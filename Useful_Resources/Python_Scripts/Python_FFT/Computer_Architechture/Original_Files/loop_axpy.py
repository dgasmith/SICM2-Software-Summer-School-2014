import subprocess as sp

total_flops = 5E8

output = []
for kernal in ['plain','blas']:
    print '\nStart %s kernal:' % kernal
    for k in range(4, 50, 4):
 
        k = 2**k
        iter_total_flops = int(total_flops/k)
 
        #If it is too big, just continue 
        if iter_total_flops == 0:
                continue
 
        cmd = './axpy '+str(k)+' '+str(iter_total_flops)+' '+kernal
        proc = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT, close_fds=True)
        data = proc.stdout.read()
        print data.strip()
