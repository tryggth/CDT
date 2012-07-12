from subprocess import Popen, PIPE
import time
import numpy

k0 = 1.0
k3 = 0.78

f = open('tuning.data','w')
f.write("Data from tuning: k0=%d, k3=%d\n" % (k0,k3))
f.close()

for muprime in range(500,1000,1):
    mu = muprime / 1000.0
    f = open('tuning.script','w')
    f.write("""(load "../utilities.lisp")
(load "globals.lisp")
(load "simplex.lisp")
(load "moves.lisp")
(load "initialization.lisp")
(load "montecarlo.lisp")

(initialize-t-slices-with-v-volume :num-time-slices 64
				   :target-volume(* 8 1024)
				   :spatial-topology"S2"
				   :boundary-conditions"PERIODIC")

(set-k0-k3-alpha-lambda-mu %f %f -1.0 1.0 %f)
(setf SAVE-EVERY-N-SWEEPS 100)
(setf NUM-SWEEPS 50000)

(generate-data-console)""" % (k0,k3,mu))
    f.close()
    print mu
    proc = Popen(["/usr/local/sbcl-1.0.49/bin/sbcl","--dynamic-space-size","1024","--script","tuning.script"], stdout=PIPE)
    time.sleep(15)
    proc.terminate()
    output = proc.stdout.read()
    print output
    output = output.split('\n')
    output = filter(lambda x: x[0:5] == "start", output)
    std = numpy.std(map(lambda x: int(x.split(' ')[14][:-1]), output))
    if output:
        f = open('tuning.data','a')
        f.write("mu: %f, std: %f" % (mu, std))
    else:
        f = open('tuning.data','a')
        f.write("mu: %f, std: error" % (mu))
    f.close()
                

