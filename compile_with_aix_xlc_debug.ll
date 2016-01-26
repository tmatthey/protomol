#!/bin/csh
#
# @ job_name         = compile
# @ total_tasks      = 1
# @ job_type         = parallel
# @ account_no       = parallab
#
# @ class            = express
# @ wall_clock_limit = 1:00:00
#
# @ resources        = ConsumableCpus(1) ConsumableMemory(250 mb)
#
# @ notification     = complete
# @ checkpoint       = no
# @ restart          = no
#
# @ error            = job.$(Host).$(jobid).err
# @ output           = job.$(Host).$(jobid).out
#
# @ queue


./realclean
aclocal; autoheader; autoconf; automake -a;
./configure --with-aix-xlc-debug
make depend
make clean
make



