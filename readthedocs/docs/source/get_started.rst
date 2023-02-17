Get started with SRoll
======================

.. _get_started:

Unit test
----------
First verify that SRoll is correctly installed by running sroll unit test.

In case of using a job scheduler (example for PBS) :

.. code-block:: bash

   $> cd path_to_sroll/srollex/run_troll/
   $> cd /unit_test/
   $> qsub sroll_unit_test_run.pbs

Otherwise it can be start using mpirun :

.. code-block:: bash

    $> cd path_to_sroll/srollex/sroll4/
    $> mpirun -np 4 ./troll_14tf params_unit_test.py

With X the number of process needed.

The output of Sroll unit test will write log in :

.. code-block:: bash

     path_to_sroll/SrollEx/run_troll/unit_test/sroll4_output_unit_test.log


SRoll parameters 
-----------------

* Create a python parameters file :

To use SRoll first create a python parameter file according to your project, with correct paths and values for your data.
Hereafter a list a the parameters for SRoll : :ref:`parameters`.


* Add a new parameter :

You can add a new parameter to python file but you will need to add it also to 'troll_param.h' with the correct type (ie PIOINT,PIODOUBLE,...).
   
   

SRoll run 
-----------

Once you have your parameters file ready you can run sroll. 

* In case of using a Job Scheduler software, write a script as below . Example
with PBS scheduler :

.. code-block:: bash

    #!/ bin / bash
    # PBS -q mpi
    # PBS -l select =3: ncpus =28: mem =115 g
    # PBS -l walltime =10:00:00
    # PBS -N unit_test_sroll
    # PBS -e unit_test_sroll_error.log
    # PBS -o unit_test_sroll_output.log
    # PBS -m n
   
    cd /path_to_sroll/srollex/
    module purge
    source srollex_setenv.sh
    cd sroll4
    module load impi /5.1.3.258
   
    mpirun - np X ./troll_Y parameters_file.py


* Else you can directly run mpirun as below : 

.. code-block:: bash

    $> cd path_to_sroll/srollex/sroll/
    $> mpirun - np X ./troll_Y parameters_file.py

With X the number of process needed and Y the version of sroll to run (857,cfosat,14tf)
