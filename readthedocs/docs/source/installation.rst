Installation
=====

.. warning::
   SRoll need to be install with a python version <= 3.6

.. _installation:

SRoll installation package
------------------------

* To use Sroll , first install the sroll library using pip:

.. code-block:: console

   (.venv) $ pip install sroll

SRoll installation package is a PyPi library that aim to create a generic local
environement for SRoll algorithm.

* Then run the installation of SRoll using the set env library contain in the sroll package : 

.. code-block:: console

    $> python
    Python 3.6.9 (default, Nov  7 2019, 10:44:02) 
    [GCC 8.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    $>>> from sroll_package import set_env
    $>>> set_env.install()


This will clone the SRoll git repository to your current path to folder 'srollex', set python
paths, update the srollex_setenv.sh script, update Makefile and create a python virtual environement named 'py_sroll'.

Set environement
--------------

Then you need to set sroll environement by modifyng the file 'srollex_setenv.sh' in '~/srollex'.

* First add the module (if needed) for compilation with a 'module load  module1 module2 ' as follow example :

.. code-block:: console

   module load intel-cmkl-18/18.0.1.163  impi/5.1.3.258 gcc/6.3.0   
 
 
* You can also modify the hostname and global variables if needed, hereafter an example of the code needed : 
.. code-block:: bash
 
 clustername[0-9]*)
    echo "clustername detected"
    export SROLLHOST=clustername

    #load modules  for compilation
    module load intel-cmkl-18/18.0.1.163  impi/5.1.3.258 gcc/6.3.0
    
    #activate python env
    source ../py_sroll/bin/activate

    #export global variables
    export PYTHONPATH=/home1/user/py_sroll/bin/python
    export LD_LIBRARY_PATH=/home1/user/py_sroll/lib
    ;;
    
    
    
Compilation
------------
To end the installation of SRoll, it need to be compile. 

Before compilation the sroll python environement need to be activate, you can either add the following line in the srollex_setenv.sh or execute it separatly :

.. code-block:: bash
   
   source sroll_dir/py_sroll/bin/activate
  
Then set your enviromnent for sroll with :

.. code-block:: bash
   
   source sroll_dir/srollex_setenv.sh
 
Once the environment is set up compile sroll using :

.. code-block:: bash

   cd ~/sroll_dir/srollex/sroll4/
   make clean all
 
 
The compilation of SRoll will generates 3 executables troll_857 troll_cfosat and troll_14tf.

