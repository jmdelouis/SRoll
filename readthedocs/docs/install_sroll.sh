
#bin/bash

{
    pip install sroll
}||{
    echo "ERROR : pip not found"
}


python -c "from sroll_package import set_env; set_env.install()"


source py_sroll/bin/activate

python -m pip install cython numpy healpy

cd srollex/
source srollex_setenv.sh
cd sroll4/
make clean all
